package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/gotree/io/nexus"
	"github.com/evolbioinfo/gotree/tree"
	"github.com/fredericlemoine/bitset"
)

func ReadTaxa(dir string) map[int][]string {
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	taxa := make(map[int][]string, 0)
	var id int
	for _, file := range files {
		if nMatch, err := fmt.Sscanf(file.Name(), "taxa_%d.txt", &id); err != nil {
			fmt.Println(file.Name(), nMatch, id)
			// panic(err)
			continue
		} else if nMatch != 1 {
			panic("expected to match 1 value")
		} else {
			content, err := os.ReadFile(filepath.Join(dir, file.Name()))
			if err != nil {
				panic(err)
			}
			taxa[id] = strings.Split(string(content), "\n")
		}
	}
	return taxa
}

func ReadPAUPResults(dir string, nPolytomy uint) []*tree.Tree {
	trees := make([]*tree.Tree, nPolytomy)
	for i := range nPolytomy {
		trees[i] = readPolytomy(dir, i)
	}
	return trees
}

func readPolytomy(dir string, i uint) *tree.Tree {
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	var j int
	candidates := make(map[int][]int) // possible trees (each with a different taxa removed)
	for _, file := range files {
		if nMatch, err := fmt.Sscanf(file.Name(), "polytomy_"+strconv.Itoa(int(i))+"_%d_scores.tsv", &j); err != nil {
			// panic(err) // TODO: catch error if this file doesn't exist as any of the files
			continue
		} else if nMatch != 1 {
			panic("expected to match 2 values")
		} else {
			fmt.Println(file.Name())
			filePointer, err := os.Open(filepath.Join(dir, file.Name()))
			if err != nil {
				panic(err)
			}
			defer filePointer.Close()
			reader := csv.NewReader(filePointer)
			reader.Comma = '\t'
			records, err := reader.ReadAll()
			if err != nil {
				panic(err)
			}
			// This loop selects which of the multiple trees *on the same taxa* outputted by PAUP* is best
			maxTree, maxVal := -1, -1
			for k, record := range records {
				if k != 0 { // skip header
					curVal, err := strconv.Atoi(record[1])
					if err != nil {
						panic(err)
					}
					if curVal < maxVal || maxVal == -1 {
						maxTree, maxVal = k-1, curVal // -1 as we don't count the heade
					}
				}
			}
			candidates[j] = []int{maxTree, maxVal}
		}
	}
	fmt.Println(candidates)
	bestTree, bestScore := -1, -1
	for k, v := range candidates {
		if v[1] < bestScore || bestScore == -1 {
			bestTree, bestScore = k, v[1]
		}
	}
	fmt.Printf("best score for polytomy %d is %d\n", i, bestScore)
	// read tree
	treeFile, err := os.Open(fmt.Sprintf("%s/polytomy_%d_%d_trees.nex", dir, i, bestTree))
	if err != nil {
		panic(err)
	}
	nxs, err := nexus.NewParser(treeFile).Parse()
	if err != nil {
		panic(err)
	}
	fmt.Println(nxs)
	var result *tree.Tree
	tIndex := 0
	nxs.IterateTrees(func(s string, t *tree.Tree) { // I don't think there's a better way to get the nth tree with this library
		if tIndex == candidates[bestTree][0] {
			result = t
		}
		tIndex++
	})
	return result
}

func getSubalignment(aln align.Alignment, taxa []string) align.Alignment {
	subaln := align.NewAlign(align.UNKNOWN)
	for _, t := range taxa {
		seq, exists := aln.SequenceByName(t)
		if exists {
			subaln.AddSequence(t, seq.Sequence(), "")
		} else {
			panic(fmt.Sprintf("sequence for taxa %s does not exist", t))
		}
	}
	return subaln
}

func CloseCycle(bestTree *tree.Tree, taxa []string, aln align.Alignment) []*tree.Tree {
	// if !slices.IsSorted(taxa) { // make sure the bitset order matches between tree alignment
	// 	panic("my assumption that taxa are sorted is wrong")
	// }
	// fmt.Println("tree")
	// for _, n := range bestTree.SortedTips() {
	// 	fmt.Println(n.Name())
	// }
	x := 0
	for i, t := range taxa {
		if !slices.Contains(bestTree.AllTipNames(), t) {
			x = i
		}
	}
	bestTree.ReinitInternalIndexes()
	// bestTree.UpdateTipIndex()
	// bestTree.UpdateBitSet()
	subaln := getSubalignment(aln, taxa)
	splits, err := CreateSplits(subaln, nil)
	if err != nil { // shouldn't happen
		panic(err)
	}
	edgeScores := preprocessEdgeScores(bestTree, splits, x)
	fmt.Println("edge scores", edgeScores)
	// postorder and preorder pass scores (assuming edges are not part of backbone)
	postorderPass := postorderScore(bestTree, edgeScores)
	preorderPass := preorderScore(bestTree, edgeScores, postorderPass)
	fmt.Println(postorderPass)
	fmt.Println(preorderPass)

	return []*tree.Tree{}
}

func preprocessEdgeScores(bestTree *tree.Tree, splits []*Split, x int) [][2]int {
	scores := make([][2]int, len(bestTree.Edges()))
	bestTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur != bestTree.Root() {
			scores[e.Id()] = scoreEdge(e, splits, x)
		}
		return true
	})
	return scores
}

func scoreEdge(e *tree.Edge, splits []*Split, x int) [2]int {
	// build bipartitions
	// if e.Left().Tip() || e.Right().Tip() {
	// 	return [2]int{0, 0}
	// }
	ogBitset := e.Bitset()
	// fmt.Println(ogBitset.String())
	if ogBitset == nil {
		panic("bitset is nil")
	}
	xLeft, xRight := bitset.New(ogBitset.Len()+1), bitset.New(ogBitset.Len()+1)
	for i := range ogBitset.Len() {
		if i < uint(x) && ogBitset.Test(i) {
			xLeft.Set(uint(i))
			xRight.Set(uint(i))
		} else if i >= uint(x) && ogBitset.Test(i) {
			xLeft.Set(uint(i + 1))
			xRight.Set(uint(i + 1))
		}
	}
	xRight.Set(uint(x))
	// fmt.Println("left", xLeft.String())
	// fmt.Println("right", xRight.String())
	// fmt.Println("x", x)
	return [2]int{CountMatches(splits, &Split{split: xLeft}), CountMatches(splits, &Split{split: xRight})}
}

func postorderScore(bestTree *tree.Tree, scores [][2]int) []int {
	result := make([]int, len(bestTree.Edges())+1)
	bestTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if bestTree.Root() == cur { // root has no corresponding edge
			return true
		}
		if (e.Left() != cur && e.Right() != cur) || (e.Left() != prev && e.Right() != prev) {
			panic("cur, prev and the edge e don't work the way I think")
		}
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if p != prev {
			panic("prev is not parent of cur")
		}
		if e.Left() == cur {
			result[e.Id()] = scores[e.Id()][1]
		} else {
			result[e.Id()] = scores[e.Id()][0]
		}
		if !cur.Tip() {
			for _, c := range children(cur) {
				result[e.Id()] += result[c.Id()]
			}
		}
		return true
	})
	return result
}

func children(node *tree.Node) []*tree.Node {
	p, err := node.Parent()
	if err != nil {
		panic(err)
	}
	children := make([]*tree.Node, 0)
	for _, n := range node.Neigh() {
		if n != p {
			children = append(children, n)
		}
	}
	if len(children) != 2 {
		panic("tree is not binary")
	}
	return children
}

func preorderScore(bestTree *tree.Tree, scores [][2]int, postorderScores []int) []int {
	result := make([]int, len(bestTree.Edges())+1)
	root := bestTree.Root()
	bestTree.PreOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if bestTree.Root() == cur { // root has no corresponding edge
			return true
		}
		if (e.Left() != cur && e.Right() != cur) || (e.Left() != prev && e.Right() != prev) {
			panic("cur, prev and the edge e don't work the way I think")
		}
		if p, err := cur.Parent(); err != nil {
			panic(err)
		} else if p != prev {
			panic("prev is not parent of cur")
		}
		if e.Left() == cur {
			result[e.Id()] = scores[e.Id()][0]
		} else {
			result[e.Id()] = scores[e.Id()][1]
		}
		edges := preorderEdges(cur, root)
		for _, edge := range edges {
			if prev == root {
				result[e.Id()] += postorderScores[edge.Id()]
			} else {
				result[e.Id()] += result[edge.Id()]
			}
		}
		return true
	})
	return result
}

// func preorderEdges(node, root *tree.Node) []*tree.Edge {
// 	result := make([]*tree.Edge, 0)
// 	p, err := node.Parent()
// 	if err != nil {
// 		panic(err)
// 	}
// 	for _, edge := range p.Edges() {
// 		if (edge.Left() != node && edge.Right() != node) || (edge.Left() != root && edge.Right() != root) {
// 			result = append(result, edge)
// 		} else if edge.Left() != root && edge.Right() != root {
// 			for _, redge := range root.Edges() {
// 				if redge != edge {
// 					result = append(result, redge)
// 				}
// 			}
// 		}
// 	}
// 	if len(result) != 2 {
// 		panic("preorder edge finder did not return two edges")
// 	}
// 	return result
// }

func preorderEdges(node, root *tree.Node) []*tree.Edge {
	result := make([]*tree.Edge, 0)
	p, err := node.Parent()
	if err != nil {
		panic(err)
	}
	for _, edge := range p.Edges() {
		if edge.Left() != node && edge.Right() != node {
			result = append(result, edge)
		}
	}
	if len(result) != 2 {
		panic(fmt.Errorf("preorder edge finder did not return the expected number of edges. Returned %d edges; root %t", len(result), p == root))
	}
	return result
}
