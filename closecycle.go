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
		if v[2] > bestScore || bestScore == -1 {
			bestTree, bestScore = k, v[2]
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

func CloseCycle(bestTree *tree.Tree, taxa []string, aln align.Alignment, n int) *tree.Tree {
	// if !slices.IsSorted(taxa) { // make sure the bitset order matches between tree alignment
	// 	panic("my assumption that taxa are sorted is wrong")
	// }
	// fmt.Println("tree")
	// for _, n := range bestTree.SortedTips() {
	// 	fmt.Println(n.Name())
	// }
	slices.Sort(taxa)
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

	backbone := findBackbone(bestTree, edgeScores, postorderPass, preorderPass)
	fmt.Println(backbone)
	attachTaxa(bestTree, backbone, taxa[x], createLabel(n))
	return bestTree
}

func preprocessEdgeScores(bestTree *tree.Tree, splits []*Split, x int) [][2]int {
	scores := make([][2]int, len(bestTree.Edges()))
	fmt.Print("splits")
	PrintSplits(splits)
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
	fmt.Println("ogbitset", ogBitset.String())
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
	fmt.Println("left", xLeft.String())
	fmt.Println("right", xRight.String())
	fmt.Println("x", x)
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
		if root == cur { // root has no corresponding edge
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
		p, err := prev.Parent()
		if err != nil && err.Error() == "The node has more than one parent" {
			panic(err)
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
			} else if edge.Left() == p || edge.Right() == p {
				result[e.Id()] += result[edge.Id()]
			} else {
				result[e.Id()] += postorderScores[edge.Id()]
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

func findBackbone(bestTree *tree.Tree, edgeScores [][2]int, post, pre []int) [2]int {
	scoresPost := make([]int, len(edgeScores)) // score if backbone ends at the index
	postStarts := make([]int, len(edgeScores)) // edge ids for where backbone starts if it ends at the index
	root := bestTree.Root()
	bestTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur == root { // no edge
			return true
		}
		score := edgeScores[e.Id()][0] + edgeScores[e.Id()][1]
		if cur.Tip() {
			scoresPost[e.Id()] = score
			postStarts[e.Id()] = e.Id()
			return true
		}
		children := childEdges(cur)
		left := score + scoresPost[children[0].Id()] + post[children[1].Id()]
		right := score + scoresPost[children[1].Id()] + post[children[0].Id()]
		startAtCur := score + post[children[1].Id()] + post[children[0].Id()]
		if startAtCur >= left && startAtCur >= right {
			scoresPost[e.Id()] = score
			postStarts[e.Id()] = e.Id()
		} else if left >= right {
			scoresPost[e.Id()] = left
			postStarts[e.Id()] = postStarts[children[0].Id()]
		} else {
			scoresPost[e.Id()] = right
			postStarts[e.Id()] = postStarts[children[1].Id()]
		}
		return true
	})
	fmt.Println("post scores", scoresPost)
	scoresPre := make([]int, len(edgeScores))
	preStarts := make([]int, len(edgeScores))
	bestTree.PreOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur == root { // no edge
			return true
		}
		score := edgeScores[e.Id()][0] + edgeScores[e.Id()][1]
		edges := preorderEdges(cur, root)
		if prev == root {
			// fmt.Printf("HERE!! %d\n", root.Nneigh())
			left := score + scoresPost[edges[0].Id()] + post[edges[1].Id()]
			right := score + scoresPost[edges[1].Id()] + post[edges[0].Id()]
			fmt.Printf("root: left %d, right %d\n", left, right)
			if left >= right {
				scoresPre[e.Id()] = left
				preStarts[e.Id()] = postStarts[edges[0].Id()]
			} else {
				scoresPre[e.Id()] = right
				preStarts[e.Id()] = postStarts[edges[1].Id()]
			}
		} else {
			p, err := prev.Parent()
			if err != nil && err.Error() == "The node has more than one parent" {
				panic(err)
			}
			if edges[0].Left() == p || edges[0].Right() == p {
				score += scoresPre[edges[0].Id()] + post[edges[1].Id()]
				preStarts[e.Id()] = preStarts[edges[0].Id()]
			} else {
				score += scoresPre[edges[1].Id()] + post[edges[0].Id()]
				preStarts[e.Id()] = preStarts[edges[1].Id()]
			}
			scoresPre[e.Id()] = score
		}
		return true
	})
	fmt.Println("pre scores", scoresPre)

	// find maximum over all possible backbones
	maxScore := 0            // best possible score for a backbone
	maxStart, maxEnd := 0, 0 // ids of best start and end point edges
	for i, s := range scoresPost {
		totalScore := s + pre[i] // need to add on subtree to end of path to get real total score
		if totalScore >= maxScore {
			maxScore = totalScore
			maxStart, maxEnd = postStarts[i], i
		}
	}
	for i, s := range scoresPre {
		totalScore := s + post[i]
		if totalScore >= maxScore {
			maxScore = totalScore
			maxStart, maxEnd = preStarts[i], i
		}
	}
	// I need to fix the bug in the preorder preprocess bit
	return [2]int{maxStart, maxEnd}
}

func childEdges(node *tree.Node) []*tree.Edge {
	children := children(node)
	result := make([]*tree.Edge, 0)
	for _, c := range children {
		for _, e := range c.Edges() {
			if e.Left() == node || e.Right() == node {
				result = append(result, e)
			}
		}
	}
	if len(result) != 2 {
		panic(fmt.Sprintf("did not find two edges, found %d", len(result)))
	}
	return result
}

func attachTaxa(bestTree *tree.Tree, backbone [2]int, taxaX, label string) {
	// bifurcate the first edge, adding labeled internal vertex aattachment point
	// then do the same for the second, but add the attachment point as a leaf
	edges := bestTree.Edges()
	x := bestTree.NewNode()
	x.SetName(taxaX)
	bestTree.GraftTipOnEdge(x, edges[backbone[0]])
	// I cannot figure out how to simply add an interanl vertex to a bestTree
	// (it might not even be possible), so I'm just adding a comment, then I'll edit the newick string

	// code to add network edge labels
	// e := x.Edges()[0]
	// e.AddComment(label)
	// if backbone[0] != edges[backbone[0]].Id() {
	// 	panic("wrong edge")
	// }
	// tip := bestTree.NewNode()
	// tip.SetName(label)
	// bestTree.GraftTipOnEdge(tip, edges[backbone[1]])
}

func createLabel(i int) string {
	return fmt.Sprintf("%s#R%d", intToAlphabet(uint(i)), i)

}

func intToAlphabet(n uint) string {
	result := ""
	for n > 0 {
		remainder := n % 26
		result = string('A'+remainder) + result
		n = n / 26
	}
	return result
}
