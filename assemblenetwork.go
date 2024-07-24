package main

import (
	"fmt"
	"slices"

	"github.com/evolbioinfo/gotree/tree"
)

func AssembleNetwork(sntree *tree.Tree, cycles []*tree.Tree) *tree.Tree {
	sntree.ReinitIndexes()
	// var cur *tree.Tree
	// fmt.Println(cycles)
	newSplits := make([][]string, 0)
	taxaNames := sntree.AllTipNames()
	slices.Sort(taxaNames)
	nameToID := make(map[string]int)
	for i, t := range taxaNames {
		nameToID[t] = i
	}
	for i, c := range cycles {
		c.ReinitIndexes()
		splits := SplitsFromTree(c)
		// PrintSplits(splits)
		for _, s := range splits {
			// fmt.Println("split", s.Clade(c.AllTipNames()))
			// fmt.Println("expand", expandSplits(sntree, nameToID, s.Clade(c.AllTipNames()), i))
			newSplits = append(newSplits, expandSplits(sntree, nameToID, s.Clade(c.AllTipNames()), i))
		}
	}
	for _, clade := range newSplits {
		poly, edges, _, err := sntree.LeastCommonAncestorUnrooted(nil, clade...)
		if err != nil {
			panic(err)
		}
		// TODO check that poly is correct vertex
		sntree.AddBipartition(poly, edges, 1, 1)
	}
	return sntree
}

func expandSplits(sntree *tree.Tree, nameToID map[string]int, ogSplit []string, n int) []string {
	//get polytomy
	//map each taxon to an edge
	//expand split to include new taxa - consider root edge case
	search, err := sntree.SelectNodes(fmt.Sprintf("polytomy_%d", n))
	// fmt.Println(fmt.Sprintf("polytomy_%d", n))
	if err != nil {
		panic(err)
	}
	if len(search) != 1 {
		// fmt.Println(sntree.Newick())
		// fmt.Println(search)
		panic("multiple (or none) nodes match polytomy name")
	}
	poly := search[0]
	result := make([]string, 0)
	for _, taxon := range ogSplit {
		subtaxa := make([]string, 0)
		var parentEdge *tree.Edge
		for _, e := range poly.Edges() {
			if e.Left() == poly && e.Bitset().Test(uint(nameToID[taxon])) {
				subtaxa = append(subtaxa, GetClade(*e.Bitset(), sntree.AllTipNames(), false)...)
				// } else if e.Bitset().Test(uint(nameToID[taxon])) {
				// 	result = append(result, GetClade(*e.Bitset(), sntree.AllTipNames(), false)...)
			}
			if e.Right() == poly {
				parentEdge = e
			}
		}
		if len(subtaxa) == 0 {
			subtaxa = append(subtaxa, GetClade(*parentEdge.Bitset(), sntree.AllTipNames(), true)...)
		}
		result = append(result, subtaxa...)
	}
	return result
}

// func getSplits(sntree, cycle *tree.Tree, n int) []string {
// 	// for each edge in the polytomy in the sntree - get the taxa names corresponding each branch
// 	// get all the splits from the cycle tree, then exapnd them to the full taxa set
// 	taxa := sntree.Tips()
// 	search, err := sntree.SelectNodes(fmt.Sprintf("polytomy_%d", n))
// 	fmt.Println(fmt.Sprintf("polytomy_%d", n))
// 	if err != nil {
// 		panic(err)
// 	}
// 	if len(search) != 1 {
// 		fmt.Println(search)
// 		panic("multiple nodes match polytomy name")
// 	}
// 	poly := search[0]
// 	for _, e := range poly.Edges() {
// 		bs := e.Bitset()
// 		for t := range
// 	}
// 	p, err := poly.Parent()
// 	if err != nil && poly != sntree.Root() {
// 		panic(err)
// 	}
// 	result := make(map[string]*tree.Tree)
// 	for _, t := range cycle.Tips() {
// 		for _, e := range poly.Edges() {
// 			if e.Left() == p || e.Right() == p {
// 				remainder := sntree.Clone()
// 				remainder.RemoveTips(true, sntree.SubTree(poly).AllTipNames()...)
// 				result[t.Name()] = remainder
// 			} else {
// 				if e.Right() != poly {
// 					result[t.Name()] = sntree.SubTree(e.Right())
// 				} else {
// 					result[t.Name()] = sntree.SubTree(e.Left())
// 				}
// 			}
// 		}
// 	}
// 	return result
// }

// func AssembleNetwork(sntree *tree.Tree, cycles []*tree.Tree) *tree.Tree {
// 	var cur *tree.Tree
// 	fmt.Println(cycles)
// 	for i, c := range cycles {
// 		cur = c
// 		subtrees := getSubtrees(sntree, c, i)
// 		fmt.Println(c.Newick())
// 		for _, t := range cur.Tips() {
// 			fmt.Print(i)
// 			cur.GraftTreeOnTip(t.Name(), subtrees[t.Name()])
// 		}
// 	}
// 	return cur
// }

// func getSubtrees(sntree, cycle *tree.Tree, n int) map[string]*tree.Tree {
// 	search, err := sntree.SelectNodes(fmt.Sprintf("polytomy_%d", n))
// 	fmt.Println(fmt.Sprintf("polytomy_%d", n))
// 	if err != nil {
// 		panic(err)
// 	}
// 	if len(search) != 1 {
// 		fmt.Println(search)
// 		panic("multiple nodes match polytomy name")
// 	}
// 	poly := search[0]
// 	p, err := poly.Parent()
// 	if err != nil && poly != sntree.Root() {
// 		panic(err)
// 	}
// 	result := make(map[string]*tree.Tree)
// 	for _, t := range cycle.Tips() {
// 		for _, e := range poly.Edges() {
// 			if e.Left() == p || e.Right() == p {
// 				remainder := sntree.Clone()
// 				remainder.RemoveTips(true, sntree.SubTree(poly).AllTipNames()...)
// 				result[t.Name()] = remainder
// 			} else {
// 				if e.Right() != poly {
// 					result[t.Name()] = sntree.SubTree(e.Right())
// 				} else {
// 					result[t.Name()] = sntree.SubTree(e.Left())
// 				}
// 			}
// 		}
// 	}
// 	return result
// }
