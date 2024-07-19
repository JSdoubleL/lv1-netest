package main

import (
	"fmt"

	"github.com/evolbioinfo/gotree/tree"
)

func AssembleNetwork(sntree *tree.Tree, cycles []*tree.Tree) *tree.Tree {
	// var cur *tree.Tree
	fmt.Println(cycles)
	newSplits := make([][]string, 0)
	for _, c := range cycles {
		c.ReinitIndexes()
		splits := SplitsFromTree(c)
		PrintSplits(splits)
		for _, s := range splits {
			fmt.Println("split", s.Clade(c.AllTipNames()))
			newSplits = append(newSplits, s.Clade(c.AllTipNames()))
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
