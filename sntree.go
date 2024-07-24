package main

import (
	// "fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/gotree/tree"
)

func SNTree(aln align.Alignment) *tree.Tree {
	// filter out sites with more than two characters
	// create tree, then add bipartitions per unique, non-conflicting site.
	aln.Sort()
	taxaNames := make([]string, len(aln.Sequences()))
	for i, seq := range aln.Sequences() {
		taxaNames[i] = seq.Name()
	}
	snSplits := snSplits(aln)
	// fmt.Printf("sn-splits %v\n", snSplits)
	// PrintSplits(snSplits)
	snTree, err := BuildTree(snSplits, taxaNames)
	if err != nil {
		panic(err)
	}
	// fmt.Println(snTree)
	return snTree
}

func snSplits(aln align.Alignment) []*Split {
	splits, err := CreateSplits(aln, nil)
	// PrintSplits(splits)
	if err != nil { // shouldn't happen
		panic(err)
	}
	conflicts := make(map[int]bool)
	// fmt.Printf("%d\n", aln.Length())
	for i := 0; i < len(splits); i++ {
		// if !conflicts[i] {
		for j := i + 1; j < len(splits); j++ {
			// fmt.Println("pair", i, j)
			if b, err := splits[i].Conflict(splits[j]); err != nil {
				panic(err)
			} else if b {
				conflicts[i] = true
				conflicts[j] = true
				// break
			}
		}
		// }
	}
	result := []*Split{}
	for i := range len(splits) {
		if !conflicts[i] {
			// fmt.Println(i)
			result = append(result, splits[i])
		}
	}
	return result
}

// func snSplits(aln align.Alignment) []int {
// 	splits, err := CreateSplits(aln, nil)
// 	if err != nil { // shouldn't happen
// 		panic(err)
// 	}
// 	conflicts := make(map[int]bool)
// 	fmt.Printf("%d\n", aln.Length())
// 	for i := 0; i < aln.Length(); i++ {
// 		if !conflicts[i] {
// 			for j := i + 1; j < aln.Length(); j++ {
// 				if conflict(i, j, aln) {
// 					conflicts[i] = true
// 					conflicts[j] = true
// 					break
// 				}
// 			}
// 		}
// 	}
// 	result := []int{}
// 	for i := range aln.Length() {
// 		if !conflicts[i] {
// 			result = append(result, i)
// 		}
// 	}
// 	return result
// }

// // could probably be optimized further
// func conflict(i, j int, aln align.Alignment) bool {
// 	sitePair, err := aln.SelectSites([]int{i, j})
// 	if err != nil { // should never happen
// 		panic(err)
// 	}
// 	comb := [4]bool{} // check whether four combinations of 00 01 10 11 are present
// 	conflict := false
// 	for _, pair := range sitePair.Sequences() {
// 		// TODO: make sure I've checked that we're using the correct character set
// 		comb[((pair.CharAt(0)-'0')*2)+(pair.CharAt(1)-'0')] = true
// 		conflict = comb[0] && comb[1] && comb[2] && comb[3]
// 		if conflict {
// 			break
// 		}
// 	}
// 	return conflict
// }
