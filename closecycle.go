package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/gotree/io/nexus"
	"github.com/evolbioinfo/gotree/tree"
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
		if nMatch, err := fmt.Sscanf(file.Name(), "polytomy_%d_%d_scores.tsv", &i, &j); err != nil {
			// panic(err) // TODO: catch error if this file doesn't exist as any of the files
			continue
		} else if nMatch != 2 {
			panic("expected to match 2 values")
		} else {
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
						maxTree, maxVal = k-1, curVal // -1 as we don't count the header
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

func CloseCycles(bestTree *tree.Tree, taxa string, aln align.Alignment) []*tree.Tree {
	return []*tree.Tree{}
}
