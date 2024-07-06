package main

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/gotree/tree"
	"github.com/fredericlemoine/bitset"
)

type Split struct {
	split *bitset.BitSet
}

// Creates list of splits from alignment.
// Selects all sites if sites is nil.
func CreateSplits(aln align.Alignment, sites []int) ([]*Split, error) {
	aln.Sort()
	if sites != nil {
		var err error
		if aln, err = aln.SelectSites(sites); err != nil {
			return nil, fmt.Errorf("cannot create splits: %w", err)
		}
	}
	splits := make([]*Split, aln.Length())
	for i := range aln.Length() {
		splits[i] = &Split{split: bitset.New(uint(len(aln.Sequences())))}
	}
	for row, seq := range aln.Sequences() {
		fmt.Println(seq.Name())
		for column := range aln.Length() {
			if seq.CharAt(column) == '1' {
				splits[column].split.Set(uint(row))
			}
		}
	}
	return splits, nil
}

func (s *Split) Length() int {
	// return len(s.split)
	return int(s.split.Len())
}

func PrintSplits(ss []*Split) {
	for _, s := range ss {
		fmt.Println(s.split.String())
	}
}

// could probably be optimized further
func (s1 *Split) Conflict(s2 *Split) (bool, error) {
	if s1.Length() != s2.Length() {
		return false, fmt.Errorf("split lengths %d and %d do not match", s1.Length(), s2.Length())
	}
	bInt := map[bool]int8{false: 0, true: 1}
	comb := [4]bool{} // check whether four combinations of 00 01 10 11 are present
	conflict := false
	for i := range s1.Length() {
		// TODO: make sure I've checked that we're using the correct character set
		comb[bInt[s1.split.Test(uint(i))]*2+bInt[s2.split.Test(uint(i))]] = true
		conflict = comb[0] && comb[1] && comb[2] && comb[3]
		if conflict {
			break
		}
	}
	return conflict, nil
}

func (s *Split) Clade(taxa []string) []string {
	clade := make([]string, 0)
	for i := range s.split.Len() {
		if s.split.Test(i) {
			clade = append(clade, taxa[i])
		}
	}
	return clade
}

// TODO: optimize this later
func CountMatches(ss []*Split, split *Split) int {
	if ss[0].Length() != split.Length() {
		panic("split lengths not equal")
	}
	// fmt.Println("split", split.split.String())
	count := 0
	for _, s := range ss {
		// fmt.Println("ss", s.split.String())
		if s.split.Equal(split.split) {
			// if s.split.EqualOrComplement(split.split) {
			count++
		}
	}
	return count
}

func BuildTree(splits []*Split, taxa []string) (*tree.Tree, error) {
	starTree, err := tree.StarTree(len(taxa))
	if err != nil { // only happens if there is less than two taxa, which shouldn't happen
		panic(err)
	}
	i := 0
	starTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur.Tip() {
			cur.SetName(taxa[i])
			i++
		}
		return true
	})
	if len(taxa) != i {
		panic("there should be as many leaves in the star tree as taxa")
	}
	for _, s := range splits {
		fmt.Println(s.Clade(taxa))
		fmt.Println(taxa)
		clade := s.Clade(taxa)
		if len(clade) == 0 || len(clade) == 1 || len(clade) == len(taxa) || len(clade) == len(taxa)-1 {
			break // not a bipartition if all taxa are on one side, also don't include trivial bipartitions (we already have them)
		}
		node, edges, monophyletic, err := starTree.LeastCommonAncestorUnrooted(nil, s.Clade(taxa)...)
		if err != nil {
			return nil, fmt.Errorf("error building tree: %w", err)
		} else if !monophyletic {
			return nil, fmt.Errorf("splits are not compatible")
		}
		starTree.AddBipartition(node, edges, 1.0, 1.0)
	}
	return starTree, nil
}
