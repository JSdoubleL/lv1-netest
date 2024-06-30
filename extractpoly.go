package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/nexus"
	"github.com/evolbioinfo/gotree/tree"
)

func ExtractPolytomies(snTree *tree.Tree) [][]string {
	snTree.ReinitInternalIndexes()
	// err := snTree.UpdateBitSet()
	// if err != nil {
	// 	panic(err)
	// }
	poly := make([][]string, 0) // list of polytomies, represented as list of their taxa labels
	taxaNames := snTree.SortedTips()
	snTree.PostOrder(func(cur, prev *tree.Node, e *tree.Edge) (keep bool) {
		if cur.Nneigh() > 3 {
			if len(cur.Edges()) != cur.Nneigh() {
				panic("assert fail") // remove if this doesn't trip
			}
			poly = append(poly, make([]string, 0))
			for _, edge := range cur.Edges() {
				split := edge.Bitset()
				if !split.Any() || split.All() {
					panic("edge isn't a split") // really don't think this should happen
				}
				for i := range split.Len() {
					if edge.Left() == cur && split.Test(i) {
						poly[len(poly)-1] = append(poly[len(poly)-1], taxaNames[i].Name())
						break
					} else if edge.Right() == cur && !split.Test(i) {
						poly[len(poly)-1] = append(poly[len(poly)-1], taxaNames[i].Name())
						break
					}
				}
			}
		}
		return true
	})
	return poly
}

func WritePolytomies(polytomies [][]string, aln align.Alignment) {
	os.Mkdir("output", 0755) // TODO: let user set prefix
	for i, polytomy := range polytomies {
		err := os.WriteFile(fmt.Sprintf("output/taxa_%d.txt", i), []byte(strings.Join(polytomy, " ")), 0644)
		if err != nil {
			panic(fmt.Errorf("could not write file: %w", err))
		}
		for j := range len(polytomy) {
			// fmt.Println(polytomy[:j], polytomy[j+1:])
			subsetTaxa := make([]string, len(polytomy)-1)
			for k := range len(polytomy) - 1 {
				if k < j {
					subsetTaxa[k] = polytomy[k]
				} else {
					subsetTaxa[k] = polytomy[k+1]
				}
			}
			// subsetTaxa := append(polytomy[:j], polytomy[j+1:]...)
			fmt.Println(subsetTaxa)
			out := align.NewAlign(align.UNKNOWN)
			for _, seqName := range subsetTaxa {
				seq, exists := aln.GetSequenceByName(seqName)
				if exists {
					out.AddSequence(seqName, seq.Sequence(), "")
				} else {
					panic("sequence does not exist")
				}
			}
			out.Alphabet()
			nexusStr := nexus.WriteAlignment(out)
			err = os.WriteFile(fmt.Sprintf("output/polytomy_%d_%d.nex", i, j), []byte(fixNexus(nexusStr, i, j)), 0644)
			if err != nil {
				panic(fmt.Errorf("could not write file: %w", err))
			}
		}
	}
}

func fixNexus(nexusStr string, i, j int) string {
	nexusStr = strings.Replace(nexusStr, "dmension", "dimension", -1) // there's a spelling error for some reason
	nexusStr = strings.Replace(nexusStr, "format datatype=dna;", "format datatype = standard gap = - missing = ? symbols = \" 0 1\";", -1)
	paupBlock := fmt.Sprintf(`
begin assumptions;
	options deftype=unord;
 end;
 
begin paup;
	set criterion=parsimony maxtrees=100 increase=no;
	hsearch start=stepwise addseq=random nreps=25 swap=tbr collapse=no;
	filter best=yes;
	describetrees 1/diag=yes;
	pscores all/ ci ri rc hi scorefile=polytomy_%d_%d_scores.tsv replace=yes;
	savetrees file=polytomy_%d_%d_trees.nex replace=yes format=nexus;
	quit;
end;`, i, j, i, j)
	return nexusStr + paupBlock
}
