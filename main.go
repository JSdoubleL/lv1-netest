package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/nexus"
	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

type args struct {
	alignmentFile string
	polytomyDir   string
	setup         bool
}

func main() {
	args := parseArgs()
	aln, err := readAlignment(args.alignmentFile)
	if err != nil {
		panic(err)
	}
	if args.setup {
		sntree := SNTree(*aln)
		fmt.Println(sntree)
		fmt.Println("SN-Tree generated...")
		polytomies := ExtractPolytomies(sntree)
		fmt.Printf("%d polytomies extracted...\n", len(polytomies))
		WritePolytomies(polytomies, *aln, args.polytomyDir)
		WriteTree(fmt.Sprintf("%s/sntree.nwk", args.polytomyDir), sntree, false)
		fmt.Println("done.")
	} else {
		taxa := ReadTaxa(args.polytomyDir)
		sntree := ReadTree(fmt.Sprintf("%s/sntree.nwk", args.polytomyDir))
		fmt.Println(taxa)
		bestTrees := ReadPAUPResults(args.polytomyDir, uint(len(taxa)))
		fmt.Println(taxa, bestTrees)
		cycles := make([]*tree.Tree, len(bestTrees))
		for i, t := range bestTrees {
			result := CloseCycle(t, taxa[i], *aln, i)
			cycles[i] = result
			fmt.Println(result)
		}
		finalNetwork := AssembleNetwork(sntree, cycles)
		WriteTree(fmt.Sprintf("%s/final_network.nwk", args.polytomyDir), finalNetwork, true)
	}
}

func parseArgs() args {
	flag.NewFlagSet("Level-1 Network", flag.ContinueOnError)
	alnFile := flag.String("a", "", "alignment file")
	polytomyDir := flag.String("d", "", "directory with polytomy (created if using setup mode")
	setup := flag.Bool("s", false, "setup mode")
	flag.Parse()
	if *alnFile == "" || *polytomyDir == "" {
		fmt.Fprintln(os.Stderr, "both -a and -d are required")
		flag.Usage()
		os.Exit(1)
	}
	return args{alignmentFile: *alnFile, polytomyDir: *polytomyDir, setup: *setup}
}

func readAlignment(alnFile string) (*align.Alignment, error) {
	f, err := os.Open(alnFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	aln, err := nexus.NewParser(f).Parse()
	return &aln, err
}

func WriteTree(name string, t *tree.Tree, network bool) {
	var nwk string
	if network {
		nwk = fixNetwork(t.Newick())
	} else {
		nwk = t.Newick()
	}
	err := os.WriteFile(name, []byte(nwk), 0644)
	if err != nil {
		panic(fmt.Errorf("could not write file: %w", err))
	}
}

func ReadTree(name string) *tree.Tree {
	f, err := os.Open(name)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	t, err := newick.NewParser(f).Parse()
	if err != nil {
		panic(err)
	}
	return t
}

func fixNetwork(nwk string) string {
	nwk = strings.ReplaceAll(nwk, "[", "")
	nwk = strings.ReplaceAll(nwk, "]", "")
	return nwk
}
