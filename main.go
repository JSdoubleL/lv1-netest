package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/nexus"
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
		fmt.Println("done.")
	} else {
		taxa := ReadTaxa(args.polytomyDir)
		fmt.Println(taxa)
		bestTrees := ReadPAUPResults(args.polytomyDir, uint(len(taxa)))
		fmt.Println(taxa, bestTrees)
		for i, t := range bestTrees {
			result := CloseCycle(t, taxa[i], *aln)
			fmt.Println(result)
		}
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
