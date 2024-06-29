package main

import (
	"fmt"
	"os"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/nexus"
)

func main() {
	alnFile := *parseArgs()
	aln, err := readAlignment(alnFile)
	if err != nil {
		panic(err)
	}
	fmt.Println(aln)
	sntree := SNTree(aln)
	fmt.Println(sntree)
	polytomies := ExtractPolytomies(sntree)
	fmt.Println(polytomies)
	WritePolytomies(polytomies, aln)

}

func parseArgs() *string {
	alnFile := os.Args[1]
	return &alnFile
}

func readAlignment(alnFile string) (align.Alignment, error) {
	f, err := os.Open(alnFile)
	if err != nil {
		return nil, err
	}
	aln, err := nexus.NewParser(f).Parse()
	return aln, err
}
