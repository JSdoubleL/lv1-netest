# Level-1 Network Estimator

Software for estimating level-1 networks from two-state characters.

## Dependencies

Requires [go](https://go.dev/doc/install), [gotree](https://github.com/evolbioinfo/gotree), and [goalign](https://github.com/evolbioinfo/goalign). Remaining dependencies can be installed with `go get .` after go has been installed. 

[PAUP*](http://phylosolutions.com/paup-test/) is also required to run the full pipeline.

## Usage

Software can be used either by running `go run . <arguments>` from inside the software's directory, or `go build .` and then use the executable. 

First, run software in "setup" mode by using the `-s` argument. `-a` is for setting the input alignment file, and `-d` the output directory name. For example:

```sh
lv1-netest -a testdata/cycle.nex -d cycle -s
```

Then run PAUP* on all in `.nex` files in the output directory. Simply, pass each file as the only argument to PAUP*; all the necessary settings to run PAUP* appropriately are included in each file. An example of a loop running PAUP* on every file  is included in `run_paup.sh`. It is important that all the output files are contained in the same directory for the next step.

To get the final output, then run:

```sh
lv1-netest -a testdata/cycle.nex -d cycle 
```

