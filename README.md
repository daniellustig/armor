# ArMOR: Architecture-independent Memory Ordering Requirement specifications

## Project Overview

For an overview of ArMOR itself, see our ISCA paper:

Daniel Lustig, Caroline Trippel, Michael Pellauer, and Margaret Martonosi,
"ArMOR: Defending Against Memory Consistency Model Mismatches in Heterogeneous Architectures",
*42nd International Symposium on Computer Architecture (ISCA)*, Portland, OR, June 2015.
[Official link](http://dl.acm.org/citation.cfm?id=2750378)
[Local link](https://github.com/daniellustig/armor/blob/master/doc/dlustig_ISCA15_TR.pdf)

If you use ArMOR (and we hope you do!), we would also appreciate if you cite the above paper.

## Contacting the Authors

Please contact Daniel Lustig (dlustig@princeton.edu) with any questions you may have.

## How to Use ArMOR

### Prerequisites

- ArMOR itself is a self-contained Python app.  The authors have tested ArMOR using Python v2.7.
- The FSM figures are drawn using [GraphViz](www.graphviz.org).
- The summary report is compiled using LaTeX, and pdfLaTeX in particular.  Other variants should work, but may require minor script changes, etc.

### Building ArMOR

No need...it's just a Python script!

### Running ArMOR

- Run `./one.py` to build some simple examples.
- Run `./build_reports.sh` to build a set of documents summarizing the shim FSMs generated for a variety of use cases.
These reports look similar to (or, in one case, are identical to) the appendix of our tech report.
- Either script can be easily adapted to any other use cases you may be interested in.

### File summary

- `armor.py`: the core ArMOR routines
- `architectures.py`: a set of pre-defined architectures
- `one.py`: a script for generating one-off shim FSMs
- `report.py`: an ArMOR script for generating a report summarizing a number of shim FSMs
- `build_reports.sh`: a script for compiling the results of `report.py` into PDF summaries
