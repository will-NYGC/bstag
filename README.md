# bstag
BStag protocol for whole genome bisulfite sequencing on the Illumina X platform. Two versions of the filter script are included, `filterFillIn2.pl` serially processes each line of the BAM. `filterFillIn3.pl` chunks the BAM and processes the chunks in parallel.

## Usage:

```bash
/filterFillIn(2|3).pl -b <BAM> -r <REF FASTA> -c <CHROMS>
		          -b  BAM alignment
	                  -r  Reference FASTA
                          -p  Output prefix [Default: trims *.bam]
                          -k  Keep CHH state strings [Outputs to <PREFIX>.chhStates.txt
                          -d  Fill-in distance [Default: 11bp]
[filterFillIn3.pl]        -n  Number of cores [Default: 1]
[filterFillIn3.pl]        -c  Number of lines to chunk to cores [Default: 1e5]
                          -v  Verbose mode
		          -h  Print this message
```

|Flag|Description|
|----|-----------|
|-b|BAM alignment file produced by BWA-meth|
|-r|Indexed reference FASTA file used for alignment|
|-p|Output prefix. By default, replaces the `*.bam` extension with `*.filterFillIn`, and produces output `*.filterFillIn.(bam\|stats)`|
|-k|Keep CHH state strings. Creates a file, PREFIX.chhStates.txt, describing the sequence of CHH calls per read with a `|` marking the expected fill-in distance`|
|-d|Fill-in distance of end repair following transposition--default is 11bp|
|-v|Verbose mode|
|-h|Print a help message|
