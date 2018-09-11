# X-WGBS fill-in filtering
X-WGBS protocol for whole genome bisulfite sequencing on the Illumina X platform. Two versions of the filter script are included, `filterFillIn2.pl` serially processes each line of the BAM. `filterFillIn3.pl` chunks the BAM and processes the chunks in parallel.

## Dependencies
|Software/Perl Module|
|--------|
|samtools v1.7|
|Bio::DB::Bam::Alignment|
|Getopt::Long|
|IPC::Cmd qw(can_run)|
|List::Util qw(min max)|
|Parallel::ForkManager|
|Fcntl qw(:flock SEEK_END)|

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
|-k|Keep CHH state strings. Creates a file, PREFIX.chhStates.txt, describing the sequence of CHH calls per read with a `\|` marking the fill-in distance position`|
|-d|Fill-in distance of end repair following transposition--default is 11bp|
|-v|Verbose mode|
|-h|Print a help message|

## Example
```
./filterFillIn2.pl -b test/na12878_10k_chr1.bam -r hg19.fa -p test/na12878_10k_chr1.filterFillIn -k -v
./filterFillIn3.pl -b test/na12878_10k_chr1.bam -r hg19.fa -p test/na12878_10k_chr1.filterFillIn -k -n 2 -c 1000 -v
```
