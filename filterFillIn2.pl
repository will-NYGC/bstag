#!/bin/env perl
# Needs samtools in the path
use strict;
use warnings;
use Bio::DB::Bam::Alignment;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw(can_run);
use List::Util qw(min max);
use POSIX;


my $script_version = 'filterFillIn2.pl v0.2';

my $cmd = join(" ",$0,@ARGV);

my($bam,$ref_fa,$prefix,$keep_string,$verbose,$help);
my $fillIn_dist = 11;
GetOptions("b=s" => \$bam,
	   "r=s" => \$ref_fa,
	   "p=s" => \$prefix,
	   "k"   => \$keep_string,     # Keep string of mC/C on reads for downstream QC
	   "d:i" => \$fillIn_dist,
	   "v"   => \$verbose,
	   "h"   => \$help);

&usage if $help or !$bam or !$ref_fa;


# Test for samtools
(print STDERR "***ERROR: \'samtools\' must be in path, exiting...\n\n" and &usage) if !(can_run('samtools'));


# Define output prefix and files
if (!$prefix) {
    $prefix = $bam;
    $prefix =~ s/\.bam$/\.filterFillIn/;
    print STDERR "Using prefix: " . basename($prefix) . "\n";
}

my $output_bam = "$prefix.bam";
open(OUTPUT,"| samtools view -hb - > $output_bam");

my $chh_file = "$prefix.chhString.txt" if $keep_string;
open(CHH_STRING,">$chh_file") if $keep_string;


# Define minimum value dependent on machine-type
my $min = (`samtools view $bam | head -1 | cut -f 1` =~ /^E/) ? "#" : "!";  # If ^E, then XTen min = 2/#; else 2500 and min = 0/!

# Read in BAM file
my $sam = Bio::DB::Sam -> new(-fasta=>$ref_fa,-bam=>$bam,-autoindex=>1);


# Print header with additional filterFillIn.pl line
my $header = $sam -> header -> text;
print OUTPUT "$header";
print OUTPUT "\@PG\tID:filterFillIn.pl\tPN:filterFillIn.pl\tVN\:$script_version\t";
print OUTPUT "DS:[Reads with 4+ mCHH outside of 1st ${fillIn_dist}bp of read are flagged as \"not passing filters\" by marking with 512/0x200 ";
print OUTPUT "bitwise flag with the CHH string included (with \"|\" after ${fillIn_dist}bp) in the \"YF:Z:C|G{#}\" tag. ";
print OUTPUT "When R1 overlaps R2 ${fillIn_dist}bp fill-in region, set R2 base qualities to 0 (2500) or 2 (X Ten) and report count in flag \"YO:i:#\"]";
print OUTPUT "\tCL\:$cmd\n";



# Parse reads by chromosome
my $mapped = 0;
my $unmapped = 0;
my $filtered = 0;

my $i = $sam -> features(-iterator=>1);

my $iter = 0 if $verbose;
while (my $read = $i->next_seq) { 
		
    if ($verbose) {
	printf("#\t%s | Filtering %s ..... %d\n",strftime("%F, %r",localtime),basename($bam),$iter) if $iter % 10e6 == 0;
	++$iter;
    }

    # Reconstruct SAM format fields	
    my @bam_line = split(/\t/,$read -> tam_line);

    if ( $bam_line[1] & 4 ) {
	++$unmapped;
	print OUTPUT join("\t",@bam_line) . "\n";
	next;
    } else {
	++$mapped;
    }

    my $isRevComp = $read -> reversed;

    # Get read from bitflag (note, actual R2's are flagged as R1 because swapped furing bwameth.sh)
    # Tried to flip but breaks PileOmeth b/c that code expects R1 to have C>T and R2 to have G>A
    my $r12 = ($bam_line[1] & 64) ? 'R1' : 'R2';
    
    # Change R1 (actually R2) overlapping base qualities to min (#|!) when insert size less than readlength add change YO:i flag to 1
    # Min base quality will cause it to get skipped on pile up
    my($insert,$seq,$qual) = @bam_line[8..10];
    my $clip_n = min(length($seq),max(0,($fillIn_dist + length($seq)) - abs($insert)));

    # Define last INCLUDED offset in read and use later to remove readthrough fillIn effects on R1 
    my $last_index = length($seq)-1 - $clip_n;

    # R1 (actually R2) may capture fillIn artifact at its 3' end which could have higher base quality and get used over R2 (actually R1)
    if ($r12 eq 'R1' and $clip_n > 0) {  # R1 is actually R2
	my @qual = ($isRevComp) ? reverse(split("",$qual)) : split("",$qual);
	$bam_line[10] = ($isRevComp) ? join("",reverse(@qual[0..$last_index],$min x min($clip_n,length($seq)))) : 
	    join("",@qual[0..$last_index],$min x min($clip_n,length($seq)));

	push(@bam_line,"YO:i:$clip_n");
    }	

    
    my($ref,$matches,$query) = $read -> padded_alignment;
    if ($isRevComp) {
	$ref = &rev_comp($ref);
	$query = &rev_comp($query);
	$matches = reverse($matches);
    }
    
    my @ref = split('',$ref);
    my @query = split('',$query);

    
    # Establish boundary index in context of query--matches--reference strings which may or may not be length of read (indels will expand it)
    my $boundary_index;    
    if ($r12 eq 'R1') {
	if ($clip_n > 0) {
	    my $counter = 0;
	    for (my $i = $#query; $i >= 0; $i--) {
		++$counter if $query[$i] ne '-';
		$boundary_index = $i and last if $counter == $clip_n;
	    }
	} else {
	    $boundary_index = $#query + 1;
	}
    } else {
        my $counter = 0;
        for (my $i = 0; $i < length($query); $i++) {
            ++$counter if $query[$i] ne '-';
            $boundary_index = $i and last if $counter == $fillIn_dist;
        }
    }

    # Declare arrays for CHH sequences
    # @chh_before: before fill-in region, <=11bp (index=10), typically
    # @chh_after: after fill-in region, >11bp (index=10)
    # @chh_filter: precise string after fill-in, including '-' for deleted reference C's
    my(@chh_before,@chh_after);
        
    my $base = ($r12 eq 'R1') ? 'C' : 'G';
    my $context = ($base eq 'C') ? '[ACT]' : '[AGT]';
    my $regex = ($r12 eq 'R1') ? "$base$context\{2\}" : "$context\{2\}$base";

    my $offset = 0;
    my $result = index($ref,$base,$offset);  # Initiate, then while loop indexes		
    
    while ($result != -1) {  #As long as there's a match and offset has room for trinucleotide
	
	my $tri = $ref[$result];
	
	if ($r12 eq 'R1') {
	    my $j = $result + 1;
	    while (length($tri)==1 and $j <= length($ref)-2) {
		$tri .= $ref[$j] if ($ref[$j] =~ /[ACGT]/);  # This is to skip "-" when deletion and get next letter
		++$j;
	    }
	    while (length($tri)==2 and $j <= length($ref)-1) {
		$tri .= $ref[$j] if ($ref[$j] =~ /[ACGT]/);
		++$j;
	    }
	} else {
	    my $j = $result - 1;
	    while (length($tri)==1 and $j >= 1) {
		$tri = "$ref[$j]$tri" if ($ref[$j] =~ /[ACGT]/);
		--$j;
	    }
	    while (length($tri)==2 and $j >= 0) {
		$tri = "$ref[$j]$tri" if ($ref[$j] =~ /[ACGT]/);
		--$j;
	    }
	}

	# Keep CHH's before and after 9th base separate	
	if ( $tri =~ /$regex/ ) {

	    # R1 can read into fill-in artifact, so must remove from 3' end, R2 is 5' end
	    # Will be reversed for R1, chh_string will read from end to front so consistent with R2 direction
	    if ($r12 eq 'R1') {
		($result >= $boundary_index) ? unshift(@chh_before,$result) : unshift(@chh_after,$result);
	    } else {
		($result <= $boundary_index) ? push(@chh_before,$result) : push(@chh_after,$result);
	    }
	}
	$offset = $result + 1;
	$result = index($ref,$base,$offset);
    }
    my $chh_string = join("",@query[@chh_before],"|",@query[@chh_after]);
    my $chh_filter_string = join("",@query[@chh_after]);
    
    push(@bam_line,"YF:Z:$chh_string");  # Adds this flag indicating CHH string, with "|" marking break before and after fill-in
    
    if ($chh_filter_string =~ /^$base{4}/) {  # 4 or more mCHH outside of fill-in is artifact
                                              # n=4; 1 - pbinom(n,p=0.3,q=n-1) = 0.0081
	                                      # highest reported meth ratios of CHH are in brain, 20-25%
	                                      # so used conservative p=0.3, although p.value<0.01 if p=0.2 or p=0.25
	$bam_line[1] += 512 if !( $bam_line[1] & 512 );
	++$filtered;
    }
    print OUTPUT join("\t",@bam_line) . "\n";    
    print CHH_STRING "$chh_string\n" if $keep_string and $chh_string;
}
close(CHH_STRING);



my $filtered_rate = ($mapped > 0) ? $filtered/$mapped : 'NA';

# Print this to stout so can redirect to $QC_DIR instead of $prefix directory
my $stats = "$prefix.stats";
open(STATS,">$stats");
print STATS join("\t","#Mapped",qw( Unmapped Map.filtered Map.filtered.rate )) . "\n";
print STATS join("\t",$mapped,$unmapped,$filtered,$filtered_rate) . "\n";




# SUBS ###########
sub usage {
    print STDERR "$0 -b <BAM> -r <REF FASTA> -c <CHROMS>\n";
    print STDERR "\t-b  BAM alignment\n";
    print STDERR "\t-r  Reference FASTA\n";
    print STDERR "\t-p  Output prefix [Default: trims *.bam]\n";
    print STDERR "\t-k  Keep CHH state strings [Outputs to <PREFIX>.chhStates.txt\n";
    print STDERR "\t-d  Fill-in distance [Default: 11bp]\n";
    print STDERR "\t-v  Verbose mode\n";
    print STDERR "\t-h  Print this message\n";
    exit;
}


sub rev_comp {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}
