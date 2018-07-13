#!/bin/env perl
use strict;
use warnings;
use Bio::DB::Bam::Alignment;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw(can_run);
use IPC::Open2;
use List::Util qw(min max);
use POSIX;
use Parallel::ForkManager;
use Fcntl qw(:flock SEEK_END);


my $script_version = 'filterFillIn2.pl v0.3';

my $cmd = join(" ",$0,@ARGV);


# Get arguments
my($bam,$ref_fa,$prefix,$keep_string,$verbose,$help);
my $fill_dist = 11;
my $cores = 1;
my $chunk_size = 1e5;
GetOptions("b=s" => \$bam,
           "r=s" => \$ref_fa,
           "p=s" => \$prefix,
           "k"   => \$keep_string,
           "d:i" => \$fill_dist,
           "n:i" => \$cores,
           "c:i" => \$chunk_size,
           "v"   => \$verbose,
           "h"   => \$help);


# Check/define variables 
&usage if $help or !$bam or !$ref_fa;


# Needs samtools in the path
(print STDERR "***ERROR: \'samtools\' must be in path, exiting...\n\n" and 
    &usage) if !(can_run('samtools'));


# Define output prefix and files
if (!$prefix) {
    $prefix = $bam;
    $prefix =~ s/\.bam$/\.filterFillIn/;
    print STDERR "Using prefix: " . basename($prefix) . "\n";
}


# Open file handles
my $output_fh_bam = "$prefix.bam";
open(my $output_fh, "|-", "samtools sort -@ $cores -O BAM - > $output_fh_bam");

my $chh_file = "$prefix.chhString.txt";
open(my $chh_string_fh, ">", $chh_file) if $keep_string;


# Define minimum value dependent on machine-type
# If ^E, then XTen min = 2/#; else 2500 and min = 0/!
my $min = (`samtools view $bam | head -1 | cut -f 1` =~ /^E/) ? "#" : "!";  


# Read in BAM file
my $sam = Bio::DB::Sam -> new(-fasta=>$ref_fa,-bam=>$bam,-autoindex=>1);


# Print header with additional filterFillIn.pl line
print $output_fh &print_header($sam -> header -> text,
                               $script_version,
                               $fill_dist,
                               $cmd);


# Initiate counting variables
my $mapped = 0;
my $unmapped = 0;
my $filtered = 0;
my @chunk;


# Create iterator for BAM file
my $i = $sam -> features(-iterator=>1);


# Create parallel environment using $cores
# Apply function to update counting variables and print filtered lines
my $pm = Parallel::ForkManager -> new($cores);
$pm -> run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal,
        $core_dump, $hash_ref) = @_;
    $mapped += $hash_ref->{mapped};
    $unmapped += $hash_ref->{unmapped};
    $filtered += $hash_ref->{filtered};
    my $chh_string = $hash_ref->{chh_string};
    my @new_chunk = @{$hash_ref->{new_chunk}};

    if ($keep_string) {
        flock($chh_string_fh, LOCK_EX);
        seek($chh_string_fh, 0, 2);    
        print $chh_string_fh "$chh_string\n";
        flock($chh_string_fh, LOCK_UN); 
    }

    flock($output_fh, LOCK_EX);
    seek($output_fh, 0, 2);    
    print $output_fh join("\n", @new_chunk) . "\n";
    flock($output_fh, LOCK_UN);
});

# Iterate through BAM (with or without verbosity)
my $iter = 0;
while (my $read = $i->next_seq) { 
    # Would be much better to pre-chunk the bam rather than reading through
    # to create chunks, then re-reading the chunks to filter artifacts
    # But with Bio::DB::Sam, would require reading pre-chunks
    # into memory, e.g. all reads by chromosome--this is too memory intense
    if ($verbose) {
        printf("#\t%s | Filtering %s ..... %d\n",
               strftime("%F, %r",localtime),
               basename($bam),$iter) if $iter % 10e6 == 0;
        ++$iter;
    }
    push(@chunk, $read);

    # Create chunks of aligned read objects to reduce flock'ing overhead
    if (scalar @chunk == $chunk_size) {
        my @pass_chunk = @chunk;
        @chunk = ();

        ## Iterator could stop before last chunk is filled so process outside of loop
       $pm -> start and next;

        my($m, $u, $f, $chh_string, @new_chunk) = &process_chunk(\@pass_chunk);
        $pm -> finish(0, { mapped => $m,
                           unmapped => $u,
                           filtered => $f,
                           chh_string => $chh_string,
                           new_chunk => \@new_chunk });
    }
}
$pm -> wait_all_children;


# Process last chunk, if necessary, and close $output_fh
if (scalar @chunk > 0) {
    my($m, $u, $f, $chh_string, @new_chunk) = &process_chunk(\@chunk);
    $mapped += $m;
    $unmapped += $u;
    $filtered += $f;
    print $output_fh join("\n", @new_chunk) . "\n";
    print $chh_string_fh "$chh_string\n" if $keep_string;
}
close($output_fh);
# END #####


# SUBS ====================================
sub usage {
    my $usage = <<"USAGE";
$0 -b <BAM> -r <REF FASTA> -c <TO_SAM_PREFIX>
        -b  BAM alignment
        -r  Reference FASTA
        -p  Output prefix [Default: trims *.bam]
        -k  Keep CHH state strings [Outputs to <PREFIX>.chhStates.txt
        -d  Fill-in distance [Default: 11bp]
        -n  Number of cores [Default: 1]
        -c  Number of lines to chunk to cores [Default: 1e5]
        -v  Verbose mode
        -h  Print this message
USAGE
    print $usage;
    exit(1);
}

sub rev_comp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub print_header {
    my($header, $version, $fill_dist, $cmd) = @_;
    chomp $header; # for clarity, remove newline and then place in here-doc
    my $new_header = <<"HEADER";
$header
\@PG\tID:filterFillIn.pl\tPN:filterFillIn.pl\tVN:$script_version\tDS:[Reads with 4+ mCHH outside of 1st ${fill_dist}bp are bit-flagged (512/0x200) and the filter status (0|1) and CHH string (before/after ${fill_dist}bp fill-in: "|") are captured by the "YF:Z:(0|1):C|G{#}" tag. When R1 overlaps the R2 ${fill_dist}bp fill-in region, set R2 base qualities to 0 (2500) or 2 (X Ten) and report count in flag "YO:i:#"]\tCL:$cmd
HEADER

    return($new_header)
}

sub process_chunk {
    my @chunk = @{$_[0]};

    my $mapped = 0;
    my $unmapped = 0;
    my $filtered = 0;
    my($chh_string, @output);
    foreach my $read (@chunk) {
        # Extract stuff from align object
        my $bam_line = $read -> tam_line;
        my @bam_line = split(/\t/, $bam_line);
     
        # Only process if mapped, otherwise just push original read to output
        my $new_bam_line;
        if ( $bam_line[1] & 4 ) {
            ++$unmapped;
            $new_bam_line = $bam_line;

        # If mapped, filter fill-in artifacts
        } else {
            ++$mapped;
            my $isRevComp = $read -> reversed;
            # Get read from bitflag (note, actual R2's are flagged as 
            # R1 because swapped during bwa-meth alignment)
            my $r12 = ($bam_line[1] & 64) ? 'R1' : 'R2';

            # Get ref/match_string/query padded sequences
            my($ref,$matches,$query) = $read -> padded_alignment;
            if ($isRevComp) {
                $ref = &rev_comp($ref);
                $query = &rev_comp($query);
                $matches = reverse($matches);
            }
            my @ref = split('',$ref);
            my @query = split('',$query);

            # Zero R1 (actual R2) BQ's and add aux tag to line if R1/R2 overlap
            (my $clip_n, @bam_line) = &zero_r2_fill_in(\@bam_line, $fill_dist, 
                                                       $isRevComp, $r12);
            
            my $boundary_index = &get_boundary_index(\@query, $clip_n, 
                                                     $fill_dist, $r12);
            (my $f, $chh_string, @bam_line) = 
                &filter_fill_in(\@bam_line, \@ref, \@query, $r12, $boundary_index);

            $filtered += $f;
            $new_bam_line = join("\t", @bam_line);
        }
        push(@output, $new_bam_line);
    }
    return($mapped, $unmapped, $filtered, $chh_string, @output);
}

sub zero_r2_fill_in {  # R2 looks like R1 in these BAMs
    my ($bam_line_array_ref, $fill_dist, $isRevComp, $r12) = @_;

    my @bam_line = @$bam_line_array_ref;
    my($insert,$seq,$qual) = @bam_line[8..10];
    my $clip_n = min(length($seq), max(0,($fill_dist + length($seq)) - abs($insert)));

    # Define last INCLUDED offset in read and use later to remove 
    # readthrough fillIn effects on R1 
    my $last_index = length($seq)-1 - $clip_n;

    # R1 (actually R2) may capture fillIn artifact at its 3' end which could have 
    # higher base quality and get used over R2 (actually R1)
    # So "zero" at R2 BQ's to prevent this
    if ($r12 eq 'R1' and $clip_n > 0) {  # R1 is actually R2
        my @qual = ($isRevComp) ? reverse(split("",$qual)) : split("",$qual);
        $bam_line[10] = ($isRevComp) ? join("",reverse(@qual[0..$last_index],
                                            $min x min($clip_n,length($seq)))) : 
                                       join("",@qual[0..$last_index],
                                            $min x min($clip_n,length($seq)));

        push(@bam_line,"YO:i:$clip_n");
    }
    return($clip_n, @bam_line);
}

sub get_boundary_index {
    my($query_array_ref, $clip_n, $fill_dist, $r12) = @_;

    my @query = @$query_array_ref;

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
        for (my $i = 0; $i < scalar @query; $i++) {
            ++$counter if $query[$i] ne '-';
            $boundary_index = $i and last if $counter == $fill_dist;
        }
    }
    return($boundary_index);
}

sub filter_fill_in {
    my($bam_line_array_ref, $ref_array_ref, $query_array_ref,
       $r12, $boundary_index) = @_;

    my @bam_line = @$bam_line_array_ref;
    my @ref = @$ref_array_ref;
    my $ref = join("", @ref);
    my @query = @$query_array_ref;

    my $base = ($r12 eq 'R1') ? 'C' : 'G';
    my $context = ($base eq 'C') ? '[ACT]' : '[AGT]';
    my $regex = ($r12 eq 'R1') ? "$base$context\{2\}" : "$context\{2\}$base";

    my $offset = 0;
    my $result = index($ref, $base, $offset);  # Initiate, then while loop indexes		
    
    my(@chh_before, @chh_after);
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
        # Keep CHH's before and after boundary	
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
    
    # 4 or more mCHH outside of fill-in is artifact
    # n=4; 1 - pbinom(n,p=0.3,q=n-1) = 0.0081
	# highest reported meth ratios of CHH are in brain, 20-25%
	# so used conservative p=0.3, although p.value<0.01 if p=0.2 or p=0.25
    my $filtered = 0;
    if ($chh_filter_string =~ /^$base{4}/) {  
        $bam_line[1] += 512 if !( $bam_line[1] & 512 );
        $filtered = 1;
    }
    return($filtered, $chh_string, @bam_line)
}