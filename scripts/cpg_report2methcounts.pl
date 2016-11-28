#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
my $suffix = '.CpG_report.txt.gz';
use File::Basename;
# ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $tag = 'methcounts';
my $mem = '1G';
my $format = 'bed';
my $gzip = "$ENV{HOME}/htslib/bin/bgzip";
my $merge_strands;

GetOptions(
           'i|input|inputfile:s'   => \$inputfile,
           'suffix:s'              => \$suffix,
           'tag:s'                 => \$tag,
           'mem:s'                 => \$mem,
           'gzip|bgzip:s'          => \$gzip,
           'format:s'              => \$format,
           'm|merge|merge_strands' => \$merge_strands,
           'debug'                 => \$debug,
           'verbose'               => \$verbose,
           'simulate'              => \$simulate,
          );

$mem .= 'G' if ($mem !~ /G|K$/);

my ($name,$path,$this_suffix) = fileparse($inputfile,$suffix);
$tag .= ".mstr" if ($merge_strands);
$tag .= ".bed" if ($format =~ /bed/);
my $outfile = "$path/$name.$tag.gz";
$cmd = "gunzip -c $inputfile";
my $buffer_size = $mem;
my $sort_opt = "sort -k 1,1 -k2,2n --buffer-size=$buffer_size --temporary-directory $path";
open(my $out_fh, "|-", "$sort_opt | $gzip -c > $outfile") or die $!;
print STDERR "$cmd | cpg_report2methcounts-transformation | $sort_opt\n";
open IN, "$cmd | " or die $!;
my ($prev_chromosome, $prev_start, $prev_position, $prev_strand, $prev_count_methylated, $prev_count_non_methylated, $prev_C_context, $prev_context, $prev_trinucleotide_context, $prev_total_reads, $count_methylated, $count_non_methylated, $prev_methlevel);
my $count = 0;
while (<IN>) {
  my $line = $_; chomp $line;
  my ($chromosome, $position, $strand, $count_methylated, $count_non_methylated, $C_context, $trinucleotide_context) = split("\t", $line);
  my $total_reads = $count_methylated+$count_non_methylated;
  my $methlevel = 0; $methlevel = $count_methylated/($total_reads) if ($count_methylated > 0);
  my $context = 'CpG';
  $context = 'CHH' if ($C_context eq 'CHH');
  $context = 'CHG' if ($C_context eq 'CHG');

  # print $out_fh "$chromosome\t$position\t$strand\t$context\t$methlevel\t$total_reads\n"  if ($format =~ /methcounts/);
  my $start = $position - 1;
  if ($merge_strands) {
    if (defined $prev_chromosome &&
	defined $prev_position &&
	defined $prev_strand &&
	defined $prev_C_context &&
	($prev_chromosome eq $chromosome) &&
	($prev_position == $position-1) &&
	($prev_strand eq '+') &&
	($strand      eq '-') &&
	($prev_strand ne $strand) &&
	($prev_C_context eq $C_context)
       ) {
      my $merged_total_reads          = $total_reads          + $prev_total_reads;
      my $merged_count_methylated     = $count_methylated     + $prev_count_methylated;
      my $merged_count_non_methylated = $count_non_methylated + $prev_count_non_methylated;
      my $merged_methlevel = 0;
      $merged_methlevel    = $merged_count_methylated/($merged_total_reads) if ($merged_count_methylated > 0);
      print $out_fh "$prev_chromosome\t$prev_start\t$prev_position\t$prev_context:$merged_total_reads\t$merged_methlevel\t$prev_strand\n" if ($format =~ /bed/);
    } else {
      # First strand count capture, nothing to be printed
    }
  } else {
    print $out_fh "$chromosome\t$start\t$position\t$context:$total_reads\t$methlevel\t$strand\n" if ($format =~ /bed/);
  }
  $prev_chromosome             = $chromosome;
  $prev_start                  = $start;
  $prev_position               = $position;
  $prev_strand                 = $strand;
  $prev_count_methylated       = $count_methylated;
  $prev_count_non_methylated   = $count_non_methylated;
  $prev_C_context              = $C_context;
  $prev_context                = $context;
  $prev_trinucleotide_context  = $trinucleotide_context;
  $prev_total_reads            = $total_reads;
  $prev_count_methylated       = $count_methylated;
  $prev_count_non_methylated   = $count_non_methylated;
  $prev_methlevel              = $methlevel;
  $count++;
}
close $out_fh;

print STDERR "outfile:\n";
print        "$outfile\n";

########################################
# From methpipe documentation

# The first column is the chromosome.

# The second is the location of the cytosine.

# The 3rd column indicates the strand, which can be either + or -.

# The 4th column is the sequence context of that site, followed by an
# x if the site has mutated in the sample away from the reference
# genome.

# The 5th column is the estimated methylation level, equal to the
# number of Cs in reads at position corresponding to the site, divided
# by the sum of the Cs and Ts mapping to that position.

# The final column is number of reads overlapping with that site.

# chr1 3000826 + CpG 0.852941 34
# chr1 3001006 + CHH 0.681818 44
# chr1 3001017 -CpG 0.609756 41
# chr1 3001276 + CpGx 0.454545 22
# chr1 3001628 -CHH 0.419753 81
# chr1 3003225 + CpG 0.357143 14
# chr1 3003338 + CpG 0.673913 46
# chr1 3003378 + CpG 0.555556 27
# chr1 3003581 -CHG 0.641026 39
# chr1 3003639 + CpG 0.285714 7

########################################

# From bismark documentation

# The genome-wide cytosine methylation output file is tab-delimited in the following format:
# ==========================================================================================
# <chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>

# CpG_report.txt.gz
# chrM    91      +       17      874     CG      CGA
# chrM    92      -       111     3949    CG      CGC
# chrM    96      +       2       953     CG      CGC
# chrM    97      -       131     4037    CG      CGT
# chrM    105     +       2       1007    CG      CGG
# chrM    106     -       98      4267    CG      CGG
# chrM    120     +       4       1064    CG      CGC
# chrM    121     -       98      4555    CG      CGA

# cpg_report2methcounts.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

cpg_report2methcounts.pl - DESCRIPTION 

perl cpg_report2methcounts.pl -i SAMPLE.CpG_report.txt.gz

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

GetOptions(
           'i|input|inputfile:s'   => \$inputfile,
           'suffix:s'              => \$suffix,
           'tag:s'                 => \$tag,
           'mem:s'                 => \$mem,
           'gzip|bgzip:s'          => \$gzip,
           'format:s'              => \$format,
           'm|merge|merge_strands' => \$merge_strands,
           'debug'                 => \$debug,
           'verbose'               => \$verbose,
           'simulate'              => \$simulate,
          );

=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut



$DB::single=1;1;
$DB::single=1;1;
