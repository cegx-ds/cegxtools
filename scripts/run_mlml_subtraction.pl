#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile_mkBS;
my $inputfile_oxBS;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
use File::Basename;
# ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $suffix = '.methcounts.bed.gz';
my $mlml = "$ENV{HOME}/methpipe/bin/mlml";
my $tag = 'mlml';
my $skip_mlml;
my $bedtools = 'bedtools';

GetOptions(
           'mkbs|mkBS|inputfile_mkBS:s' => \$inputfile_mkBS,
           'oxbs|oxBS|inputfile_oxBS:s' => \$inputfile_oxBS,
           'suffix:s'                   => \$suffix,
           'tag:s'                      => \$tag,
           'mlml:s'                     => \$mlml,
           'bedtools:s'                 => \$bedtools,
           'skip_mlml'                  => \$skip_mlml,
           'debug'                      => \$debug,
           'verbose'                    => \$verbose,
           'simulate'                   => \$simulate,
          );

if ($suffix eq '.methcounts.bed.gz') {
  die "undefined or empty file inputfile_mkBS -- $!" if (!defined($inputfile_mkBS) || !(-s $inputfile_mkBS));
  die "undefined or empty file inputfile_oxBS -- $!" if (!defined($inputfile_oxBS) || !(-s $inputfile_oxBS));
  my ($mkBS_name,$mkBS_path,$mkBS_suffix) = fileparse($inputfile_mkBS,$suffix);
  my ($oxBS_name,$oxBS_path,$oxBS_suffix) = fileparse($inputfile_oxBS,$suffix);

  my $mkBS_outmeth = "$mkBS_path/$mkBS_name.meth";
  $cmd = "gunzip -c $inputfile_mkBS > $mkBS_outmeth";
  print STDERR "# $cmd\n";
  $ret = `$cmd` unless ($simulate || $skip_mlml);

  my $oxBS_outmeth = "$oxBS_path/$oxBS_name.meth";
  $cmd = "gunzip -c $inputfile_oxBS > $oxBS_outmeth";
  print STDERR "# $cmd\n";
  $ret = `$cmd` unless ($simulate || $skip_mlml);

  my $outfile = "$mkBS_path/$mkBS_name.vs.$oxBS_name.$tag.txt";
  my $logfile = "$mkBS_path/$mkBS_name.vs.$oxBS_name.$tag.log";
  my $errfile = "$mkBS_path/$mkBS_name.vs.$oxBS_name.$tag.err";
  $cmd = "$mlml -v -bsseq $mkBS_outmeth -oxbsseq $oxBS_outmeth -output $outfile 1>$logfile 2>$errfile";
  print STDERR "# $cmd\n";
  $ret = `$cmd` unless ($simulate || $skip_mlml);

  $cmd = "rm -f $mkBS_outmeth $oxBS_outmeth";
  print STDERR "# $cmd\n";
  $ret = `$cmd` unless ($simulate || $skip_mlml);


  ########################################
  # Transform the output file into a bed file for intersection with input
  my $outfile_s__mC = "$mkBS_path/$mkBS_name.vs.$oxBS_name.$tag.5mC.CpG_report.txt.gz";
  my $outfile_s_hmC = "$mkBS_path/$mkBS_name.vs.$oxBS_name.$tag.5hmC.CpG_report.txt.gz";
  open(my $out_s__mC, "|-", " gzip -c > $outfile_s__mC") or die $!;
  open(my $out_s_hmC, "|-", " gzip -c > $outfile_s_hmC") or die $!;

  $cmd = "cat $outfile" . q{ | awk '{print $1"\t"$2"\t"$3"\t"$4":"$5":"$6":"$7"\t1\t."}' } . "| $bedtools intersect -sorted -nobuf -wb -a stdin -b $inputfile_mkBS | $bedtools intersect -sorted -nobuf -wb -a stdin -b $inputfile_oxBS";
  print STDERR "# $cmd\n";
  open ISECT, "$cmd |" or die $!;
  while (<ISECT>) {
    my $line = $_; chomp $line;
    my ($mlml_chr, $mlml_start, $mlml_end, $mlml_name, $mlml_score, $mlml_strand,
	$mkbs_chr, $mkbs_start, $mkbs_end, $mkbs_name, $mkbs_score, $mkbs_strand,
       	$oxbs_chr, $oxbs_start, $oxbs_end, $oxbs_name, $oxbs_score, $oxbs_strand) = split("\t",$line);

    my ($s__mC,$s_hmC,$s___C,$confl)  = split(":",$mlml_name);
    my ($mkbs_context,$mkbs_coverage) = split(":",$mkbs_name);
    my ($oxbs_context,$oxbs_coverage) = split(":",$oxbs_name);
    $DB::single=1;1;
    my $count_methylated__mC = 0; $count_methylated__mC = $s__mC*($mkbs_coverage+$oxbs_coverage) if ($s__mC ne 'nan' && $s__mC > 0);
    my $count_methylated_hmC = 0; $count_methylated_hmC = $s_hmC*($mkbs_coverage+$oxbs_coverage) if ($s_hmC ne 'nan' && $s_hmC > 0);
    my $count_unmetlated___C = 0; $count_unmetlated___C = $s___C*($mkbs_coverage+$oxbs_coverage) if ($s___C ne 'nan' && $s___C > 0);
    my $out__mC = "$mlml_chr\t$mlml_end\t$mkbs_strand\t$count_methylated__mC\t$count_unmetlated___C\tCG\tCGN\n";
    my $out_hmC = "$mlml_chr\t$mlml_end\t$mkbs_strand\t$count_methylated_hmC\t$count_unmetlated___C\tCG\tCGN\n";
    print $out_s__mC $out__mC;
    print $out_s_hmC $out_hmC;
    # example case
    #chr1  969041   969042   0.416667:0.483333:0.1:0          .  chr1  969041   969042   CpG:56   0.892857142857143    -  chr1  969041   969042   CpG:24   0.416666666666667    -
    # echo "0.416666666666667*(56+24)" | bc -l
    # echo "0.483333*(56+24)" | bc -l
    # echo "0.1*(56+24)" | bc -l
  }
  close ISECT;
  close $out_s__mC;
  close $out_s_hmC;
  print STDERR "outfile_s__mC:\n";
  print        "$outfile_s__mC\n";
  print STDERR "outfile_s_hmC:\n";
  print        "$outfile_s_hmC\n";

  # Clean up .txt file from mlml
  $cmd = "rm -f $outfile";
  print STDERR "# $cmd\n";
  $ret = `$cmd` unless ($simulate);


} else {
  print STDERR "subtraction not implemented for suffix $suffix\n";
  exit(0);
}

# The genome-wide cytosine methylation output file is tab-delimited in the following format:
# ==========================================================================================
# <chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>

# Usage: mlml [OPTIONS]
# 
# Options:
#   -o, -output     Name of output file (default: stdout) 
#   -u, -bsseq      Name of input BS-Seq methcounts file 
#   -h, -tabseq     Name of input TAB-Seq methcounts file 
#   -m, -oxbsseq    Name of input oxBS-Seq methcounts file 
#   -t, -tolerance  EM convergence threshold. Default 1e-10 
#   -a, -alpha      significance level of binomial test for each site. Default 
#                   0.05 
#   -v, -verbose    print run statistics 
# 
# Help options:
#   -?, -help       print this help message 
#       -about      print about message 
# 
# PROGRAM: mlml

########################################

# Output
# chrom   start   end     pm      ph      pu      #_of_conflict

# Log
# Total sites: 58806012
# Sites with overshoot: 6323762 (10.7536%)
# Sites conflicting to at least two input: 646 (0.00109853%)

# The columns are:

# chromosome name
# start position
# end position
# 5-mC level
# 5-hmC level
# unmethylated level
# number of conflicts

# To calculate the last column, a binomial test is performed for each
# input methylation level (can be 2 or 3 in total depending on
# parameters). If the estimated methylation level falls out of the
# confidence interval calculated from input coverage and methylation
# level, then such event is counted as one conflict. It is recommended
# to filter estimation results based on the number of conflicts; if
# more conflicts happens on one site then it is possible that
# information from such site is not reliable.

# run_mlml_subtraction.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

run_mlml_subtraction.pl - DESCRIPTION 

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

GetOptions(
           'mkbs|mkBS|inputfile_mkBS:s' => \$inputfile_mkBS,
           'oxbs|oxBS|inputfile_oxBS:s' => \$inputfile_oxBS,
           'suffix:s'                   => \$suffix,
           'tag:s'                      => \$tag,
           'mlml:s'                     => \$mlml,
           'bedtools:s'                 => \$bedtools,
           'skip_mlml'                  => \$skip_mlml,
           'debug'                      => \$debug,
           'verbose'                    => \$verbose,
           'simulate'                   => \$simulate,
          );

=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


1;
