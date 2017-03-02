#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $dir;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
use File::Basename;
# ($name,$path,$suffix) = fileparse($dir,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $which_read = "R1";
my $suffix     = ".fastq.gz";
my $self = bless {};
my $bsexpress         = "$ENV{HOME}/cegx-bfx-cegx-bsexpress-base-db6e49e95b74/scripts/bsExpress";
$self->{controls}{SQ} = "$ENV{HOME}/cegx-bfx-cegx-bsexpress-base-db6e49e95b74/control_reference/oxBS_controls-v1.0.fa";
$self->{controls}{DC} = "$ENV{HOME}/cegx-bfx-cegx-bsexpress-base-db6e49e95b74/digestion_reference/DC_controls-v1.0.fa";
my $type = 'SQ:DC';
my $skip_trim;
my $three_prime_clip_R1;
my $quality = 20;

GetOptions(
	   'd|dir:s'              => \$dir,
           'which_read:s'         => \$which_read,
           'suffix:s'             => \$suffix,
           'bsexpress:s'          => \$bsexpress,
           'three_prime_clip_R1:s'=> \$three_prime_clip_R1,
           'quality:s'            => \$quality,
           'sq|SQ_controls:s'     => \$self->{controls}{SQ},
           'dc|DC_controls:s'     => \$self->{controls}{DC},
           'type:s'               => \$type,
           'skip_trim'            => \$skip_trim,
           'debug'                => \$debug,
           'verbose'              => \$verbose,
           'simulate'             => \$simulate,
          );

$cmd = "find -L $dir -name \"*$which_read*$suffix\"";
print STDERR "# $cmd\n";
$ret = `$cmd`; chomp $ret;

my @files = split("\n",$ret);
my @types = split(":",$type);

if (scalar @files < 1) {
  print STDERR "No input files found.\n";
  exit(0);
} else {
  foreach my $file (sort @files) {
    foreach my $type (@types) {
      my $this_control = $self->{controls}{$type};
      $DB::single=1;1;
      die "Cannot find control. Check controls options -- $!" unless (-s $this_control);
      my ($this_name,$this_path,$this_suffix) = fileparse($file,$suffix);
      my $outdir = "$this_path/$this_name.runqc_".$type;
      my $prefix =            "$this_name.runqc_".$type;
      my $this_maxlen = 60; $this_maxlen = 100 if ($type =~ /DC/);
      $DB::single=1;1;
      my $skip_trim_opt = ''; $skip_trim_opt = '--skip_trim' if ($skip_trim);
      my $verbose_opt   = ''; $verbose_opt   = '--verbose'   if ($verbose);
      my $three_prime_clip_R1_opt   = ''; $three_prime_clip_R1_opt   = "--three_prime_clip_R1 $three_prime_clip_R1" if ($three_prime_clip_R1);
      my $trim_galore_opt = "--trim_galore_opts=\'--quality $quality $three_prime_clip_R1_opt\'";
      $cmd = "$bsexpress -i $file -r $this_control -p $prefix --outdir $outdir --maxlen $this_maxlen --skip_fastqc $skip_trim_opt $verbose_opt $trim_galore_opt";
      print STDERR "# $cmd\n";
      $ret = `$cmd`; chomp $ret;
    }
  }
}

$DB::single=1;1;
$DB::single=1;1;
