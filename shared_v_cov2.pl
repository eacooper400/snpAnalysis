#!/usr/bin/perl -w
#
# shared_v_cov2.pl
#
# A new script to look at sharing versus coverage levels
# across any number of populations
#
# July 30, 2014
# Liz Cooper

use strict;

# Get the input and output file names from the command line
my ($USAGE) = "\n$0 <input.counts> <output.shared>
\tinput.counts = A text file of shared snp types created by SNP_categories.pl
\toutput.shared = A output file with a line for every coverage value, and the proportion of shared sites at that cov\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";

my %shared = ();
my %total = ();

while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  my $cat = $info[0];
  my @covs = @info[1..(scalar @info -1)];
  
  my $sum = 0;
  foreach my $c (@covs) {
    $sum += $c;
  }
  my $mean = int($sum/(scalar @covs));

  if (exists $total{$mean}) {
    $total{$mean} += 1;
  } else {
    $total{$mean} = 1;
  }

  if ($cat == 15) {
    if (exists $shared{$mean}) {
      $shared{$mean} += 1;
    } else {
      $shared{$mean} = 1;
    }
  }
}
close(IN);

open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my @cov_means = sort {$a <= $b} keys %total;

foreach my $coverage (@cov_means) {
  if (exists $shared{$coverage}) {
    my $prop = $shared{$coverage}/$total{$coverage};
    print OUT $coverage, "\t", $prop, "\t", $total{$coverage}, "\n";
  } else {
    print OUT $coverage, "\t", "0", "\t", "0\n";
  }
}
close(OUT);
exit;
