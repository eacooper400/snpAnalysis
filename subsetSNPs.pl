#!/usr/bin/perl -w
#
# subsetSNPs.pl
#
# Get a specified subset of individuals into a new snp file
#
# October 18, 2013
# Liz Cooper

use strict;

# Take as input arguments the following:
# The input file in .snps format, the output file name,
# A list of individuals to split into pops, specied as pop.ind,pop.ind, with spaces between different output pops
my ($USAGE) = "\n$0 <input.snps> <output.snps> <pop1.ind1,pop1.ind2> <pop2.ind1,pop2.ind2>
\tinput.snps = The input file in .snps format
\toutput.snps = The name of an output file also in .snps format
\t A separate list of individuals to include in each subset population
\t\tPopulation and Individual Numbers should be 0-based
\t\tAny number of new output populations may be specified\
\t\tIndividual IDs may be re-used\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my $input = $ARGV[0];
my $output = $ARGV[1];
my @list = @ARGV[2..(scalar @ARGV - 1)];

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  print OUT $info[0], "\t", $info[1], "\t", $info[2];
  
  foreach my $newpop (@list) {
    my @inds = split(/,/, $newpop);
    my $newstring = '';
    
    foreach my $i (@inds) {
      my ($p, $n) = split(/\./, $i);
      my $genotype = substr($info[$p+3], $n, 1);
      $newstring .= $genotype;
    }

    print OUT "\t", $newstring;
  }
  print OUT "\n";
}
close(IN);
close(OUT);
exit;
