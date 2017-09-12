#!/usr/bin/perl -w
#
# permuteFst.pl
#
# Do a permutation test to get the confidence intervals for randomized Fsts
#
# October 4, 2013
# Liz Cooper

use strict;
use Math::Combinatorics;
use List::Util;
srand(time|$$);

# Take as input:
# The starting SNP file and an output file to save permuted Fsts
# The number of permutations to do
# The coverage and MAF cut-off values to use
my ($USAGE) = "\n$0 <input.snps> <output.perm> <num_iterations> <coverage_cutoff> <MAF_minimum>
\tinput.snps = The input file in .snps format
\toutput.perm = The output file to save all permuted Fst values
\tnum.iterations = The number of permutations to perform
\tcoverage_cutoff = The min. number individuals to keep a permuted population
\tMAF_minimum = The minimum MAF to keep a permuted population\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $iter, $mincov, $minmaf) = @ARGV;

# Print a header line to the output file
open (IN, $input) || die "\nUnable to open the file $input!\n";
my $first = <IN>;
chomp $first;
my @fields = split(/\t/, $first);
my @pops = @fields[3..(scalar @fields - 1)];
my $npops = scalar @pops;
my $n_comb = combine(2, @pops);
my @pop_string = ();
for (my $i = 1; $i <= $npops; $i++) {
  push (@pop_string, $i);
}
my $popstring = join(",", @pop_string);

open (OUT, ">$output") || die "\nUnable to open the file $output!\n";
print OUT "Iter";
for (my $c = 0; $c < $n_comb; $c++) {
  print OUT "\tFst", $c;
}
print OUT "\n";
close(IN);

# Set the iteration counter and start looping
while ($iter) {

  # Get the shuffled order of individuals
  my @inds = ();
    
  # Open a temporary SNP file to save shuffled genotypes
  open (SNP, ">temp.snps") || die "\nUnable to open the file temp.snps!\n";
  open (IN, $input) || die "\nUnable to open the file $input!\n";

  while (<IN>) {
    chomp $_;
    my @info = split(/\t/, $_);
    my @p = @info[3..(scalar @info - 1)];
    my @sizes = ();
    foreach my $p (@p) {
      push (@sizes, length $p);
    }
    my $string = join('', @p);
    
    if (scalar @inds == 0) {
      my @temp_inds = ();
      for (my $u = 0; $u < length $string; $u++) {
	push (@temp_inds, $u);
      }
      @inds = List::Util::shuffle(@temp_inds);
    }
    #my @temp_array = split('', $string);
    #my @shuffled = List::Util::shuffle(@temp_array);
    print SNP $info[0], "\t", $info[1], "\t", $info[2];
    my $prev = 0;
    #my $long = join('', @shuffled);
    my $long = '';
    foreach my $ind (@inds) {
      my $new_gen = substr($string, $ind, 1);
      $long .= $new_gen;
    }
    foreach my $s (@sizes) {
      my $subpop = substr($long, $prev, $s);
      #print $subpop, "\n";
      print SNP "\t", $subpop;
      $prev = $s;
    }
    print SNP "\n";
  }
  close(IN);
  close(SNP);

  # Create a temporary MAF file, and remove the temporary SNP file afterwards
  `perl /Volumes/LC_WD/NGS_Workshop/Pipeline_Scripts/snps2maf.pl temp.snps temp.maf $popstring`;
  `rm temp.snps`;

  # Filter out sites based on the user supplied cut-off values
  # remove the MAF file
  `perl /Volumes/LC_WD/NGS_Workshop/Pipeline_Scripts/filterMAF.pl temp.maf temp.filtered $mincov $minmaf $npops`;
  `rm temp.maf`;

  # Calculate Fst for every site
  # remove the temp filtered file
  `perl /Volumes/LC_WD/NGS_Workshop/Pipeline_Scripts/FST_by_MAF_Uni.pl temp.filtered temp.fst $npops`;
  `rm temp.filtered`;

  # Calculate the mean Fst for every comparison in the temp file
  `perl /Volumes/LC_WD/Monarcha_Alignment/Scripts/meanCI_Fst.pl temp.fst temp.fst_stats 0.025 .975`;
  `rm temp.fst`;

  # Process the means file and print the mean for each comparison to the out file
  open(MEANS, "temp.fst_stats") || die "\nUnable to open the file temp.fst_stats!\n";
  while (<MEANS>) {
    chomp $_;
    if ($_ =~ /^Mean/) {
      $_ =~ s/^Mean\:\s{1,}//;
      print OUT $iter, "\t", $_, "\n";
    } else {
      next;
    }
  }
  close(MEANS);
  $iter--;      
}
close(OUT);
exit;
