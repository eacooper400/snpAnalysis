#!/usr/bin/perl -w
#
# SNP_Categories_uni.pl
#
# A more robust, universal script for finding shared and private polymorphisms
# Uses the maf output file as input
#
# October 17, 2013
# Liz Cooper

use strict;
#use Math::Combinatorics;
use Data::Dumper;
#use Algorithm::Loops qw(NestedLoops);

# Take as input the name of the maf file and the name of 2 output files
# Output one is for the total counts in each category
# output 2 is for the full snp list with coverage a category info
my $USAGE = "\n$0 <input.maf> <output.cats> <output.counts>
\tinput.maf = The input file with minor allele frequencies
\toutput.cats = A text file with the total counts for each SNP category
\toutput.counts = The full SNP list with a coverage category\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $out1, $out2) = @ARGV;

# Open the input file and process the header line
# Use the pop names in the header to determine the sharing categories
open (IN, $input) || die "\nUnable to open the file $input!\n";
my $header = <IN>;
chomp $header;
my @fields = split(/\t/, $header);
my @pops = ();
foreach my $field (@fields) {
  if ($field =~ /Alleles/) {
    my $pop = ((split(/_/, $field))[0]);
    push (@pops, $pop);
  }
}
my $npops = scalar @pops;
my %pop_codes = ();
for (my $p = 0; $p < scalar @pops; $p++) {
  $pop_codes{$p} = $pops[$p];
}
my %categories = ();
my %cat_codes = ();
my $start = 1;
for my $comb (1..2**@pops-1) {
  my $s = join ':',
  map $pops[$_],
  grep $comb & (1<<$_),
  0..$#pops;
  $categories{$s} = 0;
  $cat_codes{$s} = $start;
  $start++;
}
 

# Open the ouput file for the snp site by category output
open (CATS, ">$out1") || die "\nUnable to open the file $out1!\n";

# Now read through the rest of the input file, and categorize each SNP
while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  my @mafs = @info[3..$npops+2];
  my @covs = @info[$npops+3..($npops+$npops+2)];
  my $code_string = '';
  for (my $m = 0; $m < scalar @mafs; $m++) {
    if ($mafs[$m] > 0) {
      my $name = $pop_codes{$m};
      $code_string .= $name . ':';
    }
  }
  chop $code_string;
  $categories{$code_string} += 1;
  my $ccode = $cat_codes{$code_string};
  #if ($ccode) {
    print CATS $ccode;
  #} else {
  #  exit;
  #}
  foreach my $cov (@covs) {
    if ($cov == 0) {
      print CATS "\t", $cov;
    } else {
      print CATS "\t", $cov/2;
    }
  }
  print CATS "\n";
}
close(IN);
close(CATS);

# Finally print the counts for each category to the Counts out file
open (COUNTS, ">$out2") || die "\nUnable to open the file $out2!\n";
my @cat_numbers = values %cat_codes;
my @sorted_cats = sort {$a <=> $b} @cat_numbers;
foreach my $cat (@sorted_cats) {
  foreach my $key (keys %cat_codes) {
    if ($cat_codes{$key} == $cat) {
      print COUNTS $cat, "\t", $key, "\t", $categories{$key}, "\n";
    }
  }
}
close(COUNTS);
exit;


#  my $combinat = Math::Combinatorics->new(count => $p+1,
#					data => [@pops],#
					# );
#  while(my @combo = $combinat->next_combination){
#    my $cat = join(' ', @combo);
#    $categories{$cat} = 0;
#    $cat_codes{$cat} = $start;
#    $start++;
#  }
#}
