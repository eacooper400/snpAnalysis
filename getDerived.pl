#!/usr/bin/perl -w
#
# getDerived_alleles.pl
#
# Find derived and ancestral alleles using Russells as outgroup
#
# September 30, 2013

use strict;

# Open the input and output files
my ($USAGE) = "\n$0 <input.filtered> <output.derived> <OUTGROUP_NAME>
\tinput.filtered = An input text file with minor allele frequences (.maf)
\toutput.derived = An output file with ancestral and derived frequencies
\tOUTGROUP_NAME = The name of the population to use as an outgroup\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $outgroup) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

my $head = <IN>;
chomp $head;
my @header = split(/\s{1,}/, $head);
my @popCol = ();
my $out;

for (my $p = 0; $p < scalar @header; $p++) {
  if ($header[$p] =~ /_MAF/) {
    if ($header[$p] =~ /$outgroup/) {
      $out = $p;
    } else {
      push (@popCol, $p);
    }
  }
}

print OUT "CHROM", "." , "POS", "\tANC\tDER";
foreach my $p (@popCol) {
  my $name = $header[$p];
  $name =~ s/MAF/DER/g;
  print OUT "\t", $name;
}
print OUT "\n";

while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  my $ancs = '';
  my $der = '';
  
  if ($info[$out] == 0) {
    $ancs = $info[1];
    $der = $info[2];
    print OUT $info[0], "\t", $ancs, "\t", $der;
    foreach my $p (@popCol) {
      print OUT "\t", $info[$p];
    }
    print OUT "\n";
  } elsif ($info[$out] == 1) {
    $ancs = $info[2];
    $der = $info[1];
    print OUT $info[0], "\t", $ancs, "\t", $der;
    foreach my $p (@popCol) {
      print OUT "\t", (1-$info[$p]);
    }
    print OUT "\n";
  } else {
    next;
  }
}
close(IN);
close(OUT);
exit;
  
