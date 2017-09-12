#!/usr/bin/perl -w
#
# FST_by_MAF_Uni.pl
#
# A more universal script to calculate Fst based on minor allele frequencies
# Should work on any MAF file regardless of the number of populations
#
# October 3, 2013

use strict;

# open the input maf file and the output fst file
my ($USAGE) = "\n$0 <infile.maf> <outfile.fst> <npops>
\tinfile.maf = The input file with minor allele frequencies; created by snps2maf
\toutfile.fst = The output file with Fst calculated for every pair of populations
\tnpops = The number of populations in the file\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $npops) = @ARGV;

open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Print a header line to the output file
print OUT "Locus";
for (my $n = 1; $n < $npops; $n++) {
  for (my $p = $n+1; $p <= $npops; $p++) {
    print OUT "\t", "Fst", $n, "_", $p;
  }
}
print OUT "\n";

# Process the input file one line at a time
while (<IN>) {
  chomp $_;
  if ($_ =~ /^CHROM/) {
    next;
  }
  my @info = split(/\t/, $_);

  my @exp = ();
  my @num = ();
  for (my $i = 0; $i < $npops; $i++) {
    my $maf = $info[$i+3];
    my $num1 = $info[$i+3+$npops];
    push (@num, $num1);
    my $exp1 = 1 - (((1-$maf)**2) + ($maf**2));
    push (@exp, $exp1);
  }

  my @Sub = ();
  my @qbar = ();
  my @Tot = ();
  for (my $i = 0; $i < $npops-1; $i++) {
    for (my $j = $i+1; $j < $npops; $j++) {
      my $S = 0;
      my $q = 0;
      unless (($num[$i] + $num[$j]) == 0) {
	$S = (($exp[$i] * $num[$i]) + ($exp[$j] * $num[$j]))/($num[$i] + $num[$j]);
	$q = (($info[$i+3] * $num[$i]) + ($info[$j+3] * $num[$j]))/($num[$i] + $num[$j]);
      }
      my $t = 1 - (((1 - $q)**2) + ($q**2));
      push (@Sub, $S);
      push (@qbar, $q);
      push (@Tot, $t);
    }
  }

  print OUT $info[0];

  for (my $y = 0; $y < scalar @Sub; $y++) {
    my $fst = 0;
    unless ($Tot[$y] == 0) {
      $fst = ($Tot[$y] - $Sub[$y])/$Tot[$y];
    }
    print OUT "\t", $fst;
  }
  print OUT "\n";
}
close(IN);
close(OUT);
exit;
