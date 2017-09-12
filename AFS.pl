#!/usr/bin/perl -w
#
# AFS_derived.pl
#
# Find the number of SNPs in each MAF bin for each the 3 populations, to plot as Allele Frequency Spectrum
#
# August 14, 2013
# Liz Cooper

use strict;
use Data::Dumper;

# Take the names of the input and output files from the command line
my ($USAGE) = "\n$0 <input.maf> <output_afs.bins>
\tinput.maf = A text file with minor allele frequencies, created by snps2maf or getDerived
\toutput.afs.bins = A text file with the allele counts for each freq. bin in each pop.\n\n";
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Open and get the first line of the input file to get the pop names
open (IN, $input) || die "\nUnable to open the file $input!\n";
my $head = <IN>;
chomp $head;
my @header = split(/\s{1,}/, $head);
my @popColumns = ();

for (my $p=1; $p < scalar @header; $p++) {
  if (($header[$p] =~ /_MAF/) || ($header[$p] =~ /_DER/)) {
    push (@popColumns, $p)
  }
}
my $first = $popColumns[0];
my $last = $popColumns[(scalar @popColumns - 1)];

my $empty = '0' x scalar @popColumns;
my @bins = ('0','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95');
my %afs = ();
foreach my $b (@bins) {
  @{$afs{$b}} = split('', $empty);
}

# Open the output files
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Read through the input file one line at a time, and bin each SNP
while (<IN>) {
  chomp $_;
  my @info = split(/\t/, $_);
  my @maf_values = @info[$first..$last];
  foreach my $b (@bins) {
    my @temp = @{$afs{$b}};
    for (my $p = 0; $p < scalar @maf_values; $p++) {
      if (($maf_values[$p] >= $b) && ($maf_values[$p] < ($b + 0.05))) {
	$temp[$p] += 1;
      }
    }
    @{$afs{$b}} = @temp;
  }
}
close(IN);

# Print the final counts for each bin to the output file
print OUT "BIN";
foreach my $p (@popColumns) {
  my $name = $header[$p];
  $name =~ s/_(MAF|DER)//g;
  print OUT "\t", $name;
}
print OUT "\n";
foreach my $b (@bins) {
  print OUT $b, '-', ($b+0.05);
  my @counts = @{$afs{$b}};
  foreach my $c (@counts) {
    print OUT "\t", $c;
  }
  print OUT "\n";
}
close(OUT);
exit;
