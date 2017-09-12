#!/usr/bin/perl -w
#
# filterMAF.pl
#
# Filter the MAF file based on per individual coverage
# and minimum Minor Allele Frequency
#
# September 26, 2013
# Liz Cooper

use strict;

# Get as input the MAF file name, the output file name, the min number individuals and the min allele frequency
my ($USAGE) = "\n$0 <input.maf> <output.filtered> <min_coverage> <min_alleleFreq> <npops>
\tinput.maf = The input file with maf values, created by snps2maf
\toutput.filtered = An output file with maf values for sites that pass filter threshold
\tmin_coverage = The minimum number of non-missing individuals per pop to keep a site
\tmin_alleleFreq = The minimum minor allele frequency to keep a site
\tnpops = The number of populations in the file\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $mincov, $minfreq, $npops) = @ARGV;

# Open the input and output files for printing
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

FILELOOP: while (<IN>) {
  chomp $_;
  
  # Skip the header 
  if ($_ =~ /^CHROM/) {
    print OUT $_, "\n";
    next FILELOOP;
  }
  
  my @info = split(/\t/, $_);

  # Check the Minor Allele Frequency
  my $maf = pop @info;
  if ($maf < $minfreq) {
    next FILELOOP;
  }

  # Check the coverages
  for (my $i = 0; $i < $npops; $i++) {
    my $cov = $info[$i+3+$npops];
    if (($cov/2) < $mincov) {
      next FILELOOP;
    }
  }

  print OUT $_, "\n";
}
close(IN);
close(OUT);
exit;
  
