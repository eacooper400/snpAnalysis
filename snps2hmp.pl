#!/usr/bin/perl -w
#
# snps2hmp.pl
#
# Convert snps format to hapmap format used by TASSEL
#
# May 29, 2014
# Liz Cooper

use strict;

# From the command line, get the names of the input and output files,
# a file with the individual names,
# the minimum population coverage
# and the minimum minor allele frequency
my ($USAGE) = "\n$0 <input.snps> <output.hmp> <names.txt> <min_cov> <min_MAF> <Haploid or Diploid>
\tinput.snps = The input file in .snps format
\toutput.hmp = The name of a hapmap file (used by TASSEL) to create
\tnames.txt = A file with sample names.
\t\tOne line per POPULATION, with a comma delimited list of sample names in each line
\tmin_cov = The minimum number of non-missing individuals to keep a site
\tmin_MAF = The minimum minor allele frequency to keep a site
\tHaploid or Diploid = Whether to use Heterozygot codes (diploid) or print 2 haploid strings for each individual\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $namefile, $minCov, $minMAF, $ploidy) = @ARGV;

# Open the input and output files
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Print the header line to the output, then get the names from the names file
print OUT "rs\t", "alleles\t", "chrom\t", "pos\t", "strand\t", "assembly\t", "center\t", "protLSID\t", "assayLSID\t", "panelLSID\t", "QCode";
open (NAMES, $namefile) || die "\nUnable to open the file $namefile!\n";
while (<NAMES>) {
  chomp $_;
  my @temp = split(/,/, $_);
  foreach my $n (@temp) {
    if ($ploidy =~ /H/) {
      print OUT "\t", $n, "_1", "\t", $n, "_2";
    } else {
      print OUT "\t", $n;
    }
  }
}
close(NAMES);
print OUT "\n";

# Save all of the positions for each chromosome outside of the loop, so that they can be sorted before printing
my %sites = ();
my $prev_chr = 0;

# Start reading through the input file, and convert each line to the hapmap format
LOOP:while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  my $rs = $info[0];
  $rs =~ s/^JH//g;
  $rs =~ s/\./_/g;
  my $chr = 0;
  my $pos = 0;
  if ($info[0] =~ /\./) {
    ($chr, $pos) = split(/\./, $info[0]);
  } else {
    $chr = $info[0];
    $pos = $info[1];
  }
  $chr =~ s/^JH//g;
  my @pops = ();
  if ($info[0] =~ /\./) {
    @pops = @info[3..(scalar @info - 1)];
  } else {
    @pops = @info[5..(scalar @info -1)];
  }
  # Check that each pop meets the minimum coverage requirements
  foreach my $pop (@pops) {
    my $temp_pop = $pop;
    $temp_pop =~ s/N//ig;
    if (length $temp_pop < $minCov) {
      next LOOP;
    }
  }
  my $string = join('', @pops);

  # Check that the SNP sites meets the minimum MAF requirement
  my $temp = $string;
  $temp =~ s/N//ig;
  my $total = length $temp;
  if ($info[0] =~ /\./) {
    $temp =~ s/$info[1]//ig;
  } else {
    $temp =~ s/$info[3]//ig;
  }
  my $count = length $temp;
  my $check = $count/$total;
  if (($check < $minMAF) || ($check > (1-$minMAF))) {
    next LOOP;
  }

  # If the site passes both filters, create an output string for it
  my $outstring = '';
  if ($info[0] =~ /\./) {
    $outstring = $rs . "\t" . $info[1] . "/" . $info[2] . "\t" . $chr . "\t" . $pos . "\t" . "+" . "\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "NA";
  } else {
    $outstring = $rs . "\t" . $info[3] . "/" . $info[4] . "\t" . $chr . "\t" . $pos . "\t" . "+" . "\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "NA";
  }
  my @genotypes = split('', $string);
  foreach my $gen (@genotypes) {
    if ($ploidy =~ /H/) {
      my $string1 = $gen;
      my $string2 = $gen;
      $string1 =~ tr/MRWSYK/AAACCG/;
      $string2 =~ tr/MRWSYK/CGTGTT/;
      $outstring .=  "\t" . $string1 . "\t" . $string2;
    } else {
      $outstring .= "\t" . $gen;
    }
  }

  # Check if this marks a the start of a new chromosome
  # If so, then sort the positions for the previous chromosome and print them to the outfile
  if ($chr == $prev_chr) {
    $sites{$pos} = $outstring;
  } elsif ($prev_chr == 0) {
    $sites{$pos} = $outstring;
    $prev_chr = $chr;
  } else {
    my @positions = sort {$a <=> $b} (keys %sites);
    foreach my $ordered (@positions) {
      print OUT $sites{$ordered}, "\n";
    }
    %sites = ();
    $sites{$pos} = $outstring;
    $prev_chr = $chr;
  }
}

# Print the information for the final chromosome
my @positions = sort {$a <=> $b} (keys %sites);
foreach my $ordered (@positions) {
  print OUT $sites{$ordered}, "\n";
}

close(IN);
close(OUT);


exit;
