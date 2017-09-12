#!/usr/bin/perl -w
#
# snps2nexus.pl
#
# Create a NEXUS file (for use in BEAUTi and BEAST)
#
# May 27, 2014

use strict;

# Get the names of the input and output files from the command line
# Get the name of a names file (same format as used in snps2hmp and snps2structure
my ($USAGE) = "\n$0 <input.snps> <output.nex> <names.txt>
\tinput.snps = The input file in .snps format
\toutput.nex = The name of a nexus output file to create
\tnames.txt = A file with 1 row per population, and a comma delimited list of names on each row\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $namefile) = @ARGV;

# Save the names for each individual
open (NAMES, $namefile) || die "\nUnable to open the file $namefile!\n";
my @names = ();
while (<NAMES>) {
  chomp $_;
  push (@names, $_);
}
close(NAMES);

# Store the individual sequence strings in a hash
my %sequences = ();

# Save the number of snps as they are looped through
my $num_snps = 0;

# Open the input file and start saving the sequences for each individual
open (IN, $input) || die "\nUnable to open the file $input!\n";
while (<IN>) {
  chomp $_;
  $num_snps++;
  my @info = split(/\s{1,}/, $_);
  my ($tag, $pos) = split(/\./, $info[0]);
  my $ref = $info[1];
  my $alt = $info[2];
  my @pops = @info[3..(scalar @info - 1)];
		  
  for (my $pop = 0; $pop < scalar @names; $pop++) {
    for (my $p = 0; $p < length $pops[$pop]; $p++) {
      my $base = substr($pops[$pop], $p, 1);
      my $code = '';
      if ($base =~ /$ref/) {
	$code = 0;
      } elsif ($base =~ /[MRWSYK]/) {
	$code = 1;
      } else {
	$code = 2;
      }
          
      my $ind_name = (split(/,/, $names[$pop]))[$p];
      if (exists $sequences{$ind_name}) {
	$sequences{$ind_name} .= $code;
	#$sequences{$ind_name} .= $base;
      } else {
	$sequences{$ind_name} = $code;
	##$sequences{$ind_name} = $base;
      }
    }
  }
}
close(IN);

# Determine the total number of taxa
my $long_string = join(",", @names);
my @temp = split(/,/, $long_string);
my $num_taxa = scalar @temp;

# Open the output file and print in the nexus format
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";
print OUT "#NEXUS", "\n";
print OUT "Begin data;", "\n";
print OUT "Dimensions ntax=", $num_taxa, " nchar=", $num_snps, ";", "\n";
print OUT "Format datatype=dna", " symbols=", "ACTG",  " missing=? gap=?;", "\n";
print OUT "Matrix", "\n";

foreach my $name (@temp) {
  my $string = $sequences{$name};
  print OUT $name, "    ", $string, "\n";
}
print OUT ";", "\n";
print OUT "End;";

close(OUT);
exit;


  
  
