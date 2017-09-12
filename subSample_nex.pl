#!/usr/bin/perl -w
#
# subSample_nex.pl
#
# Take random sub-samples of individuals from a larger nexus file
#
# June 11, 2014
# Liz Cooper

use strict;
srand(time|$$);

# Take as command line input:
# 1.  The name of the full nexus file
# 2.  The name of the output nexus file
# 3.  The number of individuals to subsample from each population (as a comma delimited list)
my ($USAGE) = "\n$0 <input.nex> <output.nex> <List of Subsample sizes>
\tinput.nex = The name of the full input nexus file
\toutput.nex = The name of the subsampled output nexus file
\tList of Subsample Sizes = A comma delimited list of the # individuals to subsample from each population\n\n";
  
unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output, $sizelist) = @ARGV;
my $templine = $sizelist;
$templine =~ s/,/\+/g;
my @sizes = split(/,/, $sizelist);

# Open the output file for printing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";

# Create an array to save all of the individuals in each population outside of the file loop
# Use a counter to keep track of which population is being processed
my @individuals = ();
my $pop = '';
my $counter = 0;

# Open and start reading through the input file
open (IN, $input) || die "\nUnable to open the file $input!\n";

while (<IN>) {
  chomp $_;

  # Go through the header lines and print them out as they are,
  # excpet for the ntax information, which needs to be updated
  if ($_ =~ /\#NEXUS/) {
    print OUT $_, "\n";
  } elsif ($_ =~ /^Begin/) {
    print OUT $_, "\n";
  } elsif ($_ =~ /^Dimensions/) {
    my @dim = split(/\s{1,}/, $_);
    my $ntax = eval($templine);
    print OUT $dim[0], " ", "ntax=", $ntax, " ", $dim[2], "\n";
  } elsif ($_ =~ /^Format/) {
    print OUT $_, "\n";
  } elsif ($_ =~ /^Matrix/) {
    print OUT $_, "\n";
  }

  # Check if the line is the start of a new population,
  # and then save all of the individuals from that population
  else {
    my ($name, $string) = split(/\s{1,}/, $_);
    my $popID = (split(/_/, $name))[0];
    #print $popID, "\n";
    
    if (length $pop > 0) {
      if ($pop =~ /$popID/) {
	push (@individuals, $_);
      } else {
	
	# Get a random subset of individuals (without replacement)
	my $n = $sizes[$counter];
	while ($n) {
	  my $pick = int(rand(scalar @individuals - 1));
	  print OUT $individuals[$pick], "\n";
	  splice(@individuals, $pick, 1);
	  $n--;
	}
	
	# Now assign a new pop ID, and re-set the array 
	@individuals = ();
	push (@individuals, $_);
	$pop = $popID;
	$counter++;
      }
    } else {
      $pop = $popID;
      push (@individuals, $_);
    }
  }
}
close(IN);

print OUT ";", "\n";
print OUT "End;", "\n";
close(OUT);
exit;
