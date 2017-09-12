#!/usr/bin/perl -w
#
# getSitelist.pl
#
# Get a list of SNP sites to use in making an IM or Migrate input file
#
# October 2, 2013

use strict;

# Get the name of the filtered snps file and the output file
my ($USAGE) = "\n$0 <input.snps> <output.list>
\tinput.snps = The input file in .snps format
\toutput.list = A filtered list of sites that meet dadi or IM requirements\n\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}
my ($input, $output) = @ARGV;

# Open the input and output files
open (IN, $input) || die "\nUnable to open the file $input!\n";
open (SITES, ">$output") || die "\nUnable to open the file $output!\n";

LOOP: while (<IN>) {
  chomp $_;
  my @info = split(/\s{1,}/, $_);
  my $id = $info[0];
  my $p = $info[1];
  my $q = $info[2];
  my @pops = @info[3..((scalar @info) -1 )]; 

  # Skip indels
  foreach my $string (@pops) {
    if ($string =~ /(0|1)/) {
      next LOOP;
    }
  }

  # Skip sites with missing data
  my $temp_string = join('', @pops);
  if ($temp_string =~ /N/) {
    next LOOP;
  }

  # Skip sites with overall MAF less than 10%
  $temp_string =~ s/A/AA/ig;
  $temp_string =~ s/C/CC/ig;
  $temp_string =~ s/G/GG/ig;
  $temp_string =~ s/T/TT/ig;
  
  $temp_string =~ s/M/AC/ig;
  $temp_string =~ s/R/AG/ig;
  $temp_string =~ s/W/AT/ig;
  $temp_string =~ s/S/CG/ig;
  $temp_string =~ s/Y/CT/ig;
  $temp_string =~ s/K/GT/ig;

  my $mac = 0;
  for (my $k = 0; $k < length $temp_string; $k++) {
    if (substr($temp_string, $k, 1) =~ /$q/) {
      $mac++;
    }
  }
  my $maf = $mac/(length $temp_string);
  if ($maf < 0.1) {
    next LOOP;
  }

  # Check that all (required) individuals are present
  #my @new_m = ();
  #my @m =split('', $mstring);
  #for (my $i = 0; $i < scalar @m; $i++) {
  #  unless (($i == 14) || ($i == 15) || ($i == 16)) {
  #    push (@new_m, $m[$i]);
  #  }
  #}
  #my $new_mstring = join('', @new_m);
  #my $new_mstring_temp = $new_mstring;
  #$new_mstring_temp =~ s/N//ig;
  #if (length $new_mstring_temp < 15) {
  #if ($new_mstring =~ /N/) {
  #  next LOOP;
  #}
  
  #my $sstring_temp = $sstring;
  #$sstring_temp =~ s/N//ig;
  #if (length $sstring_temp < 15) {
  #if ($sstring =~ /N/) {
  #  next LOOP;
  #}
  
  #my @new_u = ();
  #my @u =split('', $ustring);
  #for (my $i = 0; $i < scalar @u; $i++) {
   # unless (($i == 17) || ($i == 20) || ($i == 23)) {
   #   push (@new_u, $u[$i]);
  #  }
  #}
  #my $new_ustring = join('', @new_u);
  #my $new_ustring_temp = $new_ustring;
  #$new_ustring_temp =~ s/N//ig;
  #if (length $new_ustring_temp < 15) {
  #if ($new_ustring =~ /N/) {
  #  next LOOP;
  #}

  #my @new_r = ();
  #my @r =split('', $rstring);
  #for (my $i = 0; $i < scalar @r; $i++) {
  #  unless ($i == 19) {
   #   push(@new_r, $r[$i]);
 #   }
 # }
  #my $new_rstring = join('', @new_r);
  #my $new_rstring_temp = $new_rstring;
  #$new_rstring_temp =~ s/N//ig;
  #if (length $new_rstring_temp < 15) {
  #if ($new_rstring =~ /N/) {
   # next LOOP;
  #}

  # Print the new sites and the new genotypes to a file
  my ($tag, $position) = split(/\./, $id);
  print SITES $tag, "\t", $position;
  foreach my $pop (@pops) {
    print SITES "\t", $pop;
  }
  print SITES "\n";
}
close(IN);
close(SITES);
exit;
