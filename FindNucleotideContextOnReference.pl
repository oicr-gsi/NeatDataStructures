#!/usr/bin/perl

use strict;


if ($#ARGV < 1) {
   print "parameter mismatch\nTo run type this command:\nperl $0 fastahack reference input_pos_file output_file\n\n";

   print " first argument = full path to fastahack\n"; 
   print " second argument = full path to reference genome\n"; 
   print " third argument = input file with arbitrary number of columns, but 1st col=chromosome name and 2nd col=position\n"; 
   print " fourth argument = output file with three columns: chromosome name, position of the center nucleotide, and the thre-nucleotide context for that position\n\n\n"; 
   exit 1;
}


my $Fastahack=$ARGV[0];
my $Reference=$ARGV[1];
open(InputPositions,             '<', $ARGV[2]) || die("Could not open file!");
open(OutputTrinucleotideContext, '>', $ARGV[3]) || die("Could not open file!");




################ read in one coordinate at a time and execute fastahack on it

# reading the header
my $head = <InputPositions>;
$head =~ s/\n|\r//;
print OutputTrinucleotideContext "$head\tContext\n";


# reading the positional information
my $line_count = 1;
while (<InputPositions>) {
   $_ =~ s/\n|\r//;
   #print "$_\n";
   my @line = split('\t', $_);

   # getting the chromosome and coordinate fields from input file
   # fastahack will need to the chromosome and coordinate to read the information from the reference
   my $chromosome = $line[0];
   my $coordinate = $line[1];

   # get coordinates of first and last character in the context
   my $start_region = $coordinate - 1;
   my $end_region = $coordinate + 1;

   # if the coordinate is the very first letter on the chromosome, then do not read before that position
   # the context becomes 2 letter code, as opposed to a trinucleotide
   if ( $start_region == 0 ) {
      $start_region = 1;
      $end_region = 2;
   }

   #print "$Fastahack -r $chromosome:$start_region..$end_region $Reference\n";
   my $context = `$Fastahack -r $chromosome:$start_region..$end_region $Reference`;
   
   # capitalize context letters
   $context = uc($context);

   # split germline column into germline allele and mutated_to allele
   my @germline = split ('/', $line[6]);

   # if germline allele does not equal reference allele, print "start_region germline allele end_region"
   # specifically, replace the middle letter of the context with the germline allele
   #print "$germline[0], $germline[1]\n";
   if ($germline[0] ne $germline[1]) {
      print "germline/reference mismatch, line number $line_count\n";
      if ($coordinate != 1) {
         substr($context,1,1)= $germline[1];   
      }
      else {
         substr($context,0,1)= $germline[1];
      }
   }

   print OutputTrinucleotideContext "$_\t$context";
   

   # to keep track of progress
   unless ($line_count%10000) {
      print "processed $line_count lines\n";
   }
   $line_count++; 
}
