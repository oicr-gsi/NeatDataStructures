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

# creating trinucleotide context data hash
my %trinucleotide_context_data;
my %context_tally_across_mutated_to;   


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
   


   ###############################
   # new section: forming the data structure
   ###############################

   # to create N_N contexts for data structure, context_code is defined as the trinucleotide context with a blank middle allele
   my $context_code=$context;
   $context_code =~ s/\n|\r//;
   substr($context_code,1,1) = "_";
   
   # create variables for mutated_from and  mutated_to nucleotides
   my $mutated_from = $line[9];
   my $mutated_to = $line[10];

   # context_codes are totalled
   $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to} = $trinucleotide_context_data{$context_code}{$mutated_from}{$mutated_to} + 1; 
   $context_tally_across_mutated_to{$context_code}{$mutated_from} = $context_tally_across_mutated_to{$context_code}{$mutated_from} + 1; 



   # to keep track of progress
   unless ($line_count%10000) {
      print "processed $line_count lines\n";
   }
   $line_count++; 
} 
# end working through the input file
   
# print trinucleotide contexts and corresponding totals
#foreach my $context_code (keys %trinucleotide_context_data) {
 #  print "$context_code $trinucleotide_context_data{$context_code}\n";
  # open(OutputTrinucleotideContext, '>', $ARGV[3]) || die("Could not open file!");
#}


# print trinucleotide contexts and corresponding totals for every mutated_to nucleotide
foreach my $context_code (keys %trinucleotide_context_data) {

   my $outfilename = $context_code."SNPstats.matrix";
   
   open(trinucleotide_file_handle, '>', $outfilename) || die("Could not open file!");
   
   my $context_sum = 0;

   foreach my $mutated_from_nucl_key (keys %{ $trinucleotide_context_data{$context_code} }) {
   foreach my $mutated_to_nucl_key (keys %{ $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key} }) {
      print "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}\n";
      $context_sum = $context_sum + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
      }
   print trinucleotide_file_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key} / $context_tally_across_mutated_to{$context_code}{$mutated_from_nucl_key}\n";
   }
   print "\n $context_code $context_sum \n\n";
}



#foreach my $name (sort keys %grades) {
#    foreach my $subject (keys %{ $grades{$name} }) {
#        print "$name, $subject: $grades{$name}{$subject}\n";
#    }
#}




