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
   


# define the output file name for InDels and open it for writing
   my $indel_prob_file_name = "other_small_substitutions";
   open(my $indel_prob_handle, '>>', $indel_prob_file_name) || die("Could not open file!");


# define nucleotide array
my @nucleotides = ("A", "T", "C", "G");
   
foreach my $nt1 (@nucleotides) {
   foreach my $nt3 (@nucleotides) {


      # define the output file name and open it for writing
      my $trinucleotide_SNP_probability_file_name = "trinucleotide_probability.".$nt1."_".$nt3;
      open(my $trinuc_prob_handle, '>', $trinucleotide_SNP_probability_file_name) || die("Could not open file!");

      # define the output file name for InDels and open it for writing
      # my $indel_prob_file_name = "other_small_substitutions";
      # open(my $indel_prob_handle, '>', $indel_prob_file_name) || die("Could not open file!");


      # print trinucleotide contexts and corresponding totals for every mutated_to nucleotide
      #foreach my $context_code (keys %trinucleotide_context_data) {
      my $context_code=$nt1."_".$nt3;

         foreach my $mutated_from_nucl_key (keys %{ $trinucleotide_context_data{$context_code} }) {
      
            # define the "mutated_to" keys in trinuc context hash
            my $mutated_to_nucl_key;

            # the sum is only across mutated_to, and will be redefined for each mutated_from
            my $context_sum_across_mutated_to = 0;
            my $context_sum_across_indel = 0;
            print "\nRaw counts for mutated_from $mutated_from_nucl_key \n";


            foreach $mutated_to_nucl_key (keys %{ $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key} }) {
               my $mutated_from_length = length( $mutated_from_nucl_key );
               my $mutated_to_length = length( $mutated_to_nucl_key );
               if ( $mutated_from_length == 1 ) {
                  if ( $mutated_from_nucl_key ne "-" ) {
                     if ( $mutated_to_length == 1 ) {
                        if ( $mutated_to_nucl_key ne "-" ) {
                           # print "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key -- $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}\n";
                           $context_sum_across_mutated_to = $context_sum_across_mutated_to + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
                        }# end if statement

                        else {
                           $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
                        }# end else statement

                     }# end if statement

                     else {
                        $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
                     }# end else statement

                  }# end if statement

                  else {
                     $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
                     }# end else statement

               }# end if statement

               else {
                  $context_sum_across_indel = $context_sum_across_indel + $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key};
               }# end else statement

               print "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key -- $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}\n";
            }# end of loop over mutated_to
            print "\nProbabilities for mutated_from $mutated_from_nucl_key:\n";


            foreach $mutated_to_nucl_key (keys %{ $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key} }) {
               my $mutated_from_length = length( $mutated_from_nucl_key );
               my $mutated_to_length = length( $mutated_to_nucl_key );
               if ( $mutated_from_length == 1 ) {
                  if ( $mutated_from_nucl_key ne "-" ) {
                     if ( $mutated_to_length == 1 ) {
                        if ( $mutated_to_nucl_key ne "-" ) {
                           my $SNP_probability = $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}/$context_sum_across_mutated_to;
                           print "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_mutated_to=$context_sum_across_mutated_to -- $SNP_probability\n";
                           print $trinuc_prob_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_mutated_to=$context_sum_across_mutated_to -- $SNP_probability\n";
                        }# end of if statement

                        else {
                           my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}/$context_sum_across_indel;
                           print $indel_prob_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                        }# end else statement

                     }# end of if statement

                     else {
                        my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}/$context_sum_across_indel;
                        print $indel_prob_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
               }# end else statement

                  }# end of if statement

                  else {
                     my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}/$context_sum_across_indel;
                     print $indel_prob_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
                  }# end else statement

               }# end of if statement
              
               else {
                  my $indel_probability = $trinucleotide_context_data{$context_code}{$mutated_from_nucl_key}{$mutated_to_nucl_key}/$context_sum_across_indel;
                  print $indel_prob_handle "$context_code, $mutated_from_nucl_key, $mutated_to_nucl_key, context_sum_across_indel=$context_sum_across_indel -- $indel_probability\n";
               }# end else statement
 
            }# end of loop over mutated_to

         }# end of loop over mutated_from
     #}# close the foreach context_code loop

     print "\n\n";
     


  }# end loop over nt3
}# end loop over nt1





