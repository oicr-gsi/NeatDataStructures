#!/usr/bin/perl

use warnings;
use strict;


my $outGFF = 'include_nonExons.gtf';

open(GFF, '<', $ARGV[0]) || die("Could not open GFF for reading!");
open(OUT, '>', $outGFF) || die("Could not open file for writing!");

my @column;
my @nextcolumn;
my $next;

while( <GFF> ) {
   my $pos = tell(GFF);
   chomp;
   my $line = $_;
   @column = split /\t/, $line;
   print OUT "$line\n";
   my $next = <GFF>;
   chomp $next;
   @nextcolumn = split /\t/, $next;
   if ( $nextcolumn[1] != $column[2] ) {
      # my $pos = tell();
      print OUT "$column[0]\t", $column[2] + 1, "\t", $nextcolumn[1] - 1, "\t\.\t\.\t.\tnon-exon region\n";
      seek GFF, $pos, 0;
   }
   else {seek GFF, $pos, 0;}
}
