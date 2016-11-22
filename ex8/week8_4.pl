#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som
#script to print the coordinates for the genes in Infectious Diseases file from knownGene and kgXref

open FILE1, "$ARGV[0]" or die "Could not open $ARGV[0]";  #knownGene file
open FILE2, "$ARGV[1]" or die "Could not open $ARGV[1]";  #kgXref file
open FILE3, "$ARGV[2]" or die "Could not open $ARGV[2]";  #Infectious disease file

<FILE3>; #read and discard the first line of the Infectious disease file as it contains a header

print "Gene\t\tChromosome\t\tStrand\t\tTx Start\t\tTx End\n"; 
#open the Infectious disease file from the second line
while ( my $line1 = <FILE3> ){
	$line1 =~ s/\r?\n//; #takes care of lines ending with \r\n or \n (chomp not working)
	seek FILE2, 0, 0; #bring the filehandler for kgXref to the beginning after each search
	my $id; #store the known Gene ID
	while ( my $line2 = <FILE2> ){
		my @element1 = split (/\t/, $line2 );
		if ( $element1[4] eq $line1 ){ #finds the IDs of the genes present in the Infectious Disease file from kgXref file
			$id = $element1[0];		
			seek FILE1, 0, 0; #bring the filehandler for knownGene to the beginning after each search
        		while ( my $line3 = <FILE1> ) {
                		my @element2 = split (/\t/, $line3 ); #gets the coordinates for the genes from knownGene file from the ID found
				if ( $id eq $element2[0]){
                                	printf "%-7s\t\t%-7s\t\t%12s\t\t%-7s\t\t%-7s\n", $line1, $element2[1], $element2[2], $element2[3], $element2[4]; #chrom,strand,txStart and txEnd
                        	}

               		}
	
		}
	}
	print "\n";
	
}

close FILE1;
close FILE2;
close FILE3;
