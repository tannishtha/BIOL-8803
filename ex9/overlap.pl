#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som
#Script to find overlaps in a BED file

my %chrHash; #stores the unique coordinates for a chromosome in a sorted manner in a hash of array 
my %inputHash; #stores the start and end coordinates of each chromosome name from the input bed file in a hash of 2d arrays 
my %outputHash; #stores the output of the program with coverage for each base in hash of 2d arrays
die "Usage : $0 <InputFile> <OutputFile>\n" if @ARGV < 2;
open FILE, "$ARGV[0]" or die "Cannot open $ARGV[0]"; #input file
my $i = 0;
my $j = 0;
my $k = 0;
my $count = 0; #stores the coverage for the coordinates of a row
my $num1 = 0;
my $num2 = 0;
my $key;
while ( my $line = <FILE> ) {
	my @element = split(/\s+/, $line);
	chomp (@element);
	#adding elements to the hash of sorted coordinates
	if ( ( exists ( $chrHash{$element[0]} ) ) eq "" ) {
		$chrHash{$element[0]} = [];
	}
	push ( @{ $chrHash{$element[0]} }, $element[1], $element[2] );
	 
	#adding elements to the hash of input coordinates
	if ( ( exists ( $inputHash{$element[0]} ) ) eq "" ){
                $i = 0;
                $inputHash{$element[0]} = [];
	}
        push ( @{ $inputHash{$element[0]}[$i] } , $element[1], $element[2] );
        $i++;
}
close FILE;

foreach $key (keys %chrHash) {
	@{ $chrHash{$key} } = sort { $a <=> $b } @{ $chrHash{$key} }; #sorting the hash of array numerically containing all the coordinates of a chromosome  
	#storing only the unique values in the hash of array
	my @result = ( $chrHash{$key}[0] );
 	my $last =  $chrHash{$key}[0];
	foreach ( @{ $chrHash{$key} }){
   		if ($_ ne $last){
         		push (@result, $_);
         		$last = $_;
   		}
	}
	@{ $chrHash{$key} } = @result;
}
#finding out the coverage for each base of the chromosomes
foreach $key (keys %chrHash){
	for $i (0..($#{$chrHash{$key}})-1){
		($num1, $num2, $count) = ($chrHash{$key}[$i], $chrHash{$key}[$i+1], 0);  #Taking two consecutive elements from the array and storing initial count of them
		for $j (0..$#{$inputHash{$key}}) {
			if ( ($num1 >= $inputHash{$key}[$j][0]) &&  ( $num2 <= $inputHash{$key}[$j][1] ) ) { #increase count if the two numbers lie in range in the coordinates for that chromosome
				$count++;
			}
		}
		if ( ( exists ( $outputHash{$key} ) ) eq "" ){
                	$k = 0;
                	$outputHash{$key} = [];
		}
                if ( $count){
			push ( @{ $outputHash{$key}[$k] } , $num1, $num2, $count );
                	$k++;
		}
	}
}

#Printing the results in the output file
open FILE1,">$ARGV[1]";
print FILE1 "Chrom\tbase1\tbase2\tcoverage\n\n";
foreach $key (sort keys %outputHash) {
        for  $j (0..$#{$outputHash{$key}}){
                print FILE1 "$key\t",join("\t", @{$outputHash{$key}[$j]}),"\n";
        }
	print FILE1 "\n";

}
close FILE1;
