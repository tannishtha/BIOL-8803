#!/usr/bin/perl -w
use strict;
use warnings;
use 5.010;
use Getopt::Long qw(GetOptions);
use List::Util qw(min max);

#variables to store the options sent from the command line
my $file1;
my $file2;
my $percent; #percent overlap
my $join=0;
my $output;
my $condition; # condition to allow overlap only on equal strand
GetOptions('i1=s' => \$file1,
	'i2=s' => \$file2,
	'm=s' => \$percent,
	'j' => \$join,
	'k:s' => \$condition, #: denotes optional option 
	'o=s' => \$output,
) or die "Usage: $0 -i1 <file1> -i2 <file2> -m <percent Overlap> -j -k <condition> -o <output file> {j and k are optional}\n";

open FILE1, "$file1" or die "Cannot open $file1\n";
open FILE2, "$file2" or die "Cannot open $file2\n";
open OUTPUT, ">$output" or die "Cannot open $output\n";
my $key;
my %file1Start; #hash storing start coordinates of file1
my %file1Stop;  #hash storing stop coordinates of file1
my %file2Start; #hash storing start coordinates of file2
my %file2Stop;  #hash storing stop coordinates of file2
my %file1Strand; #hash to store the strand of coordinates of file1
my %file2Strand; #hash to store the strand of coordinates of file1
my $i;
my $j;
my $overlap;
my $overlapPercent;

#`sort -c -V -k1,1 -k2n,2 $file1`; #checks if file1 is sorted or not
#`sort -c -V -k1,1 -k2n,2 $file2`; #checks if file2 is sorted or not

#Adding elements to file1Start, file1Stop and file1Strand
while ( my $line = <FILE1> ) {
	my @element = split(/\s+/, $line);
	chomp (@element);
	if ( ( exists ( $file1Start{$element[0]} ) ) eq "" ) {
		$file1Start{$element[0]} = [];
	}
	push ( @{ $file1Start{$element[0]} }, $element[1]);
	 
	if ( ( exists ( $file1Stop{$element[0]} ) ) eq "" ){
                $file1Stop{$element[0]} = [];
	}
    push ( @{ $file1Stop{$element[0]} } , $element[2] );
    
    if ($#element>2){
    	if ( ( exists ( $file1Strand{$element[0]} ) ) eq "" ){
                $file1Strand{$element[0]} = [];
		}
    	push ( @{ $file1Strand{$element[0]} } , $element[3] );
    }

}

#Adding elements to file2Start, file2Stop and file2Strand
while ( my $line = <FILE2> ) {
        my @element = split(/\s+/, $line);
        chomp (@element);
        if ( ( exists ( $file2Start{$element[0]} ) ) eq "" ) {
                $file2Start{$element[0]} = [];
        }
        push ( @{ $file2Start{$element[0]} }, $element[1]);

        if ( ( exists ( $file2Stop{$element[0]} ) ) eq "" ){
                $file2Stop{$element[0]} = [];
        }
        push ( @{ $file2Stop{$element[0]} } , $element[2] );
        
        if ($#element>2){
    		if ( ( exists ( $file2Strand{$element[0]} ) ) eq "" ){
                $file2Strand{$element[0]} = [];
			}
    		push ( @{ $file2Strand{$element[0]} } , $element[3] );
    	}
}

close FILE1;
close FILE2;

#subroutine to calculate the overlap 
sub overlap{
	my ( $i, $j, $key) = @_;
	$overlap = min($file1Stop{$key}[$i] , $file2Stop{$key}[$j]) - max($file1Start{$key}[$i] , $file2Start{$key}[$j]);
	$overlapPercent =($percent/100) * ($file1Stop{$key}[$i] - $file1Start{$key}[$i]);
	if ($condition){
		if ($file1Strand{$key}[$i] eq $file2Strand{$key}[$j] ){ 	
			if ($overlap >= $overlapPercent){
				print OUTPUT "\n$key\t $file1Start{$key}[$i]\t $file1Stop{$key}[$i]\t $file1Strand{$key}[$i]";
				if ( $join) {
					print OUTPUT "\t$key\t $file2Start{$key}[$j]\t $file2Stop{$key}[$j]\t $file2Strand{$key}[$j]";
				}
			}
		}
	}
	else {
		if ($overlap >= $overlapPercent){
			print OUTPUT "\n$key\t $file1Start{$key}[$i]\t $file1Stop{$key}[$i]";
			if ( $join) {
				print OUTPUT "\t$key\t $file2Start{$key}[$j]\t $file2Stop{$key}[$j]";
			}
		}
	}	
}

foreach $key (sort keys %file1Start){
	my $file1Len =$#{$file1Start{$key}};
	my $file2Len =$#{$file2Start{$key}};
	for ($i=0, $j=0; ($file1Len < $file2Len ) ? $i<$file1Len && $j<$file2Len : $j<$file2Len && $i<$file1Len ; ){
		if ( $file2Start{$key}[$j] < $file1Start{$key}[$i] ) {
			if ( $file2Stop{$key}[$j] > $file1Start{$key}[$i] ) {
				 overlap($i,$j,$key);
			}	
			else{
				$j++;
				next;
			}
		}
		else{
			if ($file2Start{$key}[$j] < $file1Stop{$key}[$i] ){
				overlap($i, $j, $key);

			}
			else{
				$i++;
				next;
			}
		}
		for (my $k = $j+1; $k < $file2Len; $k++){
			if ($file2Start{$key}[$k] < $file1Stop{$key}[$i]){
				overlap($i, $k, $key);
				
           	}
			else{
				goto NEXT;
			}
		}
		NEXT: 
		$i++;
	}
}
