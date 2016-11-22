#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som
#Script to summarise a BED file

my $countEntries=0; #variable to count the total number of Entries in the file
my $totalLength=0; #variable to count the total length of the entries in the file
my $plusEntries=0; #variable to count the number of entries on the + strand
my $minusEntries=0; #variable to count the number of entries on the - strand
my $longest=0; #length of the longest entry
my $longestEntry; #name of the longest entry
my $shortest=10000; #length of the shortest entry
my $shortestEntry; #name of the shortest entry
my $length=0; #length of each line
my $avLength=0; #average length 
my $sdDeviation=0; #standard deviation of the length 
my $sum=0; #variable to calculate the summation of the squared differences between length and average length
my @arrayLength; #array to store the length of each entry to calculate standard deviation

open FILE, "$ARGV[0]" or die "Cannot open '$ARGV[0]'";

while ( my $line = <FILE> ){
	$countEntries++; #part (a) of the question
	my @elements = split(/\t/, $line);
	chomp(@elements);
	$length = $elements[2] - $elements[1] + 1;
	push (@arrayLength, $length);
	$totalLength += $length; #part (b) of the question

	#part (c) of the question
	if ( $elements[5] eq "+" ){
		$plusEntries++;
	}
	else{
		$minusEntries++;
	}
	
	#part (d) of the question
	if ( $length > $longest){
		$longest = $length;
		$longestEntry = $elements[3];
	}

	#part (e) of the question	
	if ( $length < $shortest){
		$shortest = $length;
		$shortestEntry = $elements[3];
	}
}
close FILE;

#part (f) of the question
$avLength = $totalLength / $countEntries; #calculating average of the length

#calculating standard deviation of the length
foreach my $value (@arrayLength){
	$sum += ( $value - $avLength )**2; 
}
$sdDeviation = sqrt ( $sum / $countEntries );

print "Total number of entries in the file is : $countEntries\n";
print "Total length of the entries in the file is : $totalLength\n";
print "Total number of the entries in + strand is : $plusEntries\n";
print "Total number of the entries in - strand is : $minusEntries\n";
print "Longest entry is : $longestEntry with a length of $longest\n";
print "Shortest entry is : $shortestEntry with a length of $shortest\n";
printf "Average of the length is : %.3f\n", $avLength;
printf "Standard deviation of the length is : %.3f\n", $sdDeviation;
