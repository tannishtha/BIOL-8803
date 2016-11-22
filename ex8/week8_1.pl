#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som
#Script to take in negative numbers at the front and positive numbers at the end of the array

#variable declaration
my $num;
my @array;
my $sum=0;

print "Please enter the numbers. Press 0 to stop entering :\n";

while ( ($num = <STDIN>) != 0 ) {
	if ( $num > 0 ) {
		push (@array, $num);  #positive numbers are added to the back of the array
	}
	elsif ( $num < 0 ) {
		unshift (@array, $num) ; #negative numbers added at the front of the array
	}
}

#calculating the sum of the array 
foreach (@array) {
	$sum += $_;
}

print "The numbers you entered are:\n";
foreach (@array) {
	chomp $_;
}
print join(".",@array); #printing the array elements with a dot

print "\nThe sum of the numbers in the array is : $sum \n";


