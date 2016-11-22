#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som 
#Script to print the count of occurrences of every repName, repFamily and repClass in RepeatMasker

#declaring the hashes to hold the three columns
my %repName;
my %repClass;
my %repFamily;

open FILE, "$ARGV[0]" or die "Could not open $ARGV[0]";
<FILE>; #read the first line and discard it as it contains the names of the columns

#function to count the number of occurrences of each repName, repFamily and repClass
while ( my $line = <FILE>){
	my @elements = split(/\t/, $line);
	chomp(@elements);
	$repName{$elements[10]}++;
	$repClass{$elements[11]}++;
	$repFamily{$elements[12]}++;
}
close FILE;

#printing the values on STDOUT
print "The number of occurrences of repName:\n";
foreach my $key1 (keys %repName){
	printf "%-20s %s\n", $key1, $repName{$key1};
}

print "\n\nThe number of occurrences of repClass:\n";
foreach my $key2 (keys %repClass){
        printf "%-20s %s\n", $key2, $repClass{$key2};
}

print "\n\nThe number of occurrences of repFamily:\n";
foreach my $key3 (keys %repFamily){
        printf "%-20s %s\n", $key3, $repFamily{$key3};
}
 
