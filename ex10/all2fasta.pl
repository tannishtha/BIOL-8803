#!/usr/bin/perl -w
use strict;

#This script converts EMBL, FASTQ, GenBank, PIR, MEGA files to the FASTA file format

die "Usage : $0 <InputFile>\n" if @ARGV < 1;
open FILE, "$ARGV[0]" or die "Cannot open $ARGV[0]"; #Input file

sub moleculeType($); #To determine if sequence contains DNA or amino acid
my @fastaDef; #Array to hold the fasta def line details

#Extracting the file name
my @fileName = split (/\./, $ARGV[0]); #$fileName[0] contains the file name

#Checking for the file format
my $line = <FILE>; 

#EMBL file
if ( $line =~ m/ID/ ){  
	my @firstLine = split (/;/, $line);  #Splitting the first line to get the primary accession number and version number
	my @accession = split (/\s+/, $firstLine[0]); #Takes the primary accession number which appears after ID
	push (@fastaDef, $accession[1]);

	while ( $line = <FILE> ){
		if ( $line =~ m/DE/ ) { #Searches for the description
			my @description = split (/\s\s\s/, $line) ; #DE and description are separated by 3 spaces
			push (@fastaDef, $description[1]);
		}
		if ( $line =~ m/SQ/ ) { #SQ denotes beginning of sequences
			last;
		}				

	}
	$line = <FILE>;
	my @sequence = split (/\s+/, $line);
	my $type = moleculeType($sequence[1]); #finding the molecule type
	if ( $type eq "DNA" ){  #opening the desired file
		open OUTPUT, ">>$fileName[0].fna";
	}
	else{
		open OUTPUT, ">>$fileName[0].faa";
	}

	print OUTPUT ">emb|$fastaDef[0]|$fastaDef[1]";
	$line =~ s/\s+|[0-9]+//g;  #stripping off any number and space from the sequence file
	$line = uc($line);  #capitalising
	print OUTPUT "$line\n";
	while ( $line = <FILE> ){
		if ($line =~ m/\/\//){
			last;
		}
		$line =~ s/\s+|[0-9]+//g;
        	$line = uc($line);
        	print OUTPUT "$line\n";
	}
	close OUTPUT;

}

#Genbank File
elsif ( $line =~ m/LOCUS/ ){
	while ( $line = <FILE> ) {
		if ( $line =~ m/VERSION/ ){  #getting the GI number and Accession+version number 
			my @version = split (/\s+|:/, $line);
			$fastaDef[0] = $version[3];
			$fastaDef[1] = $version[1];
		}
		if ($line =~ m/DEFINITION/ ){  #getting the description
			my @def = split (/\s\s/, $line); #separation by two spaces
		 	$fastaDef[2] = $def[1];
		}
		if ( $line =~ m/ORIGIN/ ){  #found the sequence start
			last;
		}
	}
	$line = <FILE>;
	my @sequence = split (/\s+/, $line);
	my $type = moleculeType($sequence[1]); #finding the molecule type
	if ( $type eq "DNA" ){  #opening the desired file
		open OUTPUT, ">>$fileName[0].fna";
	}
	else{
		open OUTPUT, ">>$fileName[0].faa";
	}

	print OUTPUT ">gi|$fastaDef[0]|$fastaDef[1]|$fastaDef[2]";
	$line =~ s/\s+|[0-9]+//g;  #stripping off any number and space from the sequence file
	$line = uc($line);  #capitalising
	print OUTPUT "$line\n";
	while ( $line = <FILE> ){
		if ($line =~ m/\/\//){
			last;
		}
		$line =~ s/\s+|[0-9]+//g;
        	$line = uc($line);
        	print OUTPUT "$line\n";
	}
	close OUTPUT;
	

}

#FASTQ file
elsif ( $line =~ m/^@/ ){
	my @fastq;
	$line = <FILE>;
	$line =~ s/\.//g;
	chomp($line);
	my $type = moleculeType($line); #finding the molecule type	
	if ( $type eq "DNA" ){  #opening the desired file
		open OUTPUT, ">>$fileName[0].fna";
	}
	else{
		open OUTPUT, ">>$fileName[0].faa";
	}

	seek FILE,0,0;
	while ( $line = <FILE> ){
		if ($line =~ m/^@/){
			$line =~ s/^@//;
			$fastq[0] = $line; #first line is id
		}
		$line = <FILE>;
		$line =~ s/\.//g;
		$fastq[1] = $line; #second line is sequence
		$fastq[2] = <FILE>; #third line is +optional
		$fastq[3] = <FILE>; #fourth line is for quality
		
        	print OUTPUT ">$fastq[0]";
		print OUTPUT "$fastq[1]";
	}	
	
	close OUTPUT;
	
}

#MEGA file
elsif ( $line =~ m/#mega/ ){
	<FILE>;
	<FILE>;
	$line = <FILE>;
	$line =~ s/^#//;
	my @defSeq = split (/\s+/, $line);
	chomp(@defSeq);
	my $type = moleculeType($defSeq[1]); #finding the molecule type
	print $type;
        if ( $type eq "DNA" ){  #opening the desired file
                open OUTPUT, ">>$fileName[0].fna";
        }
        else{
                open OUTPUT, ">>$fileName[0].faa";
        }
	$defSeq[1] = uc($defSeq[1]);
	print OUTPUT ">$defSeq[0]\n";  #each line contains the title and sequence
	print OUTPUT "$defSeq[1]\n";
	while ( $line = <FILE> ){
		if ( $line =~ m/^#/){
			$line =~ s/^#//;        		
			@defSeq = split (/\s+/, $line);
        		chomp(@defSeq);
			print OUTPUT ">$defSeq[0]\n";
			if ( $defSeq[1]){  #printing only if the line contains a sequence
				$defSeq[1] = uc($defSeq[1]);
	        		print OUTPUT "$defSeq[1]\n";
			}		
		}
		elsif ($line =~ m/^\s/){
			next;
		}
		else{
			print OUTPUT "\n";
		}
	}
	close OUTPUT;
}

#PIR file
elsif ( $line =~ m/>/ ){
	open OUTPUT, ">>$fileName[0].faa";
	seek FILE,0,0;
	while ( $line = <FILE> ){
		if ($line =~ m/^>/){ 
			$line =~ s/^>//;
			my @id = split (/;/, $line);  #Gets the identification code from the first line
			chomp(@id);
			print OUTPUT ">pir|$id[1]|";
			$line = <FILE>; #2nd line contains description
			print OUTPUT "$line";
			$line = <FILE>;
		}
		if ($line =~ m/^\s/){
			next;
		}
		if ( $line =~ m/\*/){  #The protein sequence ends with *; so remove it
			$line =~ s/\*//;
		}
		print OUTPUT "$line";
	}
}
else{
	print "Unrecognised format\n";
}

close FILE;


sub moleculeType($) {
	my ($var) = @_;
	if ( $var =~ m/[^acgtACGT]/ ){
		return "Protein";
	}
	else{
		return "DNA";
	}

}
