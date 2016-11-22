#!/usr/bin/perl -w
use strict;

#Author : Tannishtha Som
#script to find orthologs using BLAST
die "Usage : $0 <seq1-File> <seq2-File> n|p  {n=nucleotide, p=protein}\n" if @ARGV < 3;
sub makedb($$$); #function to make the database from the sequence and type  
sub querydb($$$); #function for BLAST query of a sequence against the other database, according to the type

my $file1 = "Seq1_versus_Genome2"; 
my $file2 = "Seq2_versus_Genome1";
my %hashA; #store sequence and best hits in genomeB for sequence1 
my %hashB; #store sequence and best hits in genomeA for sequence2 

makedb($ARGV[0], $ARGV[1], $ARGV[2]);
querydb($ARGV[0], $ARGV[1], $ARGV[2]);
`rm $ARGV[0].*`; #removal of files created by makeblastdb for seq1
`rm $ARGV[1].*`; #removal of files created by makeblastdb for seq2

open FILE1, "$file1" or die "cannot open";
open FILE2, "$file2" or die "cannot open";
#store the sequence and best hits in database2 for genome1 in hashA
while ( my $line1 = <FILE1> ){
	my @element1 = split (/\t/, $line1);
	chomp(@element1);
	$hashA{$element1[0]} = $element1[1];
}
#store the sequence and best hits in database1 for genome2 in hashB
while ( my $line2 = <FILE2> ){
	my @element2 = split (/\t/, $line2);
	chomp(@element2);
    	$hashB{$element2[0]} = $element2[1];
       
}

close FILE1;
close FILE2;
`rm $file1`; #removing temporary file1
`rm $file2`; #removing temporary file2

print "Reciprocal hits\n";
foreach my $key ( keys %hashA ){
	my $value1 = $hashA{$key};
	if ( exists ($hashB{$value1}) ) {
		my $value2 = $hashB{$value1};
		if ( $value2 eq $key ){
			print "$key => $hashA{$key}\n"; 
		}
	}
}

sub makedb($$$){
	my ($seq1, $seq2, $type) = @_;
	if ( $type eq "n"){
		#Creating databases for the two nucleotide sequences using makeblastdb
		`makeblastdb -in $seq1 -dbtype nucl -out $seq1`;
		`makeblastdb -in $seq2 -dbtype nucl -out $seq2`;
	}
	elsif ( $type eq "p"){
		#Creating databases for the two protein sequences using makeblastdb
		`makeblastdb -in $seq1 -dbtype prot -out $seq1`;
		`makeblastdb -in $seq2 -dbtype prot -out $seq2` ;
	}

}

sub querydb($$$){
	my ($seq1, $seq2, $type) = @_;
        if ( $type eq "n"){ 
                `blastn -query $seq1 -db $seq2 -out $file1 -max_target_seqs 1 -outfmt "6 qacc sacc"`; #max_target_seqs stores the specified number of hits
                `blastn -query $seq2 -db $seq1 -out $file2 -max_target_seqs 1 -outfmt "6 qacc sacc"`;
        }
        elsif ( $type eq "p"){
                `blastp -query $seq1 -db $seq2 -out $file1 -max_target_seqs 1 -outfmt "6 qacc sacc"`;
                `blastp -query $seq2 -db $seq1 -out $file2 -max_target_seqs 1 -outfmt "6 qacc sacc"`;
        }

}

