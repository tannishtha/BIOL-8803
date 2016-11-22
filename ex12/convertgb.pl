#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long qw(GetOptions);

#script to convert genbank file to Fasta or EMBL file using Bio::SeqIO

my $inputFile;
my $format;
my $outputFile;
my $in;
my $out;

#Get the options from the command line
GetOptions('i=s' => \$inputFile,
	'f=s' => \$format,
	'o=s' => \$outputFile) or die "Usage $0 -i <input_file1> -f <fasta|embl> -o <output_file>\n";

#input file in genbank format
$in  = Bio::SeqIO->new(-file => "$inputFile" , '-format' => 'genbank');

#output file either in fasta or EMBL format
if ($format eq "fasta"){
	$out = Bio::SeqIO->new(-file => ">$outputFile" , '-format' => 'Fasta');
}
elsif ($format eq "embl"){
	$out = Bio::SeqIO->new(-file => ">$outputFile" , '-format' => 'EMBL');
}
else{
	print "Unrecognised format\n";
	exit;
}

#write each entry in input file to the output file   
while ( my $seq = $in->next_seq() ) {
	$out->write_seq($seq);
}
