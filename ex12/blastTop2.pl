#!/usr/bin/perl -w
use strict;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Seq;
use Bio::SearchIO;
use Getopt::Long qw(GetOptions);

my $inputFile;
my $sequenceDB;
my $blastMethod;
my $outputFile;
my @params;
my $report_obj;
my $blast_obj;

#Get the options from the command line
GetOptions('i=s' => \$inputFile,
	'd=s' => \$sequenceDB,
	'm=s' => \$blastMethod,
        'o=s' => \$outputFile) or die "Usage $0 –i [input_file] –d [sequence_db.fa] -m [blast_method] –o [output_file] ; blast_method as nucl or prot\n";

#index the database sequence and stipulate the parameters used by the blastall program by populating an array, @params 
if ( $blastMethod eq 'nucl'){
	print "Indexing database sequence...\n";
	`makeblastdb -in $sequenceDB -dbtype nucl -out $sequenceDB`;
	@params = (-program  => 'blastn', -database => $sequenceDB);
}
elsif ($blastMethod eq 'prot'){
	`makeblastdb -in $sequenceDB -dbtype prot -out $sequenceDB`;
	@params = (-program  => 'blastp', -database => $sequenceDB);
}
else{
	print "Can only handle blastn and blastp\n";
}

print "Writing BLAST results to output file...\n";
open FILE, ">$outputFile";
#create a Bio::Seq object to provide the query
my $str = Bio::SeqIO->new(-file  => $inputFile, -format => 'Fasta');
while (my $input = $str->next_seq()){
	$blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
	$report_obj = $blast_obj->blastall($input);
	while( my $result = $report_obj->next_result){
		my $query = $result->query_name();
		print FILE "Query : $query\n";
		for (my $i=1; $i<=2; $i++){  #writing only the top 2 BLAST results 
			my $hit = $result->next_hit();
         	print FILE "Hit:".$hit->name. ",Evalue:".$hit->significance."\n";  		
		}
		print FILE "\n";
	}
}
close FILE;
print "Removing temporary files...\n";
`rm $sequenceDB.*`; 

