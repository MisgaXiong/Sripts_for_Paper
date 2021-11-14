#!/usr/bin/perl -w

use warnings;

use strict;

use Bio::SeqIO;

my %strain2lineage;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "GISAID Name"){
	
		next;
	
	}
	else{
	
		$strain2lineage{$line[1]} = $line[2];
	
	}

}

close FILE;

my $fasta = Bio::SeqIO -> new(

	-file => $ARGV[1],
	
	-format => 'fasta',

);

while(my $seq = $fasta -> next_seq){

	my $GISAID_name = $seq -> id;
	
	my $seq_seq = $seq -> seq;
	
	my @line = split/\|/,$GISAID_name;
	
	my @id = grep {$_ =~ m/EPI_ISL_/} @line;
	
	if($strain2lineage{$id[0]}){
	
		open(OUT, ">", "Seq.fasta");
		
		print OUT ">".$GISAID_name."|".$id[0]."|".$strain2lineage{$id[0]}."\n".$seq_seq."\n";
		
		close OUT;
		
		system(qq(bash linVar.sh));
	
	}

}
