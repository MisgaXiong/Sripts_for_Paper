#!/usr/bin/perl -w
#Author: Dongyan Xiong

use warnings;

use strict;

use Bio::SeqIO;

my %strain2lineage;

open(FILE, "<", $ARGV[0]);

for my $i (<FILE>){

	chomp($i);
	
	my @line = split/\t/,$i;
	
	if($line[0] eq "Virus Strain Name"){
	
		next;
	
	}
	else{
	
		if($line[4] eq "NA" or $line[4] eq "-" or !$line[4]){
		
			$strain2lineage{$line[1]} = "NA";
		
		}
		else{
		
			my @GISAID = grep{$_ =~ m/EPI_ISL_/} @line;
			
			if($GISAID[0] =~ m/,/){
			
				my @sub_GISAID = split/,/,$GISAID[0];
				
				my @match_GISAID = grep {$_ =~ m/EPI_ISL_/} @sub_GISAID;
				
				$strain2lineage{$match_GISAID[0]} = $line[4];
			
			}
			else{
			
				$strain2lineage{$GISAID[0]} = $line[4];
			
			}
		
		}
	
	}

}

close FILE;

my $fasta = Bio::SeqIO -> new(

	-file => $ARGV[1],
	
	-format => 'fasta',

);

my @out;

push @out, join("\t",("GISAID Name", "GISAID ID", "Lineage"));

push @out, "\n";

while(my $seq = $fasta -> next_seq){

	my $id = $seq -> id;
	
	my @line = split/\|/,$id;
	
	my @EPI = grep {$_ =~ m/EPI_ISL_/} @line;
	
	push @out, join("\t", ($id, $EPI[0], $strain2lineage{$EPI[0]}));
	
	push @out, "\n";

}

open(OUT, ">", "GISAID_ID_Lineage.tsv");

print OUT @out;

close OUT;

undef @out;

my %lineage_stat;

for my $id (keys %strain2lineage){

	$lineage_stat{$strain2lineage{$id}} += 1;

}

push @out, join("\t", ("Lineage", "Number"));

push @out, "\n";

for my $id (keys %lineage_stat){

	push @out, join("\t", ($id, $lineage_stat{$id}));
	
	push @out, "\n";

}

open(OUT, ">", "Total_Lineage_stat.tsv");

print OUT @out;

close OUT;

undef @out;

my @unmatch_seq;

my %hash;

open(FILE, "<", "GISAID_ID_Lineage.tsv");

my $unmatch = 0;

my $match = 0;

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "GISAID Name"){
	
		push @out, $_."\n";
		
		push @unmatch_seq, $_."\n";
	
	}
	else{
	
		if(!$line[2]){
		
			push @unmatch_seq, $_."\n";
			
			$unmatch += 1;
		
		}
		else{
		
			push @out, $_."\n";
			
			$hash{$line[2]} += 1;
			
			$match += 1;
		
		}
	
	}

}

close FILE;

open(OUT, ">", "Match_GISAID_ID_Lineage.tsv");

print OUT @out;

close OUT;

open(OUT, ">", "Unmatch_GISAID_ID_Lineage.tsv");

print OUT @unmatch_seq;

close OUT;

undef @out;

push @out, join("\t", ("Lineage", "Number"));

push @out, "\n";

for my $id (keys %hash){

	push @out, join("\t", ($id, $hash{$id}));
	
	push @out, "\n";

}

open(OUT, ">", "MatchSeq_Lineage_stat.tsv");

print OUT @out;

close OUT;

print "Number of ".$match." seqs from query genomes are matched to NGSDC high quality genome metadata file.\n";

print "Number of ".$unmatch." seqs from query genomes are unmatched to NGSDC high quality genome metadata file.\n";
