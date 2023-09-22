#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;

my %spabundance;
my %microbiome;
open(FILE, "<", $ARGV[0]);
my $rank = 0;
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[1] eq "Species"){
		next;
	}
	else{
		$rank += 1;
		my $sp_name = $line[1];
		$sp_name =~ s/ /-/g;
		$sp_name =~ s/\//-/g;
		if(!$sp_name){
			$sp_name = "NA";
		}
		$spabundance{$sp_name} = $rank;
	}
}
close FILE;

open(FILE, "<", $ARGV[1]);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Accession ID"){
		next;
	}
	else{
		my $sp_name = $line[1];
		$sp_name =~ s/ /-/g;
		$sp_name =~ s/\//-/g;
		if(!$sp_name){
			$sp_name = "NA";
		}
		if($spabundance{$sp_name}){
			$microbiome{$line[0]} = $sp_name;
		}
	}
}
close FILE;

my %species_contigs;
my $fasta = Bio::SeqIO -> new(-file => $ARGV[2], -format => 'fasta'); #Final_Assembly.fasta
while(my $seq = $fasta -> next_seq){
	my ($id, $seq_seq) = ($seq -> id, $seq -> seq);
	if(!$microbiome{$id}){
		next;
	}
	else{
		push @{$species_contigs{$microbiome{$id}}}, 1;
		@{$species_contigs{$microbiome{$id}}}[$#{$species_contigs{$microbiome{$id}}}] = ">".$microbiome{$id}."|seq".($#{$species_contigs{$microbiome{$id}}}+1)."|".length($seq_seq)."bp|\n".$seq_seq."\n";
	}
}
for my $sp_name (sort {$b cmp $a} keys %spabundance){
	open(OUT, ">", $ARGV[3]."/"."sp".$spabundance{$sp_name}."_".$sp_name.".fasta");
	print OUT join("", @{$species_contigs{$sp_name}});
	close OUT;
}


