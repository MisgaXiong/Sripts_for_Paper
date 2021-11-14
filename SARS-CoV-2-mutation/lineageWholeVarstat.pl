#!/usr/bin/perl -w
#Author: Dongyan Xiong

use warnings;

use strict;

my @out;

push @out, join("\t", ("Lineage", "Var type", "Site on Gene", "Site on Genome", "Var", "Gene", "AAVar", "Synonymous or non-syn", "Number of Var", "Number of lineage isolate", "Percentage of Var", "Repeat", "Correlation lineage"));

push @out, "\n";

my %lineage;

my %gene_feat;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Lineage"){
	
		next;
	
	}
	else{

		$lineage{$line[5]}{$line[0]}{$line[1]}{$line[3].":".$line[4]} = $_;

		push @{$gene_feat{$line[5]}{$line[1]}{$line[3].":".$line[4]}}, $line[0];
		
	}

}

close FILE;

for my $gene (sort keys %lineage){

	my @gene_out;
	
	for my $linea (sort keys %{$lineage{$gene}}){
	
		for my $var (sort keys %{$lineage{$gene}{$linea}}){
		
			for my $varinfo (sort keys %{$lineage{$gene}{$linea}{$var}}){
			
				my $repeat = scalar(@{$gene_feat{$gene}{$var}{$varinfo}});
				
				my $correlat = join(";", @{$gene_feat{$gene}{$var}{$varinfo}});
				
				push @gene_out, join("\t", ($lineage{$gene}{$linea}{$var}{$varinfo}, $repeat, $correlat));
				
				push @gene_out, "\n";
			
			}
		
		}
	
	}
	
	push @out, @gene_out;

}

print @out;
