#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

my %lineage;

my %lineage_feat;

my @SNPstore;

my @indelstore;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Lineage"){
	
		next;
	
	}
	else{
	
		if($line[1] eq "SNP"){
		
			push @SNPstore, $line[3].":".$line[4];

			$lineage_feat{$line[5]}{$line[0]}{$line[1]}{$line[3].":".$line[4]} = $_."\n"; 
		
		}
		else{
		
			push @indelstore, $line[3].":".$line[4];
			
			$lineage_feat{$line[5]}{$line[0]}{$line[1]}{$line[3].":".$line[4]} = $_."\n"; 
		
		}
		
		
	
	}

}

close FILE;

my @out;

push @out, join("\t", ("Lineage", "Var type", "Site on Gene", "Site on Genome", "Var", "Gene", "AAVar", "Synonymous or non-syn", "Number of Var", "Number of lineage isolate", "Percentage of Var"));

push @out, "\n";

for my $gene (sort keys %lineage_feat){

	for my $linea (sort keys %{$lineage_feat{$gene}}){
	
		for my $type (sort keys %{$lineage_feat{$gene}{$linea}}){
		
			for my $varinfo (sort keys %{$lineage_feat{$gene}{$linea}{$type}}){
		
				if($type eq "SNP"){
			
					my $match = grep{$_ eq $varinfo} @SNPstore;
				
					if($match == 1){
				
						push @out, $lineage_feat{$gene}{$linea}{$type}{$varinfo};
				
					}
			
				}
				else{
			
					my $match = grep{$_ eq $varinfo} @indelstore;
				
					if($match == 1){
				
						push @out, $lineage_feat{$gene}{$linea}{$type}{$varinfo};
				
					}
			
				}
				
			}
			
		}
	
	}

}

print @out;
