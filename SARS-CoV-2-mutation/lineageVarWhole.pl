#!/usr/bin/perl -w
#Author: Dongyan Xiong

use warnings;

use strict;

my %lineage;

my %lineage_stat;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "GISAID ID"){
	
		next;
	
	}
	else{
	
		my $info = join("\t", @line[3..$#line]);

		$lineage{$line[1]}{$line[6]}{$line[2]}{$info} += 1;

		$lineage_stat{$line[1]}{$line[0]} += 1;
	
	}

}

close FILE;

my @out;

push @out, join("\t", ("Lineage", "Var type", "Site on Gene", "Site on Genome", "Var", "Gene", "AAVar", "Synonymous or non-syn", "Number of Var", "Number of lineage isolate", "Percentage of Var"));

push @out, "\n";

my @sig_out;

push @sig_out, join("\t", ("Lineage", "Var type", "Site on Gene", "Site on Genome", "Var", "Gene", "AAVar", "Synonymous or non-syn", "Number of Var", "Number of lineage isolate", "Percentage of Var"));

push @sig_out, "\n";

for my $key (sort keys %lineage){

	for my $gene (sort keys %{$lineage{$key}}){
	
		for my $var (sort keys %{$lineage{$key}{$gene}}){
	
			for my $info (sort keys %{$lineage{$key}{$gene}{$var}}){
		
				my $total_isolate = keys %{$lineage_stat{$key}};
		
				push @out, join("\t", ($key, $var, $info, $lineage{$key}{$gene}{$var}{$info}, $total_isolate, $lineage{$key}{$gene}{$var}{$info}/$total_isolate));
			
				push @out, "\n";
			
				if($lineage{$key}{$gene}{$var}{$info}/$total_isolate >= $ARGV[1] and $total_isolate >= $ARGV[2]){
			
					push @sig_out, join("\t", ($key, $var, $info, $lineage{$key}{$gene}{$var}{$info}, $total_isolate, $lineage{$key}{$gene}{$var}{$info}/$total_isolate));
				
					push @sig_out, "\n";
			
				}
		
			}
	
		}
		
	}

}

open(OUT, ">", "LineageVar.tsv");

print OUT @out;

close OUT;

open(OUT, ">", "sig_LineageVar.tsv");

print OUT @sig_out;

close OUT;
