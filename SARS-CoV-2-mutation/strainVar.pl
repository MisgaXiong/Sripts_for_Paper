#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

my @out;

push @out, join("\t", ("GISAID ID", "Lineage", "Var type", "Site on Gene", "Site on Genome", "Var", "Gene", "AAVar", "Synonymous or non-syn"));

push @out, "\n";

my %gene_hash;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	my @gene = split/\;/,$line[$#line];
	
	for my $i (@gene){
	
		if($i =~ /Name=/){
		
			$i =~ s/Name=//;
			
			@{$gene_hash{$i}} = ($line[3], $line[4]); 
		
		}
	
	}

}

close FILE;

open(FILE, "<", "SNPAA.tsv");

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Site"){
	
		next;
	
	}
	else{
	
		for my $gene (keys %gene_hash){
		
			if($line[0] >= @{$gene_hash{$gene}}[0] and $line[0] <= @{$gene_hash{$gene}}[1]){
			
				my $siteongene = $line[0] - @{$gene_hash{$gene}}[0] + 1;

				push @out, join("\t", $ARGV[1], $ARGV[2], "SNP", $siteongene, @line);
				
				push @out, "\n";
			
			}
		
		}
	
	}

}

close FILE;

open(FILE, "<", "deletion.tsv");

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Site"){
	
		next;
	
	}
	else{
	
		my @delsites = split/\~/,$line[0];
	
		for my $gene (keys %gene_hash){
		
			if($delsites[0] >= @{$gene_hash{$gene}}[0] and $delsites[0] <= @{$gene_hash{$gene}}[1]){
			
				my $delsite_start = $delsites[0] - @{$gene_hash{$gene}}[0] + 1;
				
				my $delsite_end = $delsites[1] - @{$gene_hash{$gene}}[0] + 1;
				
				my $siteongene = join("~", ($delsite_start, $delsite_end));
				
				push @out, join("\t", $ARGV[1], $ARGV[2], "deletion", $siteongene, @line, $gene, "-", "non-synonymous");
				
				push @out, "\n";
			
			}
		
		}
	
	}

}

close FILE;

open(FILE, "<", "insertion.tsv");

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Site"){
	
		next;
	
	}
	else{
	
		for my $gene (keys %gene_hash){
		
			if($line[0] >= @{$gene_hash{$gene}}[0] and $line[0] <= @{$gene_hash{$gene}}[1]){
			
				my $siteongene = $line[0] - @{$gene_hash{$gene}}[0] + 1;
				
				push @out, join("\t", $ARGV[1], $ARGV[2], "insertion", $siteongene, @line, $gene, "-", "non-synonymous");
				
				push @out, "\n";
			
			}
		
		}
	
	}

}

close FILE;

open(OUT, ">", "Strain.var");

print OUT @out;

close OUT;
