#!/usr/bin/perl -w

use warnings;
use strict;

my @out;
push @out, join("\t", ("Lineage", "Mutation type", "Site in Gene", "Site in Genome", "Mutation", "Gene", "Amino acid change", "Synonymous or non-syn", "Number of isolates of with mutation", "Number of isolates of this lineage", "Mutation rate in this lineage", "Number of lineages with the same mutation", "Lineages that shared this same mutation"));
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
