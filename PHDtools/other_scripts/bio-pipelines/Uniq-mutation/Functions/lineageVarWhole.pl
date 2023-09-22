#!/usr/bin/perl -w

=head1 Description

	Merge the lineage mutation profile to unique mutation analysis
	
Usage

	perl lineageVarWhole.pl -WorkPath <path> -TotalVar <Total_strainVar.tsv> -p <mutation rate> -minNum <minum number of isolates>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-TotalVar	[str]	Input the Total_strainVar.tsv
	
	-p	[int]	Input the mutation rate, default 0.5
	
	-minNum	[str] Input the minum number of isolates
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
=cut

use warnings;
use strict;
use Getopt::Long;

my ($analysis_path, $tVarfile, $help) = ("", "", "");
my ($percent, $minNum) = (0.5, 5);
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'TotalVar=s' => \$tVarfile,
	'p=f' => \$percent,
	'minNum=i' => \$minNum,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$tVarfile) or ($help));

my %lineage;
my %lineage_stat;
open(FILE, "<", $analysis_path.$tVarfile);
while(<FILE>){
	chomp;
	my @line = split/\t/,$_;
	if($line[0] eq "Seq ID"){
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
push @out, join("\t", ("Lineage", "Mutation type", "Site in Gene", "Site in Genome", "Mutation", "Gene", "Amino acid change", "Synonymous or non-syn", "Number of isolates with mutation", "Number of isolates of this lineage", "Mutation rate in this lineage"));
push @out, "\n";
my @sig_out;
push @sig_out, join("\t", ("Lineage", "Mutation type", "Site in Gene", "Site in Genome", "Mutation", "Gene", "Amino acid change", "Synonymous or non-syn", "Number of isolates with mutation", "Number of lineages of this isolate", "Mutation rate in this lineage"));
push @sig_out, "\n";

for my $key (sort keys %lineage){
	for my $gene (sort keys %{$lineage{$key}}){
		for my $var (sort keys %{$lineage{$key}{$gene}}){
			for my $info (sort keys %{$lineage{$key}{$gene}{$var}}){
				my $total_isolate = keys %{$lineage_stat{$key}};
				push @out, join("\t", ($key, $var, $info, $lineage{$key}{$gene}{$var}{$info}, $total_isolate, $lineage{$key}{$gene}{$var}{$info}/$total_isolate));
				push @out, "\n";
				if($lineage{$key}{$gene}{$var}{$info}/$total_isolate >= $percent and $total_isolate >= $minNum){
					push @sig_out, join("\t", ($key, $var, $info, $lineage{$key}{$gene}{$var}{$info}, $total_isolate, $lineage{$key}{$gene}{$var}{$info}/$total_isolate));
					push @sig_out, "\n";
				}
			}
		}
	}
}
open(OUT, ">", $analysis_path."LineageVar.tsv");
print OUT @out;
close OUT;

open(OUT, ">", $analysis_path."sig_LineageVar.tsv");
print OUT @sig_out;
close OUT;
