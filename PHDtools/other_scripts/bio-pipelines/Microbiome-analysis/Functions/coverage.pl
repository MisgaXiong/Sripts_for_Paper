#!/usr/bin/perl -w

=head1 Description
	Calculate the genome coverage of mNGS results
Usage
	perl coverage.pl -microbes <microbiome_sort.tsv> -reflens <genomes.csv>
Parameters
	-microbes	[str]	Input the microbiome_sort.tsv
	-reflens	[str]	Input the genomes.csv
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.10.26 19:38 0.0.1
=cut

use warnings;
use strict;
use Getopt::Long;

my ($microbes, $reflens, $help) = ("", "", "");
GetOptions(
	'microbes=s' => \$microbes,
	'reflens=s' => \$reflens,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$microbes) or (!$reflens) or ($help));

my @out;
my %specieslen;
open(FILE, "<", $microbes);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Species"){
		push @out, join("\t", ("Kingdom", @line, "Average reference length (Mbp)", "%Coverage"))."\n";
		open(REF, "<", $reflens);
		while(<REF>){
			chomp;
			$_ =~ s/\"//g;
			my @row = split/\,/;
			my @dot = split/\;/,$row[1];
			@{$specieslen{uc($row[0])}} = ($dot[0], $row[2]);
		}
		close REF;
	}
	else{
		my @subline = split/\s/,$line[0];
		my $reflen = 0;
		my $kingdom = "NA";
		if(!$specieslen{uc($line[0])}){
			for my $i (0..$#subline){
				if(!$specieslen{uc(join(" ", (@subline[0..($#subline - $i)])))}){
					next;
				}
				else{
					$reflen = @{$specieslen{uc(join(" ", (@subline[0..($#subline - $i)])))}}[1];
					$kingdom = @{$specieslen{uc(join(" ", (@subline[0..($#subline - $i)])))}}[0];
					last;
				}
			}
		}
		else{
			$reflen = @{$specieslen{uc($line[0])}}[1];
			$kingdom = @{$specieslen{uc($line[0])}}[0];
		}
		if($reflen == 0){
			push @out, join("\t", ("NA", @line, "NA", "NA"))."\n";
		}
		else{
			if($line[1]*$line[2]/($reflen*1000*1000) > 1){
				push @out, join("\t", ($kingdom, @line, $reflen, "90.0"))."\n";
			}
			else{
				push @out, join("\t", ($kingdom, @line, $reflen, 100*$line[1]*$line[2]/($reflen*1000*1000)))."\n";
			}
		}
	}
}
close FILE;
open(OUT, ">", $microbes."2");
print OUT @out;
close OUT;
