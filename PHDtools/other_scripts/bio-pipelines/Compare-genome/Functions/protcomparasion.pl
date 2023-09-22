#!/usr/bin/perl -w

=head1 Description

	Report the comparative results according to the protein sequences
	
Usage

	perl protcomparasion.pl -WorkPath <path> -snp <snp file> -del <deletion file> -ins <insertion file>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-snp	[str]	Input the snp file

	-del	[str]	Input the deletion file
	
	-ins	[str]	Input the insertion file
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time
	
	2022.07.24 23:23 0.0.1
	
=cut

use strict;
use warnings;
use Getopt::Long;

my $analysis_path = "";
my ($snp, $del, $ins, $help) = ("", "", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'snp=s' => \$snp,
	'del=s' => \$del,
	'ins=s' => \$ins,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$snp) or (!$del) or (!$ins) or ($help));

my @out;
push @out, join("\t", ("Site", "Mutation type", "Mutation"))."\n";
open(FILE, "<", $analysis_path.$snp);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		push @out, join("\t", ($line[0], "point mutation", $line[1]))."\n";
	}
}
close FILE;

open(FILE, "<", $analysis_path.$del);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		push @out, join("\t", ($line[0], "deletion", $line[1]))."\n";
	}
}
close FILE;

open(FILE, "<", $analysis_path.$ins);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		push @out, join("\t", ($line[0], "insertion", $line[1]))."\n";
	}
}
close FILE;

open(OUT, ">", $analysis_path."comparative-results.tsv");
print OUT @out;
close OUT;
