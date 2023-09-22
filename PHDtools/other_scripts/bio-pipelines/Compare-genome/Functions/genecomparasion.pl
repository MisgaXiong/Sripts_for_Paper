#!/usr/bin/perl -w

=head1 Description

	Merge the genes comparasion results
	
Usage

	perl genecomparasion.pl -WorkPath <path> -snp <snp file> -del <deletion file> -ins <insertion file> -codon <codon table>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-snp	[str]	Input the SNP file

	-del	[str]	Input the deletion file
	
	-ins	[str]	Input the insertion file
	
	-codon	[int]	Input the translation table
		
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time
	
	2022.07.21 23:56 0.0.1
	
=cut

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use POSIX;
use Getopt::Long;

my $analysis_path = "";
my ($reference, $snp, $del, $ins, $help) = ("", "", "", "", "");
my $codon_tbl = 1;
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'snp=s' => \$snp,
	'del=s' => \$del,
	'ins=s' => \$ins,
	'codon:i' => \$codon_tbl,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$snp) or (!$del) or (!$ins) or ($help));

my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$reference, -format => 'fasta');
my $seq = $fasta -> next_seq;
my $sequence = $seq -> seq;
my @mutsequence = split//,$sequence;
open(FILE, "<", $analysis_path.$snp);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		my @mut = split/->/,$line[1];
		$mutsequence[$line[0] - 1] = $mut[1];
	}
}
close FILE;
my $mutseq = join("", @mutsequence);


my $ref_obj = Bio::Seq -> new(-seq => $sequence, -alphabet => 'dna');
my $ref_AA = $ref_obj -> translate(-codontable_id => $codon_tbl,) -> seq;
my @ref_AA_arr = split//,$ref_AA;
my $mut_obj = Bio::Seq -> new(-seq => $mutseq, -alphabet => 'dna');
my $mut_AA = $mut_obj -> translate(-codontable_id => $codon_tbl,) -> seq;
my @mut_AA_arr = split//,$mut_AA;

my @out;
push @out, join("\t", ("Site", "Mutation type", "Mutation", "Amino acid change", "Synonymous or non-syn"))."\n";

open(FILE, "<", $analysis_path.$snp);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Site"){
		next;
	}
	else{
		my ($syn, $AAmut, $site);
		if($ref_AA_arr[floor(($line[0] - 1)/3)] eq $mut_AA_arr[floor(($line[0] - 1)/3)]){
			$site = floor(($line[0] - 1)/3) + 1;
			$syn = "synonymous";
			$AAmut = $ref_AA_arr[floor(($line[0] - 1)/3)].$site.$mut_AA_arr[floor(($line[0] - 1)/3)];
		}
		else{
			$site = floor(($line[0] - 1)/3) + 1;
			$syn = "non-synonymous";
			$AAmut = $ref_AA_arr[floor(($line[0] - 1)/3)].$site.$mut_AA_arr[floor(($line[0] - 1)/3)];
		}
		push @out, join("\t", ($line[0], "SNP", $line[1], $AAmut, $syn))."\n";
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
		push @out, join("\t", ($line[0], "deletion", $line[1], "-", "non-synonymous"))."\n";
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
		push @out, join("\t", ($line[0], "insertion", $line[1], "-", "non-synonymous"))."\n";
	}
}
close FILE;

open(OUT, ">", $analysis_path."comparative-results.tsv");
print OUT @out;
close OUT;
