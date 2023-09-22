#!/usr/bin/perl -w

=head1 Description
	Extract the target taxid sequence
Usage
	perl extrataxseq.pl -WorkPath <analysis path> -taxidtbl <taxid map> -taxid <taxonomy id> -assemble <assembled genome>
Parameters
	-WorkPath	[str]	Input the analysis path
	-taxidtbl	[str]	Input the taxid map
	-taxid	[int]	Input the species level taxonomy id
	-assemble	[str]	Input the assembled genome
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2021.08.25 13:56 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my ($analysis_path, $taxidtbl, $taxid, $assemble, $help) = ("", "", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'taxidtbl=s' => \$taxidtbl,
	'taxid=i' => \$taxid,
	'assemble=s' => \$assemble,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$taxidtbl) or (!$taxid) or (!$assemble) or ($help));

my %seqids;
open(FILE, "<", $analysis_path.$taxidtbl);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Accession ID"){
		next;
	}
	else{
		if($taxid == $line[1]){
			$seqids{$line[0]} = 1;
		}
	}
}
close FILE;

my @out;
my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$assemble, -format => 'fasta');
while(my $seq = $fasta -> next_seq){
	my ($id, $seq_seq) = ($seq -> id, $seq -> seq);
	if($seqids{$id}){
		push @out, ">".$id."\n".$seq_seq."\n";
	}
}
open(OUT, ">", $analysis_path."Final_assembly.fasta");
print OUT @out;
close OUT;
