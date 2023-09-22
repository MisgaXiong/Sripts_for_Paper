#!/usr/bin/perl -w

=head1 Description
	Extract the amplicon sequence of each isolate to primers alignment
Usage
	perl extractseq2primeralign.pl -WorkPath <path> -genomes <genome of target pathogen fasta> -targetm6 <target m6 file>
Parameters
	-WorkPath	[str]	Input the analysis path
	-genomes	[str]	Input the isolate genomes fasta
	-targetm6	[str]	Input the target m6 file
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.18 22:17 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $analysis_path = "";
my ($isolates, $targetm6, $help) = ("", "", "");

GetOptions(
	'WorkPath=s' => \$analysis_path,
	'genomes=s' => \$isolates,
	'targetm6=s' => \$targetm6,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$isolates) or (!$targetm6) or ($help));

my $fasta = Bio::SeqIO -> new(
	-file => $analysis_path.$isolates,
	-format => 'fasta',
);
open(FILE, "<", $analysis_path.$targetm6);
my %blast_res;
my @unblasted;
my @out;
while(<FILE>){
	chomp;
	my @line = split/\t/,$_;
	if($line[8] < $line[9]){
		@{$blast_res{$line[1]}} = ($line[8], $line[9], "NRC");
	}
	else{
		@{$blast_res{$line[1]}} = ($line[9], $line[8], "RC");
	}
}
close FILE;
while(my $seq = $fasta -> next_seq){
	my $id = $seq -> id;
	my $seq_seq = $seq -> seq;
	my $seq_len = $seq -> length;
	if(!$blast_res{$id}){
		push @unblasted, $id;
	}
	else{
		my ($start, $end) = (0, 0);
		if(!$ARGV[2]){
			if(@{$blast_res{$id}}[0] < 30){
				$start = 1;
			}
			else{
				$start = @{$blast_res{$id}}[0] - 30;
			}
		
			if(@{$blast_res{$id}}[1] + 30 > $seq_len){
				$end = $seq_len;
			}
			else{
				$end = @{$blast_res{$id}}[1] + 30;
			}
		}
		else{
			$start = @{$blast_res{$id}}[0];
			$end = @{$blast_res{$id}}[1];
		}
		my $extracted_seq = substr($seq_seq, ($start - 1), ($end - $start + 1));
		if(@{$blast_res{$id}}[2] eq "RC"){
			$extracted_seq = reverse($extracted_seq);
			$extracted_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
		}
		push @out, ">".$id."_2_PrimerAlign\n".$extracted_seq."\n";
	}
}
open(OUT, ">", $analysis_path."Sequences2PrimerAlign.fasta");
print OUT @out;
close OUT;
open(OUT, ">", $analysis_path."lowquality-genomes-uncovered-Amplicon.txt");
print OUT join("\n", @unblasted);
close OUT;
