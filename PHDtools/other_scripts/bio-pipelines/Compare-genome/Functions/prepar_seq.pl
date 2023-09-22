#!/usr/bin/perl -w

=head1 Description
	Prepare the genomes to compariative genome analysis
=head1 Usage
	perl prepar_seq.pl -WorkPath <path> -reference <reference genome> -isolate <isolate genome>
=head1 Parameters
	-WorkPath	[str]	Input the analysis path
	-reference	[str]	Input the reference genome
	-isolate	[str]	Input the isolate genome
	-h/-help	[str]	print help
=head1 Auther
	Dongyan Xiong
=head1 Edit Time
	2021.05.51 15:22 0.0.2
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $analysis_path = "";
my ($reference, $isolate, $help) = ("", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'isolate=s' => \$isolate,
	'h|help:s' => \$help,
);

die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$isolate) or ($help));

my $refseq = Bio::SeqIO -> new(-file => $analysis_path.$reference, -format => 'fasta');
my $ref_seq = $refseq -> next_seq;
my ($ref_seq_seq, $ref_len) = ($ref_seq -> seq, $ref_seq -> length);
$ref_seq_seq = uc($ref_seq_seq);

my $kmer = 21;
my %ref_kmer;
for my $i (0..($ref_len - $kmer)){
	my $kmer_seq = substr($ref_seq_seq, $i, $kmer);
	$ref_kmer{$kmer_seq} += 1;
}

my $isoseq = Bio::SeqIO -> new(-file => $analysis_path.$isolate, -format => 'fasta');
my $iso_seq = $isoseq -> next_seq;
my ($iso_id, $iso_seq_seq, $iso_len) = ($iso_seq -> id, $iso_seq -> seq, $iso_seq -> length);
$iso_seq_seq = uc($iso_seq_seq);

my ($match, $match_rc) = (0, 0);
for my $i (0..($iso_len - $kmer)){
	my $kmer_seq = substr($iso_seq_seq, $i, $kmer);
	my $kmer_seq_rc = reverse($kmer_seq);
	$kmer_seq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
	if($ref_kmer{$kmer_seq}){
		$match += 1;
	}
	if($ref_kmer{$kmer_seq_rc}){
		$match_rc += 1;
	}
}
if($match >= $match_rc){
	$iso_seq_seq = $iso_seq_seq;
}
else{
	$iso_seq_seq = reverse($iso_seq_seq);
	$iso_seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
}

open(OUT, ">", $analysis_path."compare-sequences.fasta");
print OUT ">Ref\n".$ref_seq_seq."\n".">".$iso_id."\n".$iso_seq_seq."\n";
close OUT;
