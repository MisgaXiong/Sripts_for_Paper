#!/usr/bin/perl -w


=head1 Description

	Comparative Genomics Analysis of the pairs genome alignment
	
Usage

	perl genomewideVarCall.pl -WorkPath <path> -refname <name of reference genome > -aligned <aligned fasta>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-refname	[str]	Input the reference genome name
	
	-aligned	[str]	Input the aligned file
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.05.51 15:22 0.0.2
	
=cut


use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $analysis_path = "";
my $refname = "";
my $aligned = "";
my $help = "";

GetOptions(
	'WorkPath=s' => \$analysis_path,
	'refname=s' => \$refname,
	'aligned=s' => \$aligned,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$refname) or (!$aligned) or ($help));

my $refseq_nm = $refname;
my $fasta = Bio::SeqIO -> new(
	-file => $aligned,
	-format => 'fasta',
);

my %ref_seq;
my %query_seq;
my $total_align_len;

while(my $seq = $fasta -> next_seq){
	my $id = $seq -> id;
	my $seq_seq = $seq -> seq;
	$seq_seq = uc($seq_seq);
	$total_align_len = length($seq_seq);
	my @base_seq = split//,$seq_seq;
	if($id eq $refseq_nm){
		for my $i (0..$#base_seq){
			$ref_seq{$i+1} = $base_seq[$i];
		}
	}
	else{
		for my $i (0..$#base_seq){
			$query_seq{$i+1} = $base_seq[$i];
		}
	}
}

my (%insertion, %deletion, %snp);
my $orig_site = 1;
my $site = 1;

while($site <= $total_align_len){
	if($ref_seq{$site} eq "-" and $query_seq{$site} ne "-"){
		$insertion{$orig_site-0.5} .= $query_seq{$site};
		$site += 1;
	}
	elsif($ref_seq{$site} ne "-" and $query_seq{$site} eq "-"){
		my $step = 0;
		while(($site+$step) <= $total_align_len and $query_seq{$site+$step} eq "-"){
			$deletion{$orig_site} .= $ref_seq{$site+$step};
			$step += 1;
		};
		$site = $site + $step;
		$orig_site = $orig_site + $step;
	}
	elsif($ref_seq{$site} ne "-" and $query_seq{$site} ne "-" and $ref_seq{$site} ne $query_seq{$site}){
		$snp{$orig_site} = $ref_seq{$site}."->".$query_seq{$site};
		$orig_site += 1;
		$site += 1;
	}
	else{
		$orig_site += 1;
		$site += 1;
	}
}

my @out_insert;
push @out_insert, join("\t",("Site", "Insertion"));
push @out_insert, "\n";
for my $i (sort{$a<=>$b} keys %insertion){
	push @out_insert, join("\t", ($i, "ins:".$insertion{$i}));
	push @out_insert, "\n";
}
open(OUT, ">", $analysis_path."insertion.tsv");
print OUT @out_insert;
close OUT;
my @out_delet;
push @out_delet, join("\t", ("Site", "Deletion"));
push @out_delet, "\n";
for my $i (sort{$a<=>$b} keys %deletion){
	my $delen = length($deletion{$i});
	push @out_delet, join("\t", ($i."~".($i + $delen - 1), "del:".$deletion{$i}));
	push @out_delet, "\n";
}
open(OUT, ">", $analysis_path."deletion.tsv");
print OUT @out_delet;
close OUT;
my @out_snp;
push @out_snp, join("\t", ("Site", "SNP"));
push @out_snp, "\n";
for my $i (sort{$a<=>$b} keys %snp){
	push @out_snp, join("\t", ($i, $snp{$i}));
	push @out_snp, "\n";
}
open(OUT, ">", $analysis_path."snp.tsv");
print OUT @out_snp;
close OUT;
