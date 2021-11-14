#!/usr/bin/perl -w

#=========================================================
#
#			Programm Doc
#
#=========================================================

=head1 Description

	Comparative Genomics Analysis of the pairs genome alignment
	
=head1 Usage

	perl genomewideVarCall.pl -refname <name of reference genome > -aligned <aligned fasta>
	
=head1 Parameters

	-refname	[str]	Input the reference genome name
	
	-aligned	[str]	Input the aligned file
	
	-h/-help	[str] print help
	
=head1 Author

	Dongyan Xiong
	
=head1 Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.05.51 15:22 0.0.2
	
=cut

#=========================================================
#
#			Define Paramters
#
#=========================================================

use warnings;

use strict;

use Bio::SeqIO;

use Getopt::Long;

my $refname = "";

my $aligned = "";

my $help = "";

GetOptions(

	'refname=s' => \$refname ,
	
	'aligned=s' => \$aligned,

	'h|help:s' => \$help,

);

die `pod2text $0` if((!$refname) or (!$aligned) or ($help));

#=========================================================
#
#			Main Programm
#
#=========================================================

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

	my $gap = "-" x length($insertion{$i});
	
	push @out_insert, join("\t", ($i, $gap."->".$insertion{$i}));
	
	push @out_insert, "\n";

}

open(OUT, ">", "insertion.tsv");

print OUT @out_insert;

close OUT;

my @out_delet;

push @out_delet, join("\t", ("Site", "Deletion"));

push @out_delet, "\n";

for my $i (sort{$a<=>$b} keys %deletion){

	my $gap = "-" x length($deletion{$i});
	
	my $delen = length($deletion{$i});
	
	push @out_delet, join("\t", ($i."~".($i + $delen - 1), $deletion{$i}."->".$gap));
	
	push @out_delet, "\n";

}

open(OUT, ">", "deletion.tsv");

print OUT @out_delet;

close OUT;

my @out_snp;

push @out_snp, join("\t", ("Site", "SNP"));

push @out_snp, "\n";

for my $i (sort{$a<=>$b} keys %snp){

	push @out_snp, join("\t", ($i, $snp{$i}));
	
	push @out_snp, "\n";

}

open(OUT, ">", "snp.tsv");

print OUT @out_snp;

close OUT;
