#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

use Bio::SeqIO;

my @fastas = glob("*fasta");

for my $i (0..$#fastas){

	my @out;

	my $fasta = Bio::SeqIO -> new(

		-file => $fastas[$i],

		-format => 'fasta',

	);

	while(my $seq = $fasta -> next_seq){

		my $id = $seq -> id;

		my $seq_seq = $seq -> seq;

		push @out, ">".$id."\n".$seq_seq."\n";

	};

	open(OUT, ">>", "mergeseq.fa");

	print OUT @out;

	close OUT;

}
