#!/usr/bin/perl -w

use warnings;

use strict;

use Bio::SeqIO;

my $length_bp = $ARGV[0];

my $fasta = Bio::SeqIO -> new(

	-file => $ARGV[1],

	-format => 'fasta',

);

my @outseq;

while(my $seq =  $fasta -> next_seq){

	my $id = $seq -> id;

	my $seq_seq = $seq -> seq;

	if(length($seq_seq) >= $length_bp){

		push @outseq, ">".$id."\n".$seq_seq."\n";

	}

}

print @outseq;

