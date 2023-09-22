#!/usr/bin/perl -w

use warnings;

use strict;

use Bio::SeqIO;

my %IDs;

my @outseqs;

my %seqs;

if($ARGV[0] eq "megahit"){

	open(FILE, "<", $ARGV[1]);

	while(<FILE>){

		chomp;
	
		if($_ =~ m/>/){
	
			$_ =~ s/>//;
		
			my @line = split/ /,$_;
		
			$IDs{$line[0]} = $_;
	
		}

	}

	close FILE;

	my $fasta = Bio::SeqIO -> new(

		-file => $ARGV[1],
	
		-format => 'fasta',

	);

	while(my $seq = $fasta -> next_seq){

		my $id = $seq -> id;
	
		my $seq_seq = $seq -> seq;
	
		my $seq_info = $IDs{$id};
	
		my @line = split/ /,$seq_info;
	
		my $length = $line[$#line];
	
		$length =~ s/len=//;
	
		my $cov = $line[$#line-1];
	
		$cov =~ s/multi=//;
	
		my $sequence = ">".$id."_length_".$length."_cov_".$cov."\n".$seq_seq."\n";
	
		push @{$seqs{$length}}, $sequence;

	}


}
elsif($ARGV[0] eq "spades"){

	my $fasta = Bio::SeqIO -> new(

		-file => $ARGV[1],
	
		-format => 'fasta',

	);
	
	while(my $seq = $fasta -> next_seq){

		my $id = $seq -> id;
	
		my $seq_seq = $seq -> seq;
	
		my $length = $seq -> length;
	
		my $sequence = ">".$id."\n".$seq_seq."\n";
	
		push @{$seqs{$length}}, $sequence;

	}

}

for my $length (sort{$b<=>$a} keys %seqs){

	push @outseqs, @{$seqs{$length}};

}

print @outseqs;

