#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my @allele_sets = glob($ARGV[0]."/*sets_maf-tidy.fasta");
my @iso_sets = glob($ARGV[1]."/*fna");
open(OUT, ">", $ARGV[2]."/Total_link.fasta");
for my $isolate (@iso_sets){
	my @line = split/\//,$isolate;
	my $file_name = $line[$#line];
	my $seq_name = $file_name;
	$seq_name =~ s/\.fna//;
	my $iso_link = "";
	for my $allele (@allele_sets){
		open(FASTA, "<", $allele);
		while(my $id = <FASTA>){
			$id =~ s/\n//;
			my @ids = split/\|-_-\|/,$id;
			my $seq_seq = <FASTA>;
			$seq_seq =~ s/\n//;
			if($file_name eq $ids[$#ids]){
				$iso_link = $iso_link.$seq_seq;
				last;
			}
		}
		close FASTA;
	}
	print OUT ">".$seq_name."\n".$iso_link."\n";
}
close OUT;
