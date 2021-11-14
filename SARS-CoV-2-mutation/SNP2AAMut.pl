#!/usr/bin/perl -w
#Author: Dongyan Xiong
use warnings;

use strict;

use Bio::Seq;

use Bio::SeqIO;

use POSIX;

my %gene_hash;

open(FILE, "<", $ARGV[0]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	my @gene = split/\;/,$line[$#line];
	
	for my $i (@gene){
	
		if($i =~ /Name=/){
		
			$i =~ s/Name=//;
			
			@{$gene_hash{$i}} = ($line[3], $line[4]); 
		
		}
	
	}

}

close FILE;

my $fasta = Bio::SeqIO -> new(

	-file => $ARGV[1],
	
	-format => 'fasta',

);

my $seq = $fasta -> next_seq;

my $genome = $seq -> seq;

my @out;

push @out, join("\t", ("Site", "SNP", "number", "Gene", "AASNP", "Synonymous or non-syn"));

push @out, "\n";

open(FILE, "<", $ARGV[2]);

while(<FILE>){

	chomp;
	
	my @line = split/\t/,$_;
	
	if($line[0] eq "Site"){
	
		next;
	
	}
	else{
	
		for my $gene (keys %gene_hash){
		
			if($line[0] >= @{$gene_hash{$gene}}[0] and $line[0] <= @{$gene_hash{$gene}}[1]){
			
				my $ref_seq = substr($genome, @{$gene_hash{$gene}}[0] - 1, @{$gene_hash{$gene}}[1] - @{$gene_hash{$gene}}[0] + 1);
				
				my @mut = split/->/,$line[1];
				
				my @base = split//,$ref_seq;
				
				$base[$line[0] - @{$gene_hash{$gene}}[0]] = $mut[1];
				
				my $mut_seq = join("", @base);
				
				my $ref_obj = Bio::Seq->new( -seq => $ref_seq, -alphabet => 'dna');
				
				my $ref_AA = $ref_obj->translate(-codontable_id => 1,)->seq;
				
				my $mut_obj = Bio::Seq->new( -seq => $mut_seq, -alphabet => 'dna');
				
				my $mut_AA = $mut_obj->translate(-codontable_id => 1,)->seq;
				
				my @ref_AA_arr = split//,$ref_AA;
				
				my @mut_AA_arr = split//,$mut_AA;
				
				my $syn;
				
				my $AAmut;
				
				my $site;
				
				if($ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)] eq $mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)]){
				
					$site = floor(($line[0] - @{$gene_hash{$gene}}[0])/3) + 1;
					
					$syn = "synonymous";
					
					$AAmut = $ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)].$site.$mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)];
				
				}
				else{
				
					$site = floor(($line[0] - @{$gene_hash{$gene}}[0])/3) + 1;
					
					$syn = "non-synonymous";
					
					$AAmut = $ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)].$site.$mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)];
					
				}
				
				push @out, join("\t", (@line, $gene, $AAmut, $syn));
				
				push @out, "\n";
			
			}
		
		}
	
	}

}

close FILE;

print @out;
