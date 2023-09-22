#!/usr/bin/perl -w

=head1 Description

	Calculate the amino acid change based on the snp information
	
Usage

	perl SNP2AAMut.pl -WorkPath <path> -reference <name of reference genome > -gff <gff file> -snp <SNP file from ATCG_SNP.pl> -codon <codon table>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-reference	[str]	Input the reference file
	
	-gff	[str]	Input the gff file
	
	-snp	[str] Input the SNP file from ATCG_SNP.pl
	
	-codon	[int]	Input the translation table
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.04.22 00:53 0.0.1
	
	2021.05.51 15:22 0.0.2
	
	2022.07.21 00:57 0.0.3
	
=cut

use warnings;
use strict;
use Bio::Seq;
use Bio::SeqIO;
use POSIX;
use Getopt::Long;

my ($analysis_path, $gff, $reference, $snpfile, $help) = ("", "", "", "", "");
my $codon_tbl = 1;
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'gff=s' => \$gff,
	'reference=s' => \$reference,
	'snp=s' => \$snpfile,
	'codon:i' => \$codon_tbl,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$gff) or (!$reference) or (!$snpfile) or ($help));

my %gene_hash;

open(FILE, "<", $analysis_path.$gff);
while(<FILE>){
	chomp;
	if($_ =~ m/#/ or $_=~/^\s*$/){
		next;
	}
	else{
		my @line = split/\t/,$_;
		my @gene = split/\;/,$line[$#line];
		if($line[2] eq "gene"){
			my $name = "FALSE";
			for my $i (@gene){
				if($i =~ m/Name=/){
					$name = "TRUE";
					$i =~ s/Name=//;
					@{$gene_hash{$i}} = ($line[3], $line[4], $line[6]);
					last;
				}
			}
			if($name eq "FALSE"){
				for my $i (@gene){
					if($i !~ m/Name=/ and $i =~ m/ID=/){
						$i =~ s/ID=//;
						@{$gene_hash{$i}} = ($line[3], $line[4], $line[6]);
						last;
					}
				}
			}
		}
	}
}
close FILE;

my $fasta = Bio::SeqIO -> new(
	-file => $analysis_path.$reference,
	-format => 'fasta',
);
my $seq = $fasta -> next_seq;
my $genome = $seq -> seq;
my @out;
push @out, join("\t", ("Site", "SNP", "Gene", "AASNP", "Synonymous or non-syn"));
push @out, "\n";

open(FILE, "<", $analysis_path.$snpfile);

my @mutgenome = split//,$genome;
while(<FILE>){
	chomp;
	my @line = split/\t/,$_;
	if($line[0] eq "Site"){
		next;
	}
	else{
		my @mut = split/->/,$line[1];
		$mutgenome[$line[0] - 1] = $mut[1];
	}
}
close FILE;
my $mut_genome = join("", @mutgenome);

my %basecomp = ('A' => 'T', 'C' => 'G', 'T' => 'A', 'G' => 'C');
open(FILE, "<", $analysis_path.$snpfile);
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
				my $mut_seq = substr($mut_genome, @{$gene_hash{$gene}}[0] - 1, @{$gene_hash{$gene}}[1] - @{$gene_hash{$gene}}[0] + 1);
				if(@{$gene_hash{$gene}}[2] eq "-"){
					$ref_seq = reverse($ref_seq);
					$ref_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
					$mut_seq = reverse($mut_seq);
					$mut_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
				}
				my $ref_obj = Bio::Seq->new(-seq => $ref_seq, -alphabet => 'dna');
				my $ref_AA = $ref_obj->translate(-codontable_id => $codon_tbl,)->seq;
				my $mut_obj = Bio::Seq->new(-seq => $mut_seq, -alphabet => 'dna');
				my $mut_AA = $mut_obj->translate(-codontable_id => $codon_tbl,)->seq;
				my @ref_AA_arr = split//,$ref_AA;
				my @mut_AA_arr = split//,$mut_AA;
				my $syn;
				my $AAmut;
				my $site;
				if(@{$gene_hash{$gene}}[2] eq "+"){
					if($ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)] eq $mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)]){
						$site = floor(($line[0] - @{$gene_hash{$gene}}[0])/3) + 1;
						$syn = "synonymous";
						$AAmut = $ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)].$site.$mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)];
						$line[1] = $mut[0].($line[0] - @{$gene_hash{$gene}}[0] + 1).$mut[1];
					}
					else{
						$site = floor(($line[0] - @{$gene_hash{$gene}}[0])/3) + 1;
						$syn = "non-synonymous";
						$AAmut = $ref_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)].$site.$mut_AA_arr[floor(($line[0] - @{$gene_hash{$gene}}[0])/3)];
						$line[1] = $mut[0].($line[0] - @{$gene_hash{$gene}}[0] + 1).$mut[1];
					}
				}
				else{
					my $gen_len = @{$gene_hash{$gene}}[1] - @{$gene_hash{$gene}}[0] + 1;
					my $pos_rc = $gen_len - ($line[0] - @{$gene_hash{$gene}}[0] + 1) + 1;
					if($ref_AA_arr[floor(($pos_rc)/3)] eq $mut_AA_arr[floor(($pos_rc)/3)]){
						$site = floor(($pos_rc)/3) + 1;
						$syn = "synonymous";
						$AAmut = $ref_AA_arr[floor(($pos_rc)/3)].$site.$mut_AA_arr[floor(($pos_rc)/3)];
						$line[1] = $mut[0].$line[0].$mut[1]."(".$basecomp{$mut[0]}.($gen_len - ($line[0] - @{$gene_hash{$gene}}[0] + 1) + 1).$basecomp{$mut[1]}.")";
					}
					else{
						$site = floor(($pos_rc)/3) + 1;
						$syn = "non-synonymous";
						$AAmut = $ref_AA_arr[floor(($pos_rc)/3)].$site.$mut_AA_arr[floor(($pos_rc)/3)];
						$line[1] = $mut[0].$line[0].$mut[1]."(".$basecomp{$mut[0]}.($gen_len - ($line[0] - @{$gene_hash{$gene}}[0] + 1) + 1).$basecomp{$mut[1]}.")";
					}
				}
				push @out, join("\t", (@line, $gene, $AAmut, $syn));
				push @out, "\n";
				last;
			}
		}
	}
}
close FILE;
open(OUT, ">", $analysis_path."SNPAA.tsv");
print OUT @out;
close OUT;
