#!/usr/bin/perl -w

=head1 Description

	Report the comparative genomics results according to the NGS
	
Usage

	perl ngscomparative.pl -WorkPath <path> -reference <reference genome> -vcf <vcf file>
	
	or
	
	perl ngscomparative.pl -WorkPath <path> -reference <reference genome> -vcf <vcf file> -gff <gff file> -codon <codon table>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-reference	[str]	Input the reference genome

	-vcf	[str]	Input the vcf file
	
	-gff	[str]	Input the gff file
	
	-codon	[str] Input the translation table
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time
	
	2022.07.23 20:47 0.0.1
	
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use POSIX;
use Getopt::Long;

my $analysis_path = "";
my ($reference, $vcf, $gff, $help) = ("", "", "", "");
my $codon_tbl = 1;
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'vcf=s' => \$vcf,
	'gff:s' => \$gff,
	'codon:i' => \$codon_tbl,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$vcf) or ($help));

my @out;
if($gff){
	push @out, join("\t", ("Chr", "Site", "Mutation type", "Mutation", "Mutation in the gene", "Amino acid change", "Synonymous or non-syn"))."\n";
}
else{
	push @out, join("\t", ("Chr", "Site", "Mutation type", "Mutation"))."\n";
}

my $fasta = Bio::SeqIO -> new(-file => $analysis_path.$reference, -format => 'fasta');
my (%genome, %mutsequence, %mutgenome);
while(my $seq = $fasta -> next_seq){
	my ($id, $seq_seq) = ($seq -> id, $seq -> seq);
	$genome{$id} = $seq_seq;
	@{$mutsequence{$id}} = split//,$seq_seq;
}

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

my %variation;
open(FILE, "<", $analysis_path.$vcf);
while(<FILE>){
	chomp;
	if($_ =~ m/#/ or $_=~/^\s*$/){
		next;
	}
	else{
		my @line = split/\t/;
		if($line[5] >= 60){
			if(length($line[3]) == 1 and length($line[4]) == 1){
				if($gff){
					my $res = "FALSE";
					for my $gene (keys %gene_hash){
						if($line[1] >= @{$gene_hash{$gene}}[0] and $line[1] <= @{$gene_hash{$gene}}[1]){
							$res = "TRUE";
							@{$mutsequence{$line[0]}}[$line[1] - 1] = $line[4];
							@{$variation{"SNP"}{$line[0]}{$line[1]}} = ($line[0], $line[1], $line[3], $line[4], $gene);
							last;
						}
					}
					if($res eq "FALSE"){
						@{$mutsequence{$line[0]}}[$line[1] - 1] = $line[4];
						@{$variation{"SNP"}{$line[0]}{$line[1]}} = ($line[0], $line[1], $line[3], $line[4], "Non-coding region");
					}
				}
				else{
					@{$variation{"SNP"}{$line[0]}{$line[1]}} = ($line[0], $line[1], $line[3]."->".$line[4]);
				}
			}
			else{
				if($line[4] =~ m/\,/){
					my @subline = split/\,/,$line[4];
					$line[4] = $subline[0];
				}
				if(length($line[3]) > length($line[4])){
					my $delseq = substr($line[3], (length($line[4])), (length($line[3]) - length($line[4])));
					my @delsite = (($line[1] + length($line[4])), ($line[1] + length($line[3]) - 1));
					if($gff){
						my $res = "FALSE";
						for my $gene (keys %gene_hash){
							if($delsite[0] >= @{$gene_hash{$gene}}[0] and $delsite[0] <= @{$gene_hash{$gene}}[1]){
								$res = "TRUE";
								@{$variation{"deletion"}{$line[0]}{($delsite[0] + $delsite[1])/2}} = ($line[0], join("~", @delsite), $delseq, $gene);
								last;
							}
						}
						if($res eq "FALSE"){
							@{$variation{"deletion"}{$line[0]}{($delsite[0] + $delsite[1])/2}} = ($line[0], join("~", @delsite), $delseq, "Non-coding region");
						}
					}
					else{
						@{$variation{"deletion"}{$line[0]}{($delsite[0] + $delsite[1])/2}} = ($line[0], join("~", @delsite), $delseq);
					}
				}
				elsif(length($line[3]) < length($line[4])){
					my $insseq = substr($line[4], (length($line[3])), (length($line[4]) - length($line[3])));
					my $inssite = $line[1] + length($line[3]) - 0.5;
					if($gff){
						my $res = "FALSE";
						for my $gene (keys %gene_hash){
							if($inssite >= @{$gene_hash{$gene}}[0] and $inssite <= @{$gene_hash{$gene}}[1]){
								$res = "TRUE";
								@{$variation{"insertion"}{$line[0]}{$inssite}} = ($line[0], $inssite, $insseq, $gene);
								last;
							}
						}
						if($res eq "FALSE"){
							@{$variation{"insertion"}{$line[0]}{$inssite}} = ($line[0], $inssite, $insseq, "Non-coding region");
						}
					}
					else{
						@{$variation{"insertion"}{$line[0]}{$inssite}} = ($line[0], $inssite, $insseq);
					}
				}
			}
		}
	}
}
close FILE;

for my $id (keys %mutsequence){
	$mutgenome{$id} = join("", @{$mutsequence{$id}});
}
my %basecomp = ('A' => 'T', 'C' => 'G', 'T' => 'A', 'G' => 'C');
my $gene_len;
for my $type (("SNP", "deletion", "insertion")){
	for my $chr (keys %{$variation{$type}}){
		for my $site (sort{$a <=> $b} keys %{$variation{$type}{$chr}}){
			if($type eq "SNP"){
				if($gff){
					if(@{$variation{$type}{$chr}{$site}}[4] ne "Non-coding region"){
						my $ref_seq = substr($genome{$chr}, (@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] - 1), (@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] + 1));
						my $mut_seq = substr($mutgenome{$chr}, (@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] - 1), (@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] + 1));
						if(@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[2] eq "-"){
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
						my ($syn, $AAmut, $pos, $point_mut);
						if(@{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[2] eq "+"){
							if($ref_AA_arr[floor(($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0])/3)] eq $mut_AA_arr[floor(($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0])/3)]){
								$syn = "synonymous";
							}
							else{
								$syn = "non-synonymous";
							}
							$pos = floor(($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0])/3) + 1;
							$AAmut = $ref_AA_arr[floor(($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0])/3)].$pos.$mut_AA_arr[floor(($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0])/3)];
							$point_mut = @{$variation{$type}{$chr}{$site}}[2].($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] + 1).@{$variation{$type}{$chr}{$site}}[3];
						}
						else{
							$gene_len = @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] + 1;
							my $pos_rc = $gene_len - ($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[4]}}[0] + 1) + 1;
							if($ref_AA_arr[floor(($pos_rc)/3)] eq $mut_AA_arr[floor(($pos_rc)/3)]){
								$syn = "synonymous";
								
							}
							else{
								$syn = "non-synonymous";
							}
							$pos = floor(($pos_rc)/3) + 1;
							$AAmut = $ref_AA_arr[floor(($pos_rc)/3)].$pos.$mut_AA_arr[floor(($pos_rc)/3)];
							$point_mut = @{$variation{$type}{$chr}{$site}}[2].$site.@{$variation{$type}{$chr}{$site}}[3]."(".$basecomp{@{$variation{$type}{$chr}{$site}}[2]}.$pos_rc.$basecomp{@{$variation{$type}{$chr}{$site}}[3]}.")";
						}
						push @out, join("\t", ($chr, $site, $type, $point_mut, @{$variation{$type}{$chr}{$site}}[4], $AAmut, $syn))."\n";
					}
					else{
						push @out, join($chr, $site, $type, @{$variation{$type}{$chr}{$site}}[2]."->".@{$variation{$type}{$chr}{$site}}[3], @{$variation{$type}{$chr}{$site}}[4], "No amino acid change", "NA")."\n";
					}
				}
				else{
					push @out, join("\t", ($chr, $site, $type, @{$variation{$type}{$chr}{$site}}[2]))."\n";
				}
			}
			elsif($type eq "deletion"){
				if($gff){
					my @delsite = split/\~/,@{$variation{$type}{$chr}{$site}}[1];
					my $delseq = @{$variation{$type}{$chr}{$site}}[2];
					if(@{$variation{$type}{$chr}{$site}}[3] ne "Non-coding region"){
						my ($del_start, $del_end, $del_record);
						if(@{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[2] eq "-"){
							$gene_len = @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1;
							$del_start = $gene_len - ($delsite[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1) + 1;
							$del_end = $gene_len - ($delsite[0] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1) + 1;
							my $delseq_rc = reverse($delseq);
							$delseq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
							$del_record = "del:".$delsite[0]."-".$delseq."-".$delsite[1]."(".$del_start."-".$delseq_rc."-".$del_end.")";
						}
						else{
							$del_start = $delsite[0] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1;
							$del_end = $delsite[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1;
							$del_record = "del:".$del_start."-".$delseq."-".$del_end;
						}
						push @out, join("\t", ($chr, @{$variation{$type}{$chr}{$site}}[1], $type, $del_record, @{$variation{$type}{$chr}{$site}}[3], "-", "non-synonymous"))."\n";
					}
					else{
						push @out, join("\t", ($chr, @{$variation{$type}{$chr}{$site}}[1], $type, "del:".$delsite[0]."-".$delseq."-".$delsite[1], @{$variation{$type}{$chr}{$site}}[3], "No amino acid change", "NA"))."\n";
					}
				}
				else{
					push @out, join("\t", ($chr, @{$variation{$type}{$chr}{$site}}[1], $type, "del:".@{$variation{$type}{$chr}{$site}}[2]))."\n";
				}
			}
			elsif($type eq "insertion"){
				if($gff){
					if(@{$variation{$type}{$chr}{$site}}[3] ne "Non-coding region"){
						my $insseq = @{$variation{$type}{$chr}{$site}}[2];
						my ($ins_pos, $ins_record);
						if(@{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[2] eq "-"){
							$gene_len = @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[1] - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1;
							$ins_pos = $gene_len - ($site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1) + 1;
							my $insseq_rc = reverse($insseq);
							$insseq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
							$ins_record = "ins:".$site."-".$insseq."(".$ins_pos."-".$insseq_rc.")";
						}
						else{
							$ins_pos = $site - @{$gene_hash{@{$variation{$type}{$chr}{$site}}[3]}}[0] + 1;
							$ins_record = "ins:".$ins_pos."-".$insseq;
						}
						push @out, join("\t", ($chr, $site, $type, $ins_record, @{$variation{$type}{$chr}{$site}}[3], "-", "non-synonymous"))."\n";
					}
					else{
						push @out, join("\t", ($chr, $site, $type, "ins:".$site."-".@{$variation{$type}{$chr}{$site}}[2], @{$variation{$type}{$chr}{$site}}[3], "No amino acid change", "NA"))."\n";
					}
				}
				else{
					push @out, join("\t", ($chr, $site, $type, "ins:".@{$variation{$type}{$chr}{$site}}[2]))."\n";
				}
			}
		}
	}
}

open(OUT, ">", $analysis_path."comparative-results.tsv");
print OUT @out;
close OUT;
