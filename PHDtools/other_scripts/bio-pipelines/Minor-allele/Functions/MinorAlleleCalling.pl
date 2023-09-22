#!/usr/bin/perl -w

=head1 Description
	Minor allele identification
=head1 Usage
	perl MinorAlleleCalling.pl -reference <reference.fasta> -baseProfile <4 base file> -output <output file>
=head1 Parameters
	-reference			[str]	Input the reference genome fasta file, required
	-baseProfile		[str]	Input the 4 base file, required
	-output/out/o	[str]	Input the name of the output file, required
	-gff			[str]	Input the gene annotation file, if the amino acid change want to be identified, optional
	-mindepth		[int]	Minimum number of sequencing depth, default = 100, recommoned >= 100, optional
	-minMuF			[float]	Minimum number of minor allele frequency, default = 0.02, optional
	-codon			[int]	Translation table, if the -gff was chosen, default = 1, optional
	-h/-help		[str]	print help
=head1 Author
	Dongyan Xiong
=head1 Edit Time
	2022.06.23 14:31 0.0.1
=cut

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use POSIX;
use POSIX qw(log10);
use Getopt::Long;

my ($input_ref, $input_6base, $input_gff, $output, $help) = ("", "", "", "", "");
my ($min_dep, $min_minor, $codon_tbl) = (100, 0.02, 1);

GetOptions(
	'reference=s' => \$input_ref,
	'baseProfile=s' => \$input_6base,
	'gff:s' => \$input_gff,
	'mindepth:i' => \$min_dep,
	'minMuF:f' => \$min_minor,
	'codon:i' => \$codon_tbl,
	'output|out|o=s' => \$output,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$input_ref) or (!$input_6base) or (!$output) or ($help));

my $reference = Bio::SeqIO -> new(-file => $input_ref, -format => 'fasta');
my (%ref_map, %ref_base);
while(my $seq = $reference -> next_seq){
	my ($id, $seq_seq) = ($seq -> id, $seq -> seq);
	my @sequences = split//,$seq_seq;
	@{$ref_map{$id}} = @sequences;
	for my $i (0..$#sequences){
		$ref_base{$id}{$i+1}{"A"} = 0;
		$ref_base{$id}{$i+1}{"T"} = 0;
		$ref_base{$id}{$i+1}{"C"} = 0;
		$ref_base{$id}{$i+1}{"G"} = 0;
		$ref_base{$id}{$i+1}{"-"} = 0;
		$ref_base{$id}{$i+1}{"I"} = 0;
		$ref_base{$id}{$i+1}{$sequences[$i]} = 1;
	}
}

my @out;
my %genetbl;
if($input_gff){
	push @out, join("\t", ("Chr", "Position", "Depth", "Reference", "A", "T", "C", "G", "-", "I", "A%", "T%", "C%", "G%", "-%", "I%", "Major base", "Minor base","Minor allele frequency", "RMSD value", "Shannon value", "Type"))."\n";
	open(FILE, "<", $input_gff);
	for my $row (<FILE>){
		chomp($row);
		if($row =~ m/#/){
			next;
		}
		else{
			my @line = split/\t/,$row;
			if($line[2] eq "gene"){
				my @subline = split/\;/,$line[$#line];
				my @genename = grep{$_ =~ m/ID\=/} @subline;
				my $gene_id = $genename[0];
				$gene_id =~ s/ID\=//;
				@{$genetbl{$line[0]}{$gene_id}} = ($line[3], $line[4], $line[6]);
			}
		}
	}
	close FILE;
}
else{
	push @out, join("\t", ("Chr", "Position", "Depth", "Reference", "A", "T", "C", "G", "-", "I", "A%", "T%", "C%", "G%", "-%", "I%", "Major base", "Minor base", "Minor allele frequency", "RMSD value", "Shannon value"))."\n";
}

open(FILE, "<", $input_6base);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if($line[0] eq "Chr"){
		next;
	}
	else{
		my ($chr, $pos, $depth, $ref_base) = ($line[0], $line[1], sum(@line[2..7]), @{$ref_map{$line[0]}}[$line[1] - 1]);
		my ($A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate) = (0, 0, 0, 0, 0, 0);
		my $type;
		if($depth > 0){
			($A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate) = ($line[2]/$depth, $line[3]/$depth, $line[4]/$depth, $line[5]/$depth, $line[6]/$depth, $line[7]/$depth);
		}
		my ($RMSD, $shannon, $muaf, $major_base, $minor_base);
		if($depth >= $min_dep){
			$RMSD = ((($A_rate - $ref_base{$chr}{$pos}{"A"})**2 + ($T_rate - $ref_base{$chr}{$pos}{"T"})**2 + ($C_rate - $ref_base{$chr}{$pos}{"C"})**2 + ($G_rate - $ref_base{$chr}{$pos}{"G"})**2 + ($D_rate - $ref_base{$chr}{$pos}{"-"})**2 + ($I_rate - $ref_base{$chr}{$pos}{"I"})**2) / 6)**(1/2);
			for my $value (($A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate)){
				if($value == 0){
					next;
				}
				else{
					$shannon += $value * (log10($value) / log10(6));
				}
			}
			$shannon = 0 - $shannon;
			($muaf, $major_base, $minor_base) = MuAFcall(($A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate));
			if($muaf >= $min_minor){
				if($input_gff){
					my (@ref, @allele);
					my ($refseq, $alleleseq);
					my $find = "FALSE";
					for my $gene (sort{$a cmp $b} keys %{$genetbl{$chr}}){
						if(@{$genetbl{$chr}{$gene}}[0] <= $pos and $pos <= @{$genetbl{$chr}{$gene}}[1]){
							$find = "TRUE";
							@ref = @{$ref_map{$chr}};
							@allele = @ref;
							$ref[$pos - 1] = $major_base;
							$allele[$pos - 1] = $minor_base;
							$refseq = join("", @ref[(@{$genetbl{$chr}{$gene}}[0] - 1)..(@{$genetbl{$chr}{$gene}}[1] - 1)]);
							$alleleseq = join("", @allele[(@{$genetbl{$chr}{$gene}}[0] - 1)..(@{$genetbl{$chr}{$gene}}[1] - 1)]);
							if(@{$genetbl{$chr}{$gene}}[2] eq "-"){
								$refseq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
								$refseq = reverse($refseq);
								$alleleseq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
								$alleleseq = reverse($alleleseq);
							}
							my $ref_obj = Bio::Seq -> new(-seq => $refseq, -alphabet => 'dna');
							my $allele_obj = Bio::Seq -> new(-seq => $alleleseq, -alphabet => 'dna');
							my $ref_aa = $ref_obj -> translate(-codontable_id => $codon_tbl,) -> seq;
							my $allele_aa = $allele_obj -> translate(-codontable_id => $codon_tbl,) -> seq;
							my @ref_AA_arr = split//,$ref_aa;
							my @allele_AA_arr = split//,$allele_aa;
							if($minor_base eq "-"){
								$type = $chr."-".$gene.":".(floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3) + 1).$ref_AA_arr[floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3)].">del";
								last;
							}
							elsif($minor_base eq "I"){
								$type = $chr."-".$gene.":".(floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3) + 1).$ref_AA_arr[floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3)].">ins";
								last;
							}
							elsif($major_base eq "-"){
								$type = $chr."-".$gene.":".(floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3) + 1)."del>".$allele_AA_arr[floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3)];
							}
							elsif($major_base eq "I"){
								$type = $chr."-".$gene.":".(floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3) + 1)."ins>".$allele_AA_arr[floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3)];
							}
							else{
								my $idx = floor(($pos - @{$genetbl{$chr}{$gene}}[0])/3);
								$type = $chr."-".$gene.":".$ref_AA_arr[$idx].($idx+1).$allele_AA_arr[$idx];
								last;
							}
						}
					}
					if($find eq "FALSE"){
						$type = "non-coding region";
					}
				}
			}
			else{
				$muaf = "<".$min_minor;
				$minor_base = "homozygous";
				$type = "##";
			}
		}
		else{
			($RMSD, $shannon) = ("##", "##");
			($muaf, $major_base, $minor_base) = ("##", "##", "##");
			$type = "##";
			$depth .= " (low depth)";
		}
		if(!$input_gff){
			push @out, join("\t", ($chr, $pos, $depth, $ref_base, @line[2..7], $A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate, $major_base, $minor_base, $muaf, $RMSD, $shannon))."\n";
		}
		else{
			push @out, join("\t", ($chr, $pos, $depth, $ref_base, @line[2..7], $A_rate, $T_rate, $C_rate, $G_rate, $D_rate, $I_rate, $major_base, $minor_base, $muaf, $RMSD, $shannon, $type))."\n";
		}
	}
}
open(OUT, ">", $output);
print OUT @out;
close OUT;

sub sum{
	my $total = 0;
	foreach(@_){
		$total += $_;
	}
	return $total;
}

sub MuAFcall{
	my @base = ("A", "T", "C", "G", "-", "I");
	my @arr = @_;
	my ($muaf, $max, $mid) = (0, 0, 0);
	my ($max_id, $mid_id) = (0, 0);
	for my $i (0..$#arr){
		if($max < $arr[$i]){
			$mid = $max;
			$max = $arr[$i];
			$mid_id = $max_id;
			$max_id = $i;
		}
		elsif($mid < $arr[$i]){
			$mid = $arr[$i];
			$mid_id = $i;
		}
	}
	$muaf = $mid;
	my ($major_base, $minor_base) = ($base[$max_id], $base[$mid_id]);
	return(($muaf, $major_base, $minor_base));
}
