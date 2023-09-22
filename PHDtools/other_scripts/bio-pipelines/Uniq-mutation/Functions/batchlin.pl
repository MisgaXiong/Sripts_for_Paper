#!/usr/bin/perl -w

=head1 Description
	Batch analysis lineage mutation
=head1 Usage
	perl batchlin.pl -WorkPath <path> -reference <reference genome> -isolate <isolate genome> -gff <gff file> -lineage <lineage information file>
=head1 Parameters
	-WorkPath	[str]	Input the analysis path
	-reference	[str]	Input the reference genome
	-isolate	[str]	Input the isolate genome
	-gff	[str]	Input the gff file
	-lineage	[str]	Input the lineage information file ("Seq ID"	"lineage")
	-codon	[int]	Input the translation table
	-h/-help	[str]	print help
=head1 Auther
	Dongyan Xiong
=head1 Edit Time
	2021.08.03 23:11 0.0.1
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $mafft = "/home/dell/miniconda3/envs/wgs/bin/mafft";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $compare_genome = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/comparativegenome.pl";
my $AGCTsnp = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/ATCG_SNP.pl";
my $SNP2AA = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/SNP2AAMut.pl";
my $merge = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/mergemutation.pl";
my $strainVar = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/Functions/strainVar.pl";

my $analysis_path = "";
my ($reference, $isolate, $gff, $lintbl, $help) = ("", "", "", "", "");
my $codon_tbl = 1;
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'isolate=s' => \$isolate,
	'gff=s' => \$gff,
	'lineage=s' => \$lintbl,
	'codon:i' => \$codon_tbl,
	'h|help:s' => \$help,
);

die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$isolate) or (!$gff) or (!$lintbl) or ($help));

my %lineageInfo;
open(FILE, "<", $analysis_path.$lintbl);
while(<FILE>){
	chomp;
	my @line = split/\t/;
	if(uc($line[0]) eq "SEQ ID"){
		next;
	}
	else{
		$lineageInfo{$line[0]} = $line[1];
	}
}
close FILE;

my $refseq = Bio::SeqIO -> new(-file => $analysis_path.$reference, -format => 'fasta');
my $ref_seq = $refseq -> next_seq;
my ($ref_seq_seq, $ref_len) = ($ref_seq -> seq, $ref_seq -> length);
$ref_seq_seq = uc($ref_seq_seq);

my $kmer = 21;
my %ref_kmer;
for my $i (0..($ref_len - $kmer)){
	my $kmer_seq = substr($ref_seq_seq, $i, $kmer);
	$ref_kmer{$kmer_seq} += 1;
}

my $isoseq = Bio::SeqIO -> new(-file => $analysis_path.$isolate, -format => 'fasta');
while(my $iso_seq = $isoseq -> next_seq){
	my ($iso_id, $iso_seq_seq, $iso_len) = ($iso_seq -> id, $iso_seq -> seq, $iso_seq -> length);
	$iso_seq_seq = uc($iso_seq_seq);
	if($lineageInfo{$iso_id}){
		my ($match, $match_rc) = (0, 0);
		for my $i (0..($iso_len - $kmer)){
			my $kmer_seq = substr($iso_seq_seq, $i, $kmer);
			my $kmer_seq_rc = reverse($kmer_seq);
			$kmer_seq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
			if($ref_kmer{$kmer_seq}){
				$match += 1;
			}
			if($ref_kmer{$kmer_seq_rc}){
				$match_rc += 1;
			}
		}
		if($match >= $match_rc){
			$iso_seq_seq = $iso_seq_seq;
		}
		else{
			$iso_seq_seq = reverse($iso_seq_seq);
			$iso_seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
		}
		open(OUT, ">", $analysis_path."compare-sequences.fasta");
		print OUT ">Ref\n".$ref_seq_seq."\n".">".$iso_id."\n".$iso_seq_seq."\n";
		close OUT;
		my $command = $mafft." --thread 20 --quiet ".$analysis_path."compare-sequences.fasta > ".$analysis_path."compare-sequences_maf.fasta";
		system(qq($command));
		$command = $perl." ".$compare_genome." -WorkPath ".$analysis_path." -refname Ref -aligned ".$analysis_path."compare-sequences_maf.fasta";
		system(qq($command));
		$command = $perl." ".$AGCTsnp." ".$analysis_path."snp.tsv > ".$analysis_path."ATCG_snp.tsv";
		system(qq($command));
		$command = $perl." ".$SNP2AA." -WorkPath ".$analysis_path." -reference ".$reference." -gff ".$gff." -snp ATCG_snp.tsv -codon ".$codon_tbl;
		system(qq($command));
		$command = $perl." ".$merge." -WorkPath ".$analysis_path." -snp ATCG_snp.tsv -del deletion.tsv -ins insertion.tsv -AAchange SNPAA.tsv -gff ".$gff;
		system(qq($command));
		$command = $perl." ".$strainVar." -WorkPath ".$analysis_path." -seqid ".$iso_id." -lineage ".$lineageInfo{$iso_id}." -varProfile comparative-results.tsv";
		system(qq($command));
		$command = "cat ".$analysis_path."Strain.var >> ".$analysis_path."Total_strainVar.tsv";
		system(qq($command));
	}
	else{
		next;
	}
}
