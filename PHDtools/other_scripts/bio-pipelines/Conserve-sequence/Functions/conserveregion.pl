#!/usr/bin/perl -w

=head1 Description
	Find the most conserved sequence of a special pathogen
Usage
	perl conserveregion.pl -WorkPath <path> -reference <reference genome> -isolates <isolate genomes> -algo <max or shannon>
Parameters
	-WorkPath	[str]	Input the analysis path
	-reference	[str]	Input the reference genome
	-isolates	[str]	Input the isolate genomes
	-algo	[str]	Choose one of the method to calculate the conservatism score [max, shannon]
	-h/-help	[str] print help
Auther
	Dongyan Xiong
Edit Time
	2022.07.16 21:11 0.0.3
=cut

use warnings;
use strict;
use POSIX;
use Bio::SeqIO;
use Getopt::Long;

my ($analysis_path, $reference, $isolates, $algorithm, $help) = ("", "", "", "", "");
my $mafft = "/home/dell/miniconda3/envs/wgs/bin/mafft";
my $perl = "/home/dell/miniconda3/envs/Perl/bin/perl";
my $compar_genome = "/var/www/other_scripts/bio-pipelines/Conserve-sequence/Functions/comparativegenome.pl";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'reference=s' => \$reference,
	'isolates=s' => \$isolates,
	'algo=s' => \$algorithm,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$reference) or (!$isolates) or (!$algorithm) or ($help));
die "The name of the algorithm for computiong conservatism must be one of the max or shannon\n" if(uc($algorithm) ne uc("max") and uc($algorithm) ne uc("shannon"));


my $fasta_ref = Bio::SeqIO -> new(
	-file => $reference,
	-format => 'fasta',
);
my $ref_genome = $fasta_ref -> next_seq;
my $ref_id = $ref_genome -> id;
my $ref_seq = $ref_genome -> seq;
$ref_seq = uc($ref_seq);
my $ref_len = length($ref_seq);

my $kmer = 21;
my %ref_kmer;
for my $i (0..($ref_len - $kmer)){
	my $kmer_seq = substr($ref_seq, $i, $kmer);
	$ref_kmer{$kmer_seq} += 1;
}

my $fasta_iso = Bio::SeqIO -> new(
	-file => $isolates,
	-format => 'fasta',
);

my %seq_hash;
for my $i (0..$ref_len){
	$seq_hash{$i}{"A"} = 0;
	$seq_hash{$i}{"T"} = 0;
	$seq_hash{$i}{"C"} = 0;
	$seq_hash{$i}{"G"} = 0;
	$seq_hash{$i}{"-"} = 0;
	$seq_hash{$i}{"I"} = 0;
}

my @var_records;
push @var_records, join("\t", ("SeqID", "VarType", "Position", "Var"));
push @var_records, "\n";
open(OUT, ">", $analysis_path."Total_vars.tsv");
print OUT @var_records;
undef @var_records;
close OUT;

my $total_seq = 0;
open(RCFILE, ">", $analysis_path."Genome_dirct.log");
while(my $iso_genome = $fasta_iso -> next_seq){
	$total_seq += 1;
	my @ref_base_arr = split//,$ref_seq;
	my @ins_arr;
	my $iso_id = $iso_genome -> id;
	my $iso_seq = $iso_genome -> seq;
	my $iso_len = $iso_genome -> length;
	$iso_seq = uc($iso_seq);
	my ($match, $match_rc) = (0, 0);
	for my $i (0..($iso_len - $kmer)){
		my $kmer_seq = substr($iso_seq, $i, $kmer);
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
		$iso_seq = $iso_seq;
	}
	else{
		$iso_seq = reverse($iso_seq);
		$iso_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
		print RCFILE $iso_id.": Genome direction was reverse complementary.\n";
	}
	open(FILE, ">", $analysis_path."scpltem.fasta");
	print FILE ">".$ref_id."\n".$ref_seq."\n".">".$iso_id."\n".$iso_seq."\n";
	close FILE;
	my $command = $mafft." --thread 25 --quiet ".$analysis_path."scpltem.fasta > ".$analysis_path."scpltem_aln.fasta";
	system(qq($command));
	$command = $perl." ".$compar_genome." -WorkPath ".$analysis_path." -refname ".$ref_id." -aligned ".$analysis_path."scpltem_aln.fasta";
	system(qq($command));
	open(OUT, ">>", $analysis_path."Total_vars.tsv");
	open(FILE, "<", $analysis_path."insertion.tsv");
	for my $ins_file (<FILE>){
		chomp($ins_file);
		my @line = split/\t/,$ins_file;
		if($line[0] ne "Site"){
			print OUT $iso_id."\t"."Insertion"."\t".$line[0]."\t".$line[1]."\n";
			$ref_base_arr[($line[0] - 0.5 - 1)] = 0;
			$seq_hash{($line[0] - 0.5)}{"I"} += 1;
			push @ins_arr, ($line[0] - 0.5);
		}
	}
	close FILE;
	open(FILE, "<", $analysis_path."snp.tsv");
	for my $snp_file (<FILE>){
		chomp($snp_file);
		my @line = split/\t/,$snp_file;
		if($line[0] ne "Site"){
			print OUT $iso_id."\t"."SNP"."\t".$line[0]."\t".$line[1]."\n";
			if(grep{$_ == $line[0]} @ins_arr){
				next;
			}
			else{
				$ref_base_arr[($line[0] - 1)] = 0;
				my @mut_base = split/\-\>/,$line[1];
				if($mut_base[1] eq "R"){
					$seq_hash{$line[0]}{"A"} += 0.5;
					$seq_hash{$line[0]}{"G"} += 0.5;
				}
				elsif($mut_base[1] eq "Y"){
					$seq_hash{$line[0]}{"C"} += 0.5;
					$seq_hash{$line[0]}{"T"} += 0.5;
				}
				elsif($mut_base[1] eq "M"){
					$seq_hash{$line[0]}{"A"} += 0.5;
					$seq_hash{$line[0]}{"C"} += 0.5;
				}
				elsif($mut_base[1] eq "K"){
					$seq_hash{$line[0]}{"G"} += 0.5;
					$seq_hash{$line[0]}{"T"} += 0.5;
				}
				elsif($mut_base[1] eq "S"){
					$seq_hash{$line[0]}{"G"} += 0.5;
					$seq_hash{$line[0]}{"C"} += 0.5;
				}
				elsif($mut_base[1] eq "W"){
					$seq_hash{$line[0]}{"A"} += 0.5;
					$seq_hash{$line[0]}{"T"} += 0.5;
				}
				elsif($mut_base[1] eq "B"){
					$seq_hash{$line[0]}{"G"} += 1/3;
					$seq_hash{$line[0]}{"T"} += 1/3;
					$seq_hash{$line[0]}{"C"} += 1/3;
				}
				elsif($mut_base[1] eq "V"){
					$seq_hash{$line[0]}{"G"} += 1/3;
					$seq_hash{$line[0]}{"A"} += 1/3;
					$seq_hash{$line[0]}{"C"} += 1/3;
				}
				elsif($mut_base[1] eq "D"){
					$seq_hash{$line[0]}{"G"} += 1/3;
					$seq_hash{$line[0]}{"A"} += 1/3;
					$seq_hash{$line[0]}{"T"} += 1/3;
				}
				elsif($mut_base[1] eq "H"){
					$seq_hash{$line[0]}{"A"} += 1/3;
					$seq_hash{$line[0]}{"C"} += 1/3;
					$seq_hash{$line[0]}{"T"} += 1/3;
				}
				elsif($mut_base[1] eq "N"){
					$seq_hash{$line[0]}{"A"} += 0.25;
					$seq_hash{$line[0]}{"T"} += 0.25;
					$seq_hash{$line[0]}{"C"} += 0.25;
					$seq_hash{$line[0]}{"G"} += 0.25;
				}
				else{
					$seq_hash{$line[0]}{$mut_base[1]} += 1;
				}
			}
		}
	}
	close FILE;
	open(FILE, "<", $analysis_path."deletion.tsv");
	for my $del_file (<FILE>){
		chomp($del_file);
		my @line = split/\t/,$del_file;
		if($line[0] ne "Site"){
			print OUT $iso_id."\t"."Deletion"."\t".$line[0]."\t".$line[1]."\n";
			my @del_num = split/~/,$line[0];
			for my $site ($del_num[0]..$del_num[1]){
				$ref_base_arr[($site - 1)] = 0;
				$seq_hash{$site}{"-"} += 1;
			}
		}
	}
	close FILE;
	close OUT;
	for my $j (0..$#ref_base_arr){
		if($ref_base_arr[$j] ne 0){
			$seq_hash{($j + 1)}{$ref_base_arr[$j]} += 1;
		}
	}
	unlink $analysis_path."scpltem.fasta";
	unlink $analysis_path."scpltem_aln.fasta";
	unlink $analysis_path."snp.tsv";
	unlink $analysis_path."deletion.tsv";
	unlink $analysis_path."insertion.tsv";
};
close RCFILE;

my @ref_base_arr = split//,$ref_seq;
open(OUTS, ">", $analysis_path."stat.tsv");
print OUTS "Position"."\t"."Ref"."\t"."A"."\t"."T"."\t"."C"."\t"."G"."\t"."-"."\t"."I"."\n";
open(OUTB, ">", $analysis_path."base_score.tsv");
print OUTB "Position"."\t"."Ref"."\t"."Max Base"."\t"."Each Max Freq"."\t"."Conservatism Score"."\n";
for my $pos (sort{$a<=>$b} keys %seq_hash){
	print OUTS $pos."\t".$ref_base_arr[($pos - 1)]."\t".$seq_hash{$pos}{"A"}."\t".$seq_hash{$pos}{"T"}."\t".$seq_hash{$pos}{"C"}."\t".$seq_hash{$pos}{"G"}."\t".$seq_hash{$pos}{"-"}."\t".$seq_hash{$pos}{"I"}."\n";
	my @base_score = max_freq(($seq_hash{$pos}{"A"}, $seq_hash{$pos}{"T"}, $seq_hash{$pos}{"C"}, $seq_hash{$pos}{"G"}, $seq_hash{$pos}{"-"}, $seq_hash{$pos}{"I"}, $total_seq));
	my $def_score = $base_score[0];
	if($base_score[2] =~ m/\-/ or $base_score[2] =~ m/I/){
		$def_score = 0;
	}
	print OUTB $pos."\t".$ref_base_arr[($pos - 1)]."\t".$base_score[2]."\t".$base_score[1]."\t".$def_score."\n";
}
close OUTB;
close OUTS;

sub shannon_index{
	my @arr = @_;
	my $total = 0;
	my $i_logi = 0;
	for my $i (0..($#arr - 1)){
		if($arr[$i] > 0){
			$total += $arr[$i];
			$i_logi += $arr[$i]*log($arr[$i]);
		}
	}
	my $shannon_deversity_idx;
	if($total == 0){
		$shannon_deversity_idx = "NA";
	}
	else{
		$shannon_deversity_idx = log($total) - $i_logi/$total;
	}
	my $conserve_score;
	if($shannon_deversity_idx eq "NA"){
		$conserve_score = -1;
	}
	else{
		$conserve_score = 1 - $shannon_deversity_idx;
	}
	if($conserve_score < 0){
		return 0;
	}
	else{
		return $conserve_score;
	}
}

sub max_freq{
	my @arr = @_;
	my @base = ("A", "T", "C", "G", "-", "I");
	my $max_base;
	my $max_freq = 0;
	my $indel_freq = 0;
	for my $i (0..($#arr - 1)){
		if($base[$i] eq "-" or $base[$i] eq "I"){
			$indel_freq += $arr[$i]/$arr[$#arr];
		}
		if($max_freq < $arr[$i]/$arr[$#arr]){
			$max_freq = $arr[$i]/$arr[$#arr];
			$max_base = $base[$i]; 
		}
		else{
			next;
		}
	}
	my $conserve_score;
	if(uc($algorithm) eq uc("shannon")){
		$conserve_score = shannon_index(@arr);
	}
	else{
		$conserve_score = $max_freq;
	}
	for my $i (0..($#arr - 1)){
		if($base[$i] ne $max_base and $max_freq == $arr[$i]/$arr[$#arr]){
			$max_base = $max_base."/".$base[$i];
		}
	}
	return($conserve_score, $max_freq, $max_base, $indel_freq);
}
