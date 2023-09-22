#!/usr/bin/perl -w

=head1 Description
	Calculate the conserved sequence for conserveregion.pl
Usage
	perl findconserve.pl -WorkPath <path> -basescore <base_score.tsv> -stats <stat.tsv> -targetLen <target length>
Parameters
	-WorkPath	[str]	Input the analysis path
	-basescore	[str]	Input the base_score.tsv
	-stats	[str]	Input the stat.tsv
	-targetLen	[int]	Input the target length
	-h/-help	[str] print help
Author
	Dongyan Xiong
Edit Time
	2022.07.16 21:23 0.0.2
=cut

use warnings;
use strict;
use Getopt::Long;

my ($analysis_path, $base_score, $stat_file, $tar_len, $help) = ("", "", "", 200, "");

GetOptions(
	'WorkPath=s' => \$analysis_path,
	'basescore=s' => \$base_score,
	'stats=s' => \$stat_file,
	'targetLen:i' => \$tar_len,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$analysis_path) or (!$base_score) or (!$stat_file) or ($help));

my @score_mat;
my %score_rcd;

open(FILE, "<", $base_score);
while(<FILE>){
	chomp;
	my @line = split/\t/,$_;
	if($line[0] ne "Position" and $line[0] != 0){
		push @score_mat, $line[$#line];
		@{$score_rcd{$line[0]}} = @line[1..4];
	}
}
close FILE;
open(FILE, "<", $stat_file);
while(<FILE>){
	chomp;
	my @line = split/\t/,$_;
	if($line[0] ne "Position" and $line[0] != 0){
		push @{$score_rcd{$line[0]}}, freq_count(@line[2..7]);
	}
}

my %seq_score;
my $max_score = 0;
for my $i (0..($#score_mat - $tar_len + 1)){
	if($max_score < sum(@score_mat[$i..($i + $tar_len - 1)])){
		$max_score = sum(@score_mat[$i..($i + $tar_len - 1)]);
	}
	push @{$seq_score{sum(@score_mat[$i..($i + $tar_len - 1)])}}, ($i + 1)."..".($i + $tar_len);
}

my %union_hash;
for my $i (0..$#{$seq_score{$max_score}}){
	my @pos = split/\.\./,@{$seq_score{$max_score}}[$i];
	my @pos_mat = $pos[0]..$pos[1];
	foreach(@pos_mat){
		$union_hash{$_} = 1;
	}
}

my @conserve_mat = sort{$a<=>$b} keys %union_hash;
my @flag_pos;
for my $i (1..$#conserve_mat){
	if($conserve_mat[$i] - $conserve_mat[($i - 1)] > 1){
		push @flag_pos, $i;
	}
}
my @conserve_region;
if($#flag_pos < 0){
	push @conserve_region, $conserve_mat[0]."..".$conserve_mat[$#conserve_mat];
}
elsif($#flag_pos == 0){
	push @conserve_region, $conserve_mat[0]."..".$conserve_mat[$flag_pos[0] - 1];
	push @conserve_region, $conserve_mat[$flag_pos[0]]."..".$conserve_mat[$#conserve_mat];
}
else{
	for my $i (0..$#flag_pos){
		if($i == 0){
			push @conserve_region, $conserve_mat[0]."..".$conserve_mat[$flag_pos[$i] - 1];
		}
		elsif($i == $#flag_pos){
			push @conserve_region, $conserve_mat[$flag_pos[$i - 1]]."..".$conserve_mat[$flag_pos[$i] - 1];
			push @conserve_region, $conserve_mat[$flag_pos[$i]]."..".$conserve_mat[$#conserve_mat];
		}
		else{
			push @conserve_region, $conserve_mat[$flag_pos[$i - 1]]."..".$conserve_mat[$flag_pos[$i] - 1];
		}
	}
}

my (@out, @conseqs);
if(scalar(@conserve_region) == 1){
	push @out, "Target sequence length inputed: ".$tar_len."\n"."Most conserved region score: ".$max_score."\n"."A total of ".scalar(@conserve_region)." was found to be most conserved.\n";
}
else{
	push @out, "Target sequence length inputed: ".$tar_len."\n"."Most conserved region score: ".$max_score."\n"."A total of ".scalar(@conserve_region)." were found to be most conserved.\n";
}
for my $i (0..$#conserve_region){
	push @out, ("=" x 60)."\n";
	my @pos = split/\.\./,$conserve_region[$i];
	push @out, "Region ".($i + 1)." position: ".$pos[0]."~".$pos[1].":\n";
	my $consensus_seq;
	my @site_info;
	push @site_info, join("\t", ("Position", "Ref Base", "Max Freq Base", "Max Freq", "A", "T", "C", "G", "-", "I", "Score"));
	push @site_info, "\n";
	for my $site ($pos[0]..$pos[1]){
		if(@{$score_rcd{$site}}[1] =~ m/\//){
			$consensus_seq .= "(".@{$score_rcd{$site}}[1].")";
		}
		else{
			$consensus_seq .= @{$score_rcd{$site}}[1];
		}
		push @site_info, join("\t", ($site, @{$score_rcd{$site}}[0..2], @{$score_rcd{$site}}[4..$#{$score_rcd{$site}}], $score_mat[$site - 1]));
		push @site_info, "\n";
	}
	push @out, $consensus_seq."\n";
	push @conseqs, ">Conserved_Seq_".($i+1)."\n".$consensus_seq."\n";
	push @out, @site_info;
	push @out, ("=" x 60)."\n";
}
open(OUT, ">", $analysis_path."Conserved_region_info.tsv");
print OUT @out;
close OUT;
open(OUT, ">", $analysis_path."conserved-sequence.fasta");
print OUT @conseqs;
close OUT;

sub sum{
	my @arr = @_;
	my $total = 0;
	for my $i (0..$#arr){
		$total += $arr[$i];
	}
	return($total);
}

sub freq_count{
	my @arr = @_;
	my $total = sum(@arr);
	for my $i (0..$#arr){
		$arr[$i] /= $total;
	}
	return(@arr);
}
