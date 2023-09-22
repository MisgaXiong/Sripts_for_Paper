#!/usr/bin/perl -w

=head1 Description

	Primers/Probe alignment
	
Usage

	perl PrimerAlign.pl -WorkPath <path> -forward <5'->3' forward seq> -reverse <5'->3' reverse seq> -probe <5'->3' probe seq> -target <5'->3' target seq fasta format>
	perl PrimerAlign.pl -WorkPath <path> -pmlist <primer/probe all should 5'->3'> -target <5'->3' target seq fasta format>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-pmlist	[str]	primers-probe list, sequence of primers or probe must be in 5'->3' direction

	-forward	[str]	5'->3' forward seq
	
	-reverse	[str]	5'->3' reverse seq
	
	-probe	[str]	5'->3' probe seq
	
	-target	[str]	5'->3' target seq
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.10.01 16:22 0.0.1
	2021.10.11 19:08 0.0.2
	2021.10.12 18:13 0.0.3
	2021.10.24 16:34 0.0.4
	
=cut

use warnings;
use strict;
use Getopt::Long;

my $analysis_path = "";
my $list = "";
my $forward = "";
my $reverse = "";
my $probe = "";
my $targetfile = "";
my $help = "";
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'pmlist:s' => \$list,
	'forward:s' => \$forward,
	'reverse:s' => \$reverse,
	'probe:s' => \$probe,
	'target:s' => \$targetfile,
	'h|help:s' => \$help,
);

my %direction;
my %RCbase = ("A" => "T", "T" => "A", "C" => "G", "G" => "C", "R" => "Y", "Y" => "R", "M" => "K", "K" => "M", "S" => "S", "W" => "W", "B" => "V", "V" => "B", "D" => "H", "H" => "D", "N" => "N", "I" => "I", "-" => "-");
if($list){
	open(PMFILE, "<", $analysis_path.$list);
	while(<PMFILE>){
		chomp;
		my @line = split/[\t|\s|,]/,$_;
		if(!$line[2]){
			die "Error! The direction of the primer F or R should be given in the file!\n";
		}
		if(uc($line[0]) eq uc("forward")){
			$forward = uc($line[1]);
			@{$direction{"forward"}} = ($forward, $line[2]);
			if(uc($line[2]) eq "R"){
				$forward =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
				$forward = reverse($forward);
			}
		}
		elsif(uc($line[0]) eq uc("reverse")){
			$reverse = uc($line[1]);
			@{$direction{"reverse"}} = ($reverse, $line[2]);
			if(uc($line[2]) eq "R"){
				$reverse =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
				$reverse = reverse($reverse);
			}
		}
		elsif(uc($line[0]) eq uc("probe")){
			$probe = uc($line[1]);
			@{$direction{"probe"}} = ($probe, $line[2]);
			if(uc($line[2]) eq "R"){
				$probe =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
				$probe = reverse($probe);
			}
		}
	}
}
close PMFILE;
die `pod2text $0` if((!$analysis_path) or (((!$forward) or (!$reverse) or (!$targetfile)) and (!$list)) or ($help));

my ($primer_direct, $primer_align, $target_align, $match_link, $end_site) = ("", "", "", "", 0);
my ($F_identi, $R_identi, $P_identi) = ("", "", "");
my (@forward_mut, @reverse_mut, @probe_mut);
my @align_res;
my %Primermut_stat;
open(FILE, "<", $analysis_path.$targetfile);
open(OUT1, ">", $analysis_path."CompletedMatch.tsv");
open(OUT2, ">", $analysis_path."UncompletedMatch.tsv");
open(OUT3, ">", $analysis_path."PrimerStat.tsv");
close OUT1;
close OUT2;
close OUT3;
while(<FILE>){
	chomp;
	if($_ =~ m/>/){
		my $id = $_;
		my $target = <FILE>;
		chomp($target);
		if($forward and $reverse and $probe){
			($forward, $reverse, $probe, $target) = (uc($forward), uc($reverse), uc($probe), uc($target));
			@align_res = PrimerTargetAlign($forward, $target);
			@forward_mut = var_stat($align_res[0], $align_res[1], "forward");
			$F_identi = $align_res[$#align_res];
			if(uc(@{$direction{"forward"}}[1]) eq "R"){
				$primer_direct = " " x ($align_res[3] - 1)."3' "." " x (length($align_res[0]) - 6)." 5'";
			}
			else{
				$primer_direct = " " x ($align_res[3] - 1)."5' "." " x (length($align_res[0]) - 6)." 3'";
			}
			$primer_align = " " x ($align_res[3] - 1).$align_res[0];
			$target_align = substr($target, 0, ($align_res[3] - 1)).$align_res[1];
			$match_link = " " x ($align_res[3] - 1).$align_res[6];
			$end_site = $align_res[5];
			undef @align_res;
			@align_res = PrimerTargetAlign($probe, $target);
			@probe_mut = var_stat($align_res[0], $align_res[1], "probe");
			$P_identi = $align_res[$#align_res];
			if(uc(@{$direction{"probe"}}[1]) eq "R"){
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."3' "." " x (length($align_res[0]) - 6)." 5'";
			}
			else{
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."5' "." " x (length($align_res[0]) - 6)." 3'";
			}
			$primer_align = $primer_align." " x (($align_res[3] - 1) - $end_site).$align_res[0];
			$target_align = $target_align.substr($target, $end_site, ($align_res[3] - $end_site - 1)).$align_res[1];
			$match_link = $match_link." " x (($align_res[3] - 1) - $end_site).$align_res[6];
			$end_site = $align_res[5];
			undef @align_res;
			@align_res = PrimerTargetAlign($reverse, $target);
			@reverse_mut = var_stat($align_res[0], $align_res[1], "reverse");
			$R_identi = $align_res[$#align_res];
			if(uc(@{$direction{"reverse"}}[1]) eq "R"){
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."3' "." " x (length($align_res[0]) - 6)." 5'"." " x (length($target) - $align_res[5]);
			}
			else{
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."5' "." " x (length($align_res[0]) - 6)." 3'"." " x (length($target) - $align_res[5]);
			}
			$primer_align = $primer_align." " x (($align_res[3] - 1) - $end_site).$align_res[0]." " x (length($target) - $align_res[5]);
			$target_align = $target_align.substr($target, $end_site, ($align_res[3] - $end_site - 1)).$align_res[1].substr($target, $align_res[5], (length($target) - $align_res[5]));
			$match_link = $match_link." " x (($align_res[3] - 1) - $end_site).$align_res[6]." " x (length($target) - $align_res[5]);
			open(OUT1, ">>", $analysis_path."CompletedMatch.tsv");
			open(OUT2, ">>", $analysis_path."UncompletedMatch.tsv");
			open(OUT3, ">>", $analysis_path."PrimerStat.tsv");
			if($F_identi == 1 and $R_identi == 1 and $P_identi == 1){
				print OUT1 $id."\n";
				print OUT1 $primer_direct."\n".$primer_align."\n".$match_link."\n".$target_align."\n";
			}
			else{
				print OUT2 $id."\n";
				print OUT2 $primer_direct."\n".$primer_align."\n".$match_link."\n".$target_align."\n";
			}
			print OUT3 $id."\n";
			print OUT3 "Forward primer identity: ".$F_identi."\n";
			print OUT3 "Reverse primer identity: ".$R_identi."\n";
			print OUT3 "Probe primer identity: ".$P_identi."\n";
			print OUT3 @forward_mut;
			print OUT3 @reverse_mut;
			print OUT3 @probe_mut;
			close OUT1;
			close OUT2;
			close OUT3;
		}
		elsif($forward and $reverse and !$probe){
			($forward, $reverse, $target) = (uc($forward), uc($reverse), uc($target));
			@align_res = PrimerTargetAlign($forward, $target);
			@forward_mut = var_stat($align_res[0], $align_res[1], "forward");
			$F_identi = $align_res[$#align_res];
			if(uc(@{$direction{"forward"}}[1]) eq "R"){
				$primer_direct = " " x ($align_res[3] - 1)."3' "." " x (length($align_res[0]) - 6)." 5'";
			}
			else{
				$primer_direct = " " x ($align_res[3] - 1)."5' "." " x (length($align_res[0]) - 6)." 3'";
			}
			$primer_align = " " x ($align_res[3] - 1).$align_res[0];
			$target_align = substr($target, 0, ($align_res[3] - 1)).$align_res[1];
			$match_link = " " x ($align_res[3] - 1).$align_res[6];
			$end_site = $align_res[5];
			undef @align_res;
			@align_res = PrimerTargetAlign($reverse, $target);
			@reverse_mut = var_stat($align_res[0], $align_res[1], "reverse");
			$R_identi = $align_res[$#align_res];
			if(uc(@{$direction{"reverse"}}[1]) eq "R"){
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."3' "." " x (length($align_res[0]) - 6)." 5'"." " x (length($target) - $align_res[5]);
			}
			else{
				$primer_direct = $primer_direct." " x (($align_res[3] - 1) - $end_site)."5' "." " x (length($align_res[0]) - 6)." 3'"." " x (length($target) - $align_res[5]);
			}
			$primer_align = $primer_align." " x (($align_res[3] - 1) - $end_site).$align_res[0]." " x (length($target) - $align_res[5]);
			$target_align = $target_align.substr($target, $end_site, ($align_res[3] - $end_site - 1)).$align_res[1].substr($target, $align_res[5], (length($target) - $align_res[5]));
			$match_link = $match_link." " x (($align_res[3] - 1) - $end_site).$align_res[6]." " x (length($target) - $align_res[5]);
			open(OUT1, ">>", $analysis_path."CompletedMatch.tsv");
			open(OUT2, ">>", $analysis_path."UncompletedMatch.tsv");
			open(OUT3, ">>", $analysis_path."PrimerStat.tsv");
			if($F_identi == 1 and $R_identi == 1){
				print OUT1 $id."\n";
				print OUT1 $primer_direct."\n".$primer_align."\n".$match_link."\n".$target_align."\n";
			}
			else{
				print OUT2 $id."\n";
				print OUT2 $primer_direct."\n".$primer_align."\n".$match_link."\n".$target_align."\n";
			}
			print OUT3 $id."\n";
			print OUT3 "Forward primer identity: ".$F_identi."\n";
			print OUT3 "Reverse primer identity: ".$R_identi."\n";
			print OUT3 @forward_mut;
			print OUT3 @reverse_mut;
			close OUT1;
			close OUT2;
			close OUT3;
		}
	}
}
close FILE;

my (@mergeout, @mergerealout);
my @base_stage = ("A", "T", "C", "G", "-", "I", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N");
if($forward and $reverse and $probe){
	my ($merge, $realmerge) = mergemut($forward, "forward");
	push @mergeout, @{$merge};
	push @mergeout, "\n";
	push @mergerealout, @{$realmerge};
	push @mergerealout, "\n";
	($merge, $realmerge) = mergemut($probe, "probe");
	push @mergeout, @{$merge};
	push @mergeout, "\n";
	push @mergerealout, @{$realmerge};
	push @mergerealout, "\n";
	($merge, $realmerge) = mergemut($reverse, "reverse");
	push @mergeout, @{$merge};
	push @mergerealout, @{$realmerge};
}
elsif($forward and $reverse and !$probe){
	my ($merge, $realmerge) = mergemut($forward, "forward");
	push @mergeout, @{$merge};
	push @mergeout, "\n";
	push @mergerealout, @{$realmerge};
	push @mergerealout, "\n";
	($merge, $realmerge) = mergemut($reverse, "reverse");
	push @mergeout, @{$merge};
	push @mergerealout, @{$realmerge};
}
open(OUT4, ">", $analysis_path."PrimerMutation-Matrix.tsv");
print OUT4 @mergeout;
close OUT4;
open(OUT4, ">", $analysis_path."RealPrimerMutation-Matrix.tsv");
print OUT4 @mergerealout;
close OUT4;

sub max{
	my $max = $_[0];
	for my $i (1..$#_){
		if($max < $_[$i]){
			$max = $_[$i];
		}
	}
	return $max;
}

sub degbase{
	my ($primer_base, $seq_base) = ($_[0], $_[1]);
	my %deg_pr_hash = ("R" => "AG", "Y" => "CT", "M" => "AC", "K" => "GT", "S" => "GC", "W" => "AT", "B" => "GTC", "V" => "GAC", "D" => "GAT", "H" => "ACT", "N" => "AGCT");
	my %deg_sq_hash = ("R" => "AG", "Y" => "CT", "M" => "AC", "K" => "GT", "S" => "GC", "W" => "AT", "B" => "GTC", "V" => "GAC", "D" => "GAT", "H" => "ACT");
	if((!$deg_pr_hash{$primer_base} and !$deg_sq_hash{$seq_base}) or ($deg_pr_hash{$primer_base} and $deg_sq_hash{$seq_base})){
		if($primer_base eq $seq_base){
			return "TRUE";
		}
		else{
			return "FALSE";
		}
	}
	elsif(!$deg_pr_hash{$primer_base} and $deg_sq_hash{$seq_base}){
		my $num = $deg_sq_hash{$seq_base} =~ s/$primer_base/$primer_base/;
		if($num > 0){
			return "TRUE";
		}
		else{
			return "FALSE";
		}
	}
	elsif($deg_pr_hash{$primer_base} and !$deg_sq_hash{$seq_base}){
		my $num = $deg_pr_hash{$primer_base} =~ s/$seq_base/$seq_base/;
		if($num > 0){
			return "TRUE";
		}
		else{
			return "FALSE";
		}
	}
}

sub PrimerTargetAlign{
	my ($match, $mismatch, $gap) = (1, -0.45, -2);
	my ($seq1, $seq2, $type) = ($_[0], $_[1], $_[2]);
	my @seq1 = split//,"-".$seq1;
	my @seq2 = split//,"-".$seq2;
	my (%score, %track);
	my ($max_point, $max_i, $max_j) = (0, 0, 0);
	for my $i (0..$#seq1){
		for my $j (0..$#seq2){
			if($i == 0 and $j == 0){
				$score{$i}{$j} = 0;
				$track{$i}{$j} = 0;
			}
			elsif($i != 0 and $j ==0){
				$score{$i}{$j} = 0;
				$track{$i}{$j} = 2;
			}
			elsif($i == 0 and $j != 0){
				$score{$i}{$j} = 0;
				$track{$i}{$j} = 3;
			}
			else{
				my $current_score;
				if(degbase($seq1[$i], $seq2[$j]) eq "TRUE"){
					$current_score = $match;
				}
				else{
					$current_score = $mismatch;
				}
				$score{$i}{$j} = max((($score{($i - 1)}{($j - 1)} + $current_score), ($score{($i - 1)}{$j} + $gap), ($score{$i}{($j - 1)} + $gap), 0));
				if($max_point <= $score{$i}{$j}){
					($max_point, $max_i, $max_j) = ($score{$i}{$j}, $i, $j);
				}
				if($score{$i}{$j} == $score{($i - 1)}{($j - 1)} + $current_score){
					$track{$i}{$j} = 1;
					next;
				}
				if($score{$i}{$j} == $score{($i - 1)}{$j} + $gap){
					$track{$i}{$j} = 2;
					next;
				}
				if($score{$i}{$j} == $score{$i}{($j - 1)} + $gap){
					$track{$i}{$j} = 3;
					next;
				}
			}
		}
	}
	my ($seq1_i, $seq2_j) = ($max_i, $max_j);
	my ($align1, $align2, $match_link, $match_base_num) = ("", "", "", 0);
	while($score{$seq1_i}{$seq2_j} > 0){
		if($track{$seq1_i}{$seq2_j} == 1){
			$align1 = $seq1[$seq1_i].$align1;
			$align2 = $seq2[$seq2_j].$align2;
			if(degbase($seq1[$seq1_i], $seq2[$seq2_j]) eq "TRUE"){
				$match_link = "|".$match_link;
				$match_base_num += 1;
			}
			else{
				$match_link = " ".$match_link;
			}
			($seq1_i, $seq2_j) = ($seq1_i - 1, $seq2_j - 1);
		}
		elsif($track{$seq1_i}{$seq2_j} == 2){
			$align1 = $seq1[$seq1_i].$align1;
			$align2 = "-".$align2;
			$match_link = " ".$match_link;
			($seq1_i, $seq2_j) = ($seq1_i - 1, $seq2_j);
		}
		elsif($track{$seq1_i}{$seq2_j} == 3){
			$align1 = "-".$align1;
			$align2 = $seq2[$seq2_j].$align2;
			$match_link = " ".$match_link;
			($seq1_i, $seq2_j) = ($seq1_i, $seq2_j - 1);
		}
	}
	my ($min_i, $min_j) = (($seq1_i + 1), ($seq2_j + 1));
	if($min_i <= 2 and $min_i > 1){
		$align1 = substr($seq1, 0, ($min_i - 1)).$align1;
		$match_link = " " x ($min_i - 1).$match_link;
		$align2 = substr($seq2, ($min_j - ($min_i - 1) - 1), ($min_i - 1)).$align2;
		$min_j = $min_j - ($min_i - 1);
		$min_i = 1;
	}
	if($max_i < $#seq1 and $max_i >= ($#seq1 - 1)){
		$align1 .= substr($seq1, $max_i, ($#seq1 - $max_i));
		$match_link .= " " x ($#seq1 - $max_i);
		$align2 .= substr($seq2, $max_j, ($#seq1 - $max_i));
		$max_j = $max_j + ($#seq1 - $max_i);
		$max_i = $#seq1;
	}
	if($min_i > 2){
		$align1 = substr($seq1, 0, ($min_i - 1)).$align1;
		$align2 = "-" x ($min_i - 1).$align2;
		$match_link = " " x ($min_i - 1).$match_link;
	}
	if($max_i < ($#seq1 - 1)){
		$align1 .= substr($seq1, $max_i, ($#seq1 - $max_i));
		$align2 .= "-" x ($#seq1 - $max_i);
		$match_link .= " " x ($#seq1 - $max_i);
	}
	return ($align1, $align2, $min_i, $min_j, $max_i, $max_j, $match_link, $match_base_num/length($match_link));
}

sub var_stat{
	my ($align1, $align2, $total_align_len) = ($_[0], $_[1], length($_[0]));
	my (%ref_seq, %query_seq);
	my @base_seq = split//,$align1;
	for my $i (0..$#base_seq){
		$ref_seq{$i+1} = $base_seq[$i];
	}
	undef @base_seq;
	@base_seq = split//,$align2;
	for my $i (0..$#base_seq){
		$query_seq{$i+1} = $base_seq[$i];
	}
	my @var_rcd;
	push @var_rcd, join("\t", ("Primer/Probe type", "Mutation site", "Mutation"));
	push @var_rcd, "\n";
	my (%insertion, %deletion, %snp, %total_var);
	my ($orig_site, $site) = (1, 1);
	while($site <= $total_align_len){
		if($ref_seq{$site} eq "-" and $query_seq{$site} ne "-"){
			$insertion{$orig_site-0.5} .= $query_seq{$site};
			$site += 1;
		}
		elsif($ref_seq{$site} ne "-" and $query_seq{$site} eq "-"){
			my $step = 0;
			while(($site+$step) <= $total_align_len and $query_seq{$site+$step} eq "-"){
				$deletion{$orig_site} .= $ref_seq{$site+$step};
				$step += 1;
			};
			$site = $site + $step;
			$orig_site = $orig_site + $step;
		}
		elsif($ref_seq{$site} ne "-" and $query_seq{$site} ne "-" and degbase($ref_seq{$site}, $query_seq{$site}) eq "FALSE"){
			$snp{$orig_site} = $ref_seq{$site}.$orig_site.$query_seq{$site};
			$orig_site += 1;
			$site += 1;
		}
		else{
			$orig_site += 1;
			$site += 1;
		}
	}
	for my $site (keys %snp){
		$total_var{$site} = $snp{$site};
		my @mutation = split/\d+/,$snp{$site};
		$Primermut_stat{$_[2]}{$site}{$mutation[1]} += 1;
	}
	for my $site (keys %insertion){
		$total_var{$site} = "ins:".$insertion{$site};
		$Primermut_stat{$_[2]}{$site}{"I"} += 1;
	}
	for my $site (keys %deletion){
		$total_var{$site} = "del:".$deletion{$site};
		for my $i ($site..($site+length($deletion{$site})-1)){
			$Primermut_stat{$_[2]}{$site}{"-"} += 1;
		}
	}
	for my $site (sort{$a<=>$b} keys %total_var){
		if($total_var{$site} =~ m/del/){
			push @var_rcd, join("\t", ($_[2], $site."~".($site+length($total_var{$site})-5), $total_var{$site}));
			push @var_rcd, "\n";
		}
		else{
			push @var_rcd, join("\t", ($_[2], $site, $total_var{$site}));
			push @var_rcd, "\n";
		}
	}
	return @var_rcd;
}

sub mergemut{
	my ($primer_seq, $primer_type) = ($_[0], $_[1]);
	my @primersite = sort{$a<=>$b} keys %{$Primermut_stat{$primer_type}};
	my @sites = (1..length($primer_seq));
	my %count;
	my @sites_uniq = grep{ ++$count{$_} < 2; } (@primersite, @sites);
	my @sites_uniq_sort = sort{$a<=>$b} @sites_uniq;
	my @primer;
	for my $i (@sites_uniq_sort){
		if(int($i) == $i){
			push @primer, substr($primer_seq, ($i-1), 1);
		}
		else{
			push @primer, "I";
		}
	}
	my @mergeinfo;
	push @mergeinfo, join("\t", ($primer_type, @sites_uniq_sort));
	push @mergeinfo, "\n";
	push @mergeinfo, join("\t", ("base", @primer));
	push @mergeinfo, "\n";
	for my $base (@base_stage){
		my @line;
		push @line, $base;
		for my $site (@sites_uniq_sort){
			if(!$Primermut_stat{$primer_type}{$site}){
				push @line, 0;
			}
			elsif(!$Primermut_stat{$primer_type}{$site}{$base}){
				push @line, 0;
			}
			else{
				push @line, $Primermut_stat{$primer_type}{$site}{$base};
			}
		}
		push @mergeinfo, join("\t", @line);
		push @mergeinfo, "\n";
	}
	my @mergerealinfo;
	if(@{$direction{$primer_type}}[1] eq "R"){
		my @real_sites_sort;
		for my $site (@sites_uniq_sort){
			unshift @real_sites_sort, length(@{$direction{$primer_type}}[0]) - $site + 1;
		}
		my $RCprimer = join("", @primer);
		$RCprimer =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
		$RCprimer = reverse($RCprimer);
		my @RC_primer = split//,$RCprimer;
		push @mergerealinfo, join("\t", ($primer_type, @real_sites_sort));
		push @mergerealinfo, "\n";
		push @mergerealinfo, join("\t", ("base", @RC_primer));
		push @mergerealinfo, "\n";
		for my $base (@base_stage){
			my @line;
			push @line, $base;
			for my $site (@real_sites_sort){
				if(!$Primermut_stat{$primer_type}{length(@{$direction{$primer_type}}[0]) - $site + 1}){
					push @line, 0;
				}
				elsif(!$Primermut_stat{$primer_type}{length(@{$direction{$primer_type}}[0]) - $site + 1}{$RCbase{$base}}){
					push @line, 0;
				}
				else{
					push @line, $Primermut_stat{$primer_type}{length(@{$direction{$primer_type}}[0]) - $site + 1}{$RCbase{$base}};
				}
			}
			push @mergerealinfo, join("\t", @line);
			push @mergerealinfo, "\n";
		}
	}
	else{
		@mergerealinfo = @mergeinfo;
	}
	return (\@mergeinfo, \@mergerealinfo);
}
