#!/usr/bin/perl -w

=head1 Description

	Primer design
	
Usage

	perl PrimerDesign.pl -WorkPath <path> -sequence <target sequence> -type <nPCR or qPCR>
	
	or

	perl PrimerDesign.pl -WorkPath <path> -sequence <target sequence> -type <nPCR or qPCR> -consFILE <conserved score file>
	
Parameters

	-WorkPath	[str]	Input the analysis path

	-sequence	[str]	Input the target nucleotide sequence but no more than 300 bp.
	
	-type	[str]	Input the PCR type (nPCR or qPCR)
	
	-consFILE	[str] Input the converse score file, output from conserved sequence determine
	
	-h/-help	[str] print help
	
Auther

	Dongyan Xiong
	
Edit Time

	2021.08.22 00:53 0.0.1
	
	2022.07.21 00:57 0.0.2
	
=cut

use warnings;
use strict;
use Getopt::Long;

my ($analysis_path, $PCR_type, $tar_seq, $cons_file, $help) = ("", "", "", "", "");
GetOptions(
	'WorkPath=s' => \$analysis_path,
	'sequence=s' => \$tar_seq,
	'type=s' => \$PCR_type,
	'consFILE:s' => \$cons_file,
	'h|help:s' => \$help,
);
$tar_seq = uc($tar_seq);
die `pod2text $0` if((!$analysis_path) or (!$tar_seq) or (!$PCR_type) or ($help));
die `pod2text $0` if(length($tar_seq) > 360);
my $checkseq = $tar_seq;
$checkseq =~ s/A//g;
$checkseq =~ s/T//g;
$checkseq =~ s/C//g;
$checkseq =~ s/G//g;
die `pod2text $0` if(length($checkseq) >= 1);
die `pod2text $0` if($PCR_type ne "qPCR" and $PCR_type ne "nPCR");

my ($minPCR_len, $maxPCR_len) = (70, 300);
my @cons_score;
if($cons_file){
	open(CONSFILE, "<", $analysis_path.$cons_file);
	while(<CONSFILE>){
		next if($. <= 7);
		chomp;
		if($_ =~ m/====/){
			last;
		}
		else{
			my @line = split/\t/,$_;
			push @cons_score, $line[$#line];
		}
	}
}
my $seq_len = length($tar_seq);
if($seq_len <= $maxPCR_len){
	$maxPCR_len = $seq_len;
}
my (%primers, %primers_cons_sort);
if($PCR_type eq "nPCR"){
	for my $F_len (19..23){
		for my $R_len (19..23){
			for my $produce ($minPCR_len..$maxPCR_len){
				for my $start (0..($seq_len - $produce + 1)){
					my $F_primer = substr($tar_seq, $start, $F_len);
					my $R_primer = reverse(substr($tar_seq, ($start + ($produce - 1) - $R_len + 1), $R_len));
					$R_primer =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
					my ($F_tm, $R_tm) = (pm_Tm_calcu($F_primer), pm_Tm_calcu($R_primer));
					my @dimerres;
					push @dimerres, dimer_calcu($F_primer, $F_primer, "N");
					push @dimerres, dimer_calcu($R_primer, $R_primer, "N");
					push @dimerres, dimer_calcu($F_primer, $R_primer, "Y");
					if(grep{$_ eq "nonpass"} @dimerres){
						$F_tm = 0;
						$R_tm = 0;
					}
					if((abs($F_tm - $R_tm) <= 2) and ($F_tm >= 50 and $F_tm <= 75) and ($R_tm >= 50 and $R_tm <= 75)){
						if($cons_file){
							my $avg_score = sum(@cons_score[($start..($start+$F_len-1), ($start + ($produce - 1) - $R_len + 1)..($start+$produce-1))])/($F_len+$R_len);
							$primers_cons_sort{$avg_score}{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nConserve socre:\n".$avg_score."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n";
						}
						else{
							$primers{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n";
						}
					}
				}
			}
		}
	}
	open(OUT, ">", $analysis_path."Primers.tsv");
	if($cons_file){
		for my $avg_score (sort{$b<=>$a} keys %primers_cons_sort){
			for my $primer_score (sort{$a<=>$b} keys %{$primers_cons_sort{$avg_score}}){
				print OUT $primers_cons_sort{$avg_score}{$primer_score};
			}
		}
	}
	for my $primer_score (sort{$a<=>$b} keys %primers){
		print OUT $primers{$primer_score};
	}
	close OUT;
}
elsif($PCR_type eq "qPCR"){
	my %primer_info;
	my $num = 0;
	for my $F_len (18..22){
		for my $R_len (18..22){
			for my $produce ($minPCR_len..$maxPCR_len){
				for my $start (0..($seq_len - $produce + 1)){
					my $F_primer = substr($tar_seq, $start, $F_len);
					my $R_primer = reverse(substr($tar_seq, ($start + ($produce - 1) - $R_len + 1), $R_len));
					$R_primer =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
					my ($F_tm, $R_tm) = (pm_Tm_calcu($F_primer), pm_Tm_calcu($R_primer));
					my @dimerres;
					push @dimerres, dimer_calcu($F_primer, $F_primer, "N");
					push @dimerres, dimer_calcu($R_primer, $R_primer, "N");
					push @dimerres, dimer_calcu($F_primer, $R_primer, "Y");
					if(grep{$_ eq "nonpass"} @dimerres){
						$F_tm = 0;
						$R_tm = 0;
					}
					if((abs($F_tm - $R_tm) <= 2) and ($F_tm >= 48 and $F_tm <= 65) and ($R_tm >= 48 and $R_tm <= 65)){
						$num += 1;
						@{$primer_info{$num}} = ($start, ($start + $F_len - 1), ($start + ($produce - 1)), ($start + ($produce - 1) - $R_len + 1), $F_primer, $F_tm, $R_primer, $R_tm);
					}
				}
			}
		}
	}
	for my $num (sort{$a<=>$b} keys %primer_info){
		my ($F_start, $F_end, $R_start, $R_end, $F_primer, $F_tm, $R_primer, $R_tm) = @{$primer_info{$num}};
		if(($R_start - $F_start + 1 - ($F_end - $F_start + 1) - ($R_start - $R_end + 1)) < (12*2 + 23)){
			for my $pb_len (20..23){
				for my $pb_F_start (($F_end + 1)..($R_end - $pb_len)){
					my $probe;
					my $pb_Tm;
					if(($pb_F_start - $F_end) < ($R_end - ($pb_F_start + $pb_len -1)) and ($pb_F_start - $F_end) < 13){
						$probe = substr($tar_seq, $pb_F_start, $pb_len);
						$pb_Tm = pb_Tm_calcu($probe);
						my @dimerres;
						push @dimerres, dimer_calcu($probe, $probe, "N");
						push @dimerres, dimer_calcu($F_primer, $probe, "Y");
						push @dimerres, dimer_calcu($R_primer, $probe, "Y");
						if(grep{$_ eq "nonpass"} @dimerres){
							$pb_Tm = 0;
						}
						if(($pb_Tm - $F_tm) <= 12 and ($pb_Tm - $F_tm) >= 8 and ($pb_Tm - $R_tm) <= 12 and ($pb_Tm - $R_tm) >= 8){
							if($cons_file){
								my $avg_score = sum(@cons_score[($F_start..$F_end, $R_end..$R_start, $pb_F_start..($pb_F_start+$pb_len-1))])/(($F_end - $F_start+1)+($R_start - $R_end+1)+$pb_len);
								$primers_cons_sort{$avg_score}{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nConserve socre:\n".$avg_score."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
							else{
								$primers{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
						}
					}
					elsif(($R_end - ($pb_F_start + $pb_len - 1)) < ($pb_F_start - $F_end) and ($R_end - ($pb_F_start + $pb_len - 1)) < 13){
						$probe = substr($tar_seq, $pb_F_start, $pb_len);
						$probe = reverse($probe);
						$probe =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
						$pb_Tm = pb_Tm_calcu($probe);
						my @dimerres;
						push @dimerres, dimer_calcu($probe, $probe, "N");
						push @dimerres, dimer_calcu($F_primer, $probe, "Y");
						push @dimerres, dimer_calcu($R_primer, $probe, "Y");
						if(grep{$_ eq "nonpass"} @dimerres){
							$pb_Tm = 0;
						}
						if(($pb_Tm - $F_tm) <= 12 and ($pb_Tm - $F_tm) >= 8 and ($pb_Tm - $R_tm) <= 12 and ($pb_Tm - $R_tm) >= 8){
							if($cons_file){
								my $avg_score = sum(@cons_score[($F_start..$F_end, $R_end..$R_start, $pb_F_start..($pb_F_start+$pb_len-1))])/(($F_end - $F_start+1)+($R_start - $R_end+1)+$pb_len);
								$primers_cons_sort{$avg_score}{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nConserve socre:\n".$avg_score."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
							else{
								$primers{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
						}
					}
				}
			}
		}
		else{
			for my $pb_len (20..23){
				for my $pb_start ((($F_end + 1)..($F_end + 13), ($R_end - 13)..($R_end - 1))){
					my $probe;
					my $pb_Tm;
					if($pb_start - $F_end < 13){
						$probe = substr($tar_seq, $pb_start, $pb_len);
						$pb_Tm = pb_Tm_calcu($probe);
						my @dimerres;
						push @dimerres, dimer_calcu($probe, $probe, "N");
						push @dimerres, dimer_calcu($F_primer, $probe, "Y");
						push @dimerres, dimer_calcu($R_primer, $probe, "Y");
						if(grep{$_ eq "nonpass"} @dimerres){
							$pb_Tm = 0;
						}
						if(($pb_Tm - $F_tm) <= 12 and ($pb_Tm - $F_tm) >= 8 and ($pb_Tm - $R_tm) <= 12 and ($pb_Tm - $R_tm) >= 8){
							if($cons_file){
								my $avg_score = sum(@cons_score[($F_start..$F_end, $R_end..$R_start, $pb_start..($pb_start+$pb_len-1))])/(($F_end - $F_start+1)+($R_start - $R_end+1)+$pb_len);
								$primers_cons_sort{$avg_score}{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nConserve socre:\n".$avg_score."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
							else{
								$primers{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
						}
					}
					elsif($R_end - $pb_start < 13){
						$probe = substr($tar_seq, ($pb_start - $pb_len + 1), $pb_len);
						$probe = reverse($probe);
						$probe =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
						$pb_Tm = pb_Tm_calcu($probe);
						my @dimerres;
						push @dimerres, dimer_calcu($probe, $probe, "N");
						push @dimerres, dimer_calcu($F_primer, $probe, "Y");
						push @dimerres, dimer_calcu($R_primer, $probe, "Y");
						if(grep{$_ eq "nonpass"} @dimerres){
							$pb_Tm = 0;
						}
						if(($pb_Tm - $F_tm) <= 12 and ($pb_Tm - $F_tm) >= 8 and ($pb_Tm - $R_tm) <= 12 and ($pb_Tm - $R_tm) >= 8){
							if($cons_file){
								my $avg_score = sum(@cons_score[($F_start..$F_end, $R_end..$R_start, ($pb_start - $pb_len + 1)..($pb_start))])/(($F_end - $F_start+1)+($R_start - $R_end+1)+$pb_len);
								$primers_cons_sort{$avg_score}{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nConserve socre:\n".$avg_score."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
							else{
								$primers{(abs($F_tm - $R_tm))} .= ("=" x 60)."\nForward: ".$F_primer."\tTm = ".$F_tm."\n"."Reverse: ".$R_primer."\tTm = ".$R_tm."\n"."Probe: ".$probe."\tTm = ".$pb_Tm."\n";
							}
						}
					}
				}
			}
		}
	}
	open(OUT, ">", $analysis_path."Primers.tsv");
	if($cons_file){
		for my $avg_score (sort{$b<=>$a} keys %primers_cons_sort){
			for my $primer_score (sort{$a<=>$b} keys %{$primers_cons_sort{$avg_score}}){
				print OUT $primers_cons_sort{$avg_score}{$primer_score};
			}
		}
	}
	else{
		for my $primer_score (sort{$a<=>$b} keys %primers){
			print OUT $primers{$primer_score};
		}
	}
	close OUT;
}

sub sum{
	my @arr = @_;
	my $total = 0;
	for my $i (0..$#arr){
		$total += $arr[$i];
	}
	return $total;
}

sub base_count{
	my @seq = @_;
	my ($count_A, $count_T, $count_C, $count_G) = (0, 0, 0, 0);
	for my $i (0..$#seq){
		if($seq[$i] eq "A"){
			$count_A += 1;
		}
		elsif($seq[$i] eq "T"){
			$count_T += 1;
		}
		elsif($seq[$i] eq "C"){
			$count_C += 1;
		}
		elsif($seq[$i] eq "G"){
			$count_G += 1;
		}
		elsif($seq[$i] eq "R"){
			$count_A += 0.5;
			$count_G += 0.5;
		}
		elsif($seq[$i] eq "Y"){
			$count_C += 0.5;
			$count_T += 0.5;
		}
		elsif($seq[$i] eq "M"){
			$count_A += 0.5;
			$count_C += 0.5;
		}
		elsif($seq[$i] eq "K"){
			$count_G += 0.5;
			$count_T += 0.5;
		}
		elsif($seq[$i] eq "S"){
			$count_G += 0.5;
			$count_C += 0.5;
		}
		elsif($seq[$i] eq "W"){
			$count_A += 0.5;
			$count_T += 0.5;
		}
		elsif($seq[$i] eq "B"){
			$count_G += 1/3;
			$count_T += 1/3;
			$count_C += 1/3;
		}
		elsif($seq[$i] eq "V"){
			$count_G += 1/3;
			$count_A += 1/3;
			$count_C += 1/3;
		}
		elsif($seq[$i] eq "D"){
			$count_G += 1/3;
			$count_A += 1/3;
			$count_T += 1/3;
		}
		elsif($seq[$i] eq "H"){
			$count_A += 1/3;
			$count_C += 1/3;
			$count_T += 1/3;
		}
		elsif($seq[$i] eq "N"){
			$count_A += 0.25;
			$count_T += 0.25;
			$count_C += 0.25;
			$count_T += 0.25;
		}
	}
	return($count_A, $count_T, $count_C, $count_G);
}

sub pm_Tm_calcu {
	my @arr = @_;
	my @seq = split//,$arr[0];
	my ($count_A, $count_T, $count_C, $count_G) = base_count(@seq);
	my $Tm_value = (81.5 - 18.182380 + 0.41*($count_C + $count_G)*100/length($arr[0]) - 600/length($arr[0]));
	if(($count_C + $count_G)/length($arr[0]) >= 0.8 or ($count_C + $count_G)/length($arr[0]) <= 0.30){
		$Tm_value = 0;
	}
	if($seq[$#seq - 2] =~ m/[G|C]/ and $seq[$#seq - 1] =~ m/[G|C]/ and $seq[$#seq] =~ m/[G|C]/){
		$Tm_value = 0;
	}
	return($Tm_value);
}

sub pb_Tm_calcu{
	my @arr = @_;
	my @seq = split//,$arr[0];
	my ($count_A, $count_T, $count_C, $count_G) = base_count(@seq);
	my $Tm_value = (81.5 - 18.182380 + 0.41*($count_C + $count_G)*100/length($arr[0]) - 600/length($arr[0]));
	if($seq[0] eq "G"){
		$Tm_value = 0;
	}
	if($count_G > $count_C){
		$Tm_value = 0;
	}
	if(($count_C + $count_G)/length($arr[0]) >= 0.7 or ($count_C + $count_G)/length($arr[0]) <= 0.4){
		$Tm_value = 0;
	}
	if($arr[0] =~ m/AAAAAA/ or join("", @seq[($#seq - 3)..$#seq]) =~ m/GGG/ or join("", @seq[1..($#seq - 1)]) =~ m/CCC/ or join("", @seq) =~ m/GGGG/){
		$Tm_value = 0;
	}
	return($Tm_value);
}

sub dimer_calcu{
	my @match;
	my ($seq1, $seq2_raw, $three) = @_[0..2];
	my $seq2 = reverse($seq2_raw);
	$seq2 =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
	if($three eq "Y"){
		my $endseq = substr($seq1, (length($seq1) - 3), 3);
		if($seq2 =~ m/$endseq/){
			my ($count_A, $count_T, $count_C, $count_G) = base_count(split//,$endseq);
			if($count_C + $count_G >= 2){
				push @match, "strong";
			}
		}
	}
	my @arr1 = split//,$seq1;
	my @arr2 = split//,$seq2;
	my %score;
	my $max_score = 0;
	for my $i (0..$#arr1){
		for my $j (0..$#arr2){
			if($arr1[$i] eq $arr2[$j]){
				if(!$score{($i - 1)}{($j - 1)}){
					$score{$i}{$j} += 1;
				}
				else{
					$score{$i}{$j} = $score{($i - 1)}{($j - 1)} + 1;
				}
				if($max_score < $score{$i}{$j}){
					$max_score = $score{$i}{$j};
				}
			}
			else{
				$score{$i}{$j} = 0;
			}
		}
	}
	my @max_position;
	my $link;
	for my $i (0..$#arr1){
		for my $j (0..$#arr2){
			if($max_score == $score{$i}{$j}){
				$max_score = $score{$i}{$j};
				push @max_position, $i."-".$j;
			}
		}
	}
	for my $align_site (0..$#max_position){
		my @position = split/-/,$max_position[$align_site];
		my ($align_seq1, $align_seq2, $align_seq2_raw) = ($seq1, $seq2, $seq2);
		$align_seq2_raw =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
		my ($seq1_prefix, $seq1_suffix, $seq2_prefix, $seq2_suffix) = ($position[0], ($#arr1 - $position[0]), $position[1], ($#arr2 - $position[1]));
		if($seq1_prefix < $seq2_prefix){
			$align_seq1 = ("-" x ($seq2_prefix - $seq1_prefix)).$align_seq1;
		}
		elsif($seq1_prefix > $seq2_prefix){
			$align_seq2 = ("-" x ($seq1_prefix - $seq2_prefix)).$align_seq2;
			$align_seq2_raw = ("-" x ($seq1_prefix - $seq2_prefix)).$align_seq2_raw;
		}
		if($seq1_suffix < $seq2_suffix){
			$align_seq1 = $align_seq1.("-" x ($seq2_suffix - $seq1_suffix));
		}
		elsif($seq1_suffix > $seq2_suffix){
			$align_seq2 = $align_seq2.("-" x ($seq1_suffix - $seq2_suffix));
			$align_seq2_raw = $align_seq2_raw.("-" x ($seq1_suffix - $seq2_suffix));
		}
		my @alignseq1 = split//,$align_seq1;
		my @alignseq2 = split//,$align_seq2;
		for my $i (0..$#alignseq1){
			if($alignseq1[$i] eq $alignseq2[$i]){
				$link .= "|";
			}
			else{
				$link .= " ";
			}
		}
		#print OUT ("=" x 60)."\n";
		#print OUT $align_seq1."\n";
		#print OUT $link."\n";
		#print OUT $align_seq2_raw."\n";
	}
	my $link_num = $link =~ tr/|/|/;
	if($max_score >= 5 or ($max_score >= 4 and $link_num >= 12)){
		push @match, "strong";
	}
	if(grep{$_ eq "strong"} @match){
		return("nonpass");
	}
	else{
		return("pass");
	}
}
