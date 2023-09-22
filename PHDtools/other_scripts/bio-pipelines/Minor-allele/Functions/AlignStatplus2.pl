#!/usr/bin/perl -w

=head1 Description
	Four base alignment results statistics
=head1 Usage
	perl AlignStat.pl -genome <path to genome file> -input <Readable bam file> -output <Output file name>
=head1 Parameters
	-genome			[str]	Input the completed path and the genome fasta file, required
	-input/in/i		[str]	Input the name of the readable bam file, required
	-output/out/o	[str]	Input the name of the output file, required
	-PCRdel			[str]	Delete the PCR repeat sequence, default = FALSE, optional
	-SecAlign		[str]	Delete the secondary alignment record, default = FALSE, optional
	-h/-help		[str]	print help
=head1 Author
	Dongyan Xiong
=head1 Edit Time
	2019.07.24 17:17 0.0.1
	2021.09.04 21:21 1.0.0
	2021.09.05 17:56 1.1.0
=cut

use warnings;
use strict;
use POSIX;
use Bio::SeqIO;
use Getopt::Long;

my $genomef = "";
my $input = "";
my $output = "";
my $PCRdel = "FALSE";
my $SecAlign = "FALSE";
my $help = "";
GetOptions(
	'genome=s' => \$genomef,
	'input|in|i=s' => \$input,
	'output|out|o=s' => \$output,
	'PCRdel:s' => \$PCRdel,
	'SecAlign:s' => \$SecAlign,
	'h|help:s' => \$help,
);
die `pod2text $0` if((!$genomef) or (!$input) or (!$output) or ($help));

my %genome_align;
my %flag_hash = (
1 => "Pair End", 2 => "Normal Alignment", 3 => "Not Alignment", 4 => "PE Another not Alignment", 5 => "Reversecomplementary Alignment",
6 => "PE Another Reversecomplementary Alignment", 7 => "PE Read1", 8 => "PE Read2", 9 => "Secondaty Alignment", 10 => "Below the Threshold Value",
11 => "PCR Repeat", 12 => "Supplementary alignment",
);

my %base_map = ("A" => 0, "T" => 1, "C" => 2, "G" => 3, "-" => 4, "I" => 5);
my $genome_file = Bio::SeqIO -> new(
	-file => $genomef,
	-format => 'fasta',
);
my %chr_map;
while(my $fasta = $genome_file -> next_seq){
	my $id = $fasta -> id;
	my $seqlen = $fasta -> length;
	$chr_map{$id} = $seqlen;
}

open(OUT, ">", $output);
my @header;
push @header, join("\t", ("Chr", "Pos", "A", "T", "C", "G", "-", "I"));
push @header, "\n";
print OUT @header;
close OUT;

open(BAMFILE, "<", $input);
my $line_num = 0;
my $current_chr;
while(<BAMFILE>){
	$line_num += 1;
	chomp;
	my @line = split/\t/,$_;
	if($line[5] !~ m/\*/){
		my ($flag, $chr, $start_pos, $cigar_value, $sequence) = ($line[1], $line[2], $line[3], $line[5], uc($line[9]));
		if($line_num == 1){
			$current_chr = $chr;
			for my $i (1..$chr_map{$current_chr}){
				@{$genome_align{$i}} = (0, 0, 0, 0, 0, 0); #(A, T, C, G, -, I)
			}
		}
		elsif($chr ne $current_chr){
			open(OUT, ">>", $output);
			my @out;
			for my $pos (sort{$a <=> $b} keys %genome_align){
				push @out, join("\t", ($current_chr, $pos, @{$genome_align{$pos}}));
				push @out, "\n";
			}
			print OUT @out;
			close OUT;
			undef @out;
			undef %genome_align;
			$current_chr = $chr;
			for my $i (1..$chr_map{$current_chr}){
				@{$genome_align{$i}} = (0, 0, 0, 0, 0, 0); #(A, T, C, G, -, I)
			}
		}
		my $flag_res = flag_trans($flag);
		if($flag_res eq "TRUE"){
			my @sequence = split//,$sequence;
			my @cigar_dig = split/\D/,$cigar_value;
			my @cigar_let = split/\d+/,$cigar_value;
			@cigar_let = @cigar_let[1..$#cigar_let];
			my $read_start = 1;
			for my $i (0..$#cigar_dig){
				if($i == 0){
					if($cigar_let[$i] eq "H"){
						next;
					}
					elsif($cigar_let[$i] eq "S"){
						$read_start += $cigar_dig[$i];
					}
					elsif($cigar_let[$i] eq "I"){
						@{$genome_align{($start_pos - 1)}}[$base_map{"I"}] += 1;
						$read_start += $cigar_dig[$i];
					}
					elsif($cigar_let[$i] eq "M"){
						for my $read_site (1..$cigar_dig[$i]){
							if($sequence[($read_site - 1)] !~ m/[A|T|C|G]/){
								next;
							}
							else{
								@{$genome_align{($start_pos + ($read_site - 1))}}[$base_map{$sequence[($read_site - 1)]}] += 1;
							}
						}
						$read_start += $cigar_dig[$i];
						$start_pos += $cigar_dig[$i];
					}
				}
				else{
					if($cigar_let[$i] eq "H"){
						next;
					}
					elsif($cigar_let[$i] eq "S"){
						$read_start += $cigar_dig[$i];
					}
					elsif($cigar_let[$i] eq "I"){
						@{$genome_align{($start_pos - 1)}}[$base_map{"I"}] += 1;
						$read_start += $cigar_dig[$i];
					}
					elsif($cigar_let[$i] eq "D"){
						for my $read_site (1..$cigar_dig[$i]){
							@{$genome_align{($start_pos + ($read_site - 1))}}[$base_map{"-"}] += 1;
						}
						$start_pos += $cigar_dig[$i];
					}
					elsif($cigar_let[$i] eq "M"){
						for my $read_site (1..$cigar_dig[$i]){
							if($sequence[($read_start + ($read_site - 1)) - 1] !~ m/[A|T|C|G]/){
								next;
							}
							else{
								@{$genome_align{($start_pos + ($read_site - 1))}}[$base_map{$sequence[($read_start + ($read_site - 1)) - 1]}] += 1;
							}
						}
						$read_start += $cigar_dig[$i];
						$start_pos += $cigar_dig[$i];
					}
				}
			}
		}
	}
}
close BAMFILE;
print "BAM file has been successfully processed!\n";
open(OUT, ">>", $output);
my @out;
for my $pos (sort{$a <=> $b} keys %genome_align){
	push @out, join("\t", ($current_chr, $pos, @{$genome_align{$pos}}));
	push @out, "\n";
}
print OUT @out;
close OUT;

sub flag_trans{
	my ($num_10, $num_2) = ($_[0], "");
	while($num_10 > 0){
		$num_2 = ($num_10 % 2).$num_2;
		$num_10 = floor($num_10 / 2);
		if($num_10 == 1){
			$num_2 = "1".$num_2;
			last;
		}
	}
	my @num2_arr = split//,reverse($num_2); #
	my $result = "TRUE";
	if($PCRdel =~ m/T/){
		unless((!$num2_arr[(11 - 1)] or $num2_arr[(11 - 1)] eq "0")){
			$result = "FALSE";
		}
	}
	if($SecAlign =~ m/T/){
		unless((!$num2_arr[(9 - 1)] or $num2_arr[(9 - 1)] eq "0")){
			$result = "FALSE";
		}
	}
	return $result;
}
