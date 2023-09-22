#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;

my $makeblastdb = "/home/dell/software/ncbi-blast/bin/makeblastdb";
my $blastn = "/home/dell/software/ncbi-blast/bin/blastn";
my $blastp = "/home/dell/software/ncbi-blast/bin/blastp";
my $blastx = "/home/dell/software/ncbi-blast/bin/blastx";
my $tblastn = "/home/dell/software/ncbi-blast/bin/tblastn";
my $mafft = "/home/dell/miniconda3/envs/wgs/bin/mafft";

my $analysis_path = $ARGV[0];
my ($model, $fasta) = @ARGV[1..2];
my $direact = $ARGV[3];
my $codontable = $ARGV[3];
my $extract_model = $ARGV[3];
my $sub_extract = $ARGV[4];
my $motif_seq = $ARGV[5];
my $motif_type = $ARGV[6];
my $extract_type = $ARGV[7];
my $batch_file = $ARGV[5];

if($model eq "reversecomplent"){
	my @out;
	my $sequence = Bio::SeqIO -> new(-file => $fasta, -format => 'fasta');
	while(my $seq = $sequence -> next_seq){
		my ($id, $desc, $seq_seq) = ($seq -> id, $seq -> desc, $seq -> seq);
		$seq_seq = uc($seq_seq);
		if($direact eq "reverse"){
			$seq_seq = reverse($seq_seq);
			push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
		}
		elsif($direact eq "complementation"){
			$seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
			push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
		}
		elsif($direact eq "reverse-complement"){
			$seq_seq = reverse($seq_seq);
			$seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
			push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
		}
	}
	print @out;
}
elsif($model eq "seqtranslate"){
	my @out;
	my $sequence = Bio::SeqIO -> new(-file => $fasta, -format => 'fasta');
	while(my $seq = $sequence -> next_seq){
		my ($id, $desc, $seq_seq) = ($seq -> id, $seq -> desc, $seq -> seq);
		my $seq_aa;
		$seq_seq = uc($seq_seq);
		if(length($seq_seq) <= 3){
			$seq_aa = "please input the sequence with at least 3bp";
		}
		else{
			my $seq_obj = Bio::Seq -> new(-seq => $seq_seq, -alphabet => 'dna');
			$seq_aa = $seq_obj -> translate(-codontable_id => $codontable,) -> seq;
			$seq_aa =~ s/\*//g;
		}
		push @out, ">".$id." ".$desc."\n".$seq_aa."\n";
	}
	print @out;
}
elsif($model eq "seqextract"){
	my @out;
	if($extract_model eq "motif-extract"){
		my $command;
		if($motif_type eq "nucleicacid" and $extract_type eq "nucleicacid"){
			$command = $makeblastdb." -in ".$analysis_path.$motif_seq." -dbtype nucl -logfile ".$analysis_path."makedb.log";
			system(qq($command));
			$command = $blastn." -query ".$analysis_path.$fasta." -db ".$analysis_path.$motif_seq." -num_threads 12 -evalue 1e-8 -max_hsps 1 -outfmt 6 -out ".$analysis_path."blast_res.m6";
			system(qq($command));
			push @out, extract_seq($analysis_path."blast_res.m6", $analysis_path.$fasta);
		}
		elsif($motif_type eq "protein" and $extract_type eq "nucleicacid"){
			$command = $makeblastdb." -in ".$analysis_path.$motif_seq." -dbtype prot -logfile ".$analysis_path."makedb.log";
			system(qq($command));
			$command = $blastx." -query ".$analysis_path.$fasta." -db ".$analysis_path.$motif_seq." -num_threads 12 -evalue 1e-8 -max_hsps 1 -outfmt 6 -out ".$analysis_path."blast_res.m6";
			system(qq($command));
			push @out, extract_seq($analysis_path."blast_res.m6", $analysis_path.$fasta);
		}
		elsif($motif_type eq "protein" and $extract_type eq "protein"){
			$command = $makeblastdb." -in ".$analysis_path.$motif_seq." -dbtype prot -logfile ".$analysis_path."makedb.log";
			system(qq($command));
			$command = $blastp." -query ".$analysis_path.$fasta." -db ".$analysis_path.$motif_seq." -num_threads 12 -evalue 1e-8 -max_hsps 1 -outfmt 6 -out ".$analysis_path."blast_res.m6";
			system(qq($command));
			push @out, extract_seq($analysis_path."blast_res.m6", $analysis_path.$fasta);
		}
		elsif($motif_type eq "nucleicacid" and $extract_type eq "protein"){
			$command = $makeblastdb." -in ".$analysis_path.$motif_seq." -dbtype nucl -logfile ".$analysis_path."makedb.log";
			system(qq($command));
			$command = $tblastn." -query ".$analysis_path.$fasta." -db ".$analysis_path.$motif_seq." -num_threads 12 -evalue 1e-8 -max_hsps 1 -outfmt 6 -out ".$analysis_path."blast_res.m6";
			system(qq($command));
			push @out, extract_seq($analysis_path."blast_res.m6", $analysis_path.$fasta);
		}
		print @out;
	}
	elsif($extract_model eq "site-extract"){
		my %site_hash;
		open(FILE, "<", $analysis_path.$batch_file);
		while(<FILE>){
			chomp;
			my @line = split/\t/;
			@{$site_hash{$line[0]}} = ($line[1], $line[2]);
		}
		my $sequence = Bio::SeqIO -> new(-file => $analysis_path.$fasta, -format => 'fasta');
		while(my $seq = $sequence -> next_seq){
			my ($id, $desc, $seq_seq) = ($seq -> id, $seq -> desc, $seq -> seq);
			my $seq_seq_seq = substr($seq_seq, (@{$site_hash{$id}}[0] - 1), (@{$site_hash{$id}}[1] - @{$site_hash{$id}}[0] + 1));
			my ($start, $end) = @{$site_hash{$id}}[0..1];
			push @out, ">".$id." ".$desc." from ".$start." to ".$end."\n".$seq_seq_seq."\n";
		}
		print @out;
	}
}
elsif($model eq "seqalignment"){
	my @out;
	my $sequence = Bio::SeqIO -> new(-file => $fasta, -format => 'fasta');
	my ($number, $kmer) = (0, 17);
	my %ref_kmer;
	while(my $seq = $sequence -> next_seq){
		$number += 1;
		my ($id, $desc, $seq_seq, $seq_len) = ($seq -> id, $seq -> desc, $seq -> seq, $seq -> length);
		my ($match, $match_rc) = (0, 0);
		for my $i (0..($seq_len - $kmer)){
			my $kmer_seq = substr($seq_seq, $i, $kmer);
			if($number == 1){
				$ref_kmer{$kmer_seq} += 1;
			}
			else{
				my $kmer_seq_rc = reverse($kmer_seq);
				$kmer_seq_rc =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
				if($ref_kmer{$kmer_seq}){
					$match += 1;
				}
				if($ref_kmer{$kmer_seq_rc}){
					$match_rc += 1;
				}
			}
		}
		if($match >= $match_rc){
			push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
		}
		else{
			$seq_seq = reverse($seq_seq);
			$seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
			push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
		}
	}
	open(OUT, ">", $analysis_path."sequence_prepar.fasta");
	print OUT @out;
	close OUT;
	my $command = $mafft." --thread 12 --quiet ".$analysis_path."sequence_prepar.fasta > ".$analysis_path."sequence_aln.fasta";
	system(qq($command));
	my $align = Bio::SeqIO -> new(-file => $analysis_path."sequence_aln.fasta", -format => 'fasta');
	undef @out;
	while(my $seq = $align -> next_seq){
		my ($id, $desc, $seq_seq) = ($seq -> id, $seq -> desc, uc($seq -> seq));
		push @out, ">".$id." ".$desc."\n".$seq_seq."\n";
	}
	unlink $analysis_path."sequence_aln.fasta";
	print @out;
}

sub extract_seq{
	my $out;
	my ($blast_file, $file_path) = ($_[0], $_[1]);
	my %extract_hash;
	open(FILE, "<", $blast_file);
	while(<FILE>){
		chomp;
		my @line = split/\t/;
		if($sub_extract eq "TRUE"){
			@{$extract_hash{$line[0]}} = ($line[6], $line[7]);
		}
		else{
			$extract_hash{$line[0]} = 1;
		}
	}
	close FILE;
	my $sequence = Bio::SeqIO -> new(-file => $file_path, -format => 'fasta');
	while(my $seq = $sequence -> next_seq){
		my ($id, $desc, $seq_seq) = ($seq -> id, $seq -> desc, $seq -> seq);
		$seq_seq = uc($seq_seq);
		if($sub_extract eq "TRUE"){
			if($extract_hash{$id}){
				if(@{$extract_hash{$id}}[0] <= @{$extract_hash{$id}}[1]){
					$out .= ">".$id." ".$desc." aligned region\n".substr($seq_seq, (@{$extract_hash{$id}}[0] - 1), (@{$extract_hash{$id}}[1] - @{$extract_hash{$id}}[0] + 1))."\n";
				}
				else{
					my $seq_seq_seq = substr($seq_seq, (@{$extract_hash{$id}}[1] - 1), (@{$extract_hash{$id}}[0] - @{$extract_hash{$id}}[1] + 1));
					$seq_seq_seq = reverse($seq_seq_seq);
					$seq_seq_seq =~ tr/ATCGRYMKSWBVDHN/TAGCYRKMSWVHDN/;
					$out .= ">".$id." ".$desc." aligned region\n".$seq_seq_seq."\n";
				}
			}
		}
		else{
			if($extract_hash{$id}){
				$out .= ">".$id." ".$desc." with target sequence\n".$seq_seq."\n";
			}
		}
	}
	return $out;
}
