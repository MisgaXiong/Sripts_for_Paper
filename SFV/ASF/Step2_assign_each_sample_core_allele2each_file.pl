#!/usr/bin/perl
use Bio::SeqIO;

=head1 Name

	Step2_assign_each_sample_core_allele2each_file.pl 

=head1 Description

	After do the Step1_merge_all_allele2total_fasta.r, a total.fasta will be output. 
	A fa_name.txt file is required. It contains all the names from genomes.fa, such as AY261361.1.fasta, MH766894.1.fasta...
	for i in `cat fa_name.txt`; do perl Step2_assign_each_sample_core_allele2each_file.pl total.fasta ${i} > ./result/${i}; done.
    The output file from this perl script need to be checked, for example, the "_" of NC_001659.2.fasta need to be removed.  
    sed -i 's/\NC_//g' NC_001659.2.fasta

=head2 Usage

	perl Step2_assign_each_sample_core_allele2each_file.pl total.fasta > ./path/yourfile

=head2 Author

	Xiong Dongyan

=cut


my $total_seq = $ARGV[0];
my $exctract_id = $ARGV[1];
my $fasta_total_seq = Bio::SeqIO->new(
	-file => $total_seq,
	-format => 'fasta',
);
my %Each_seq;
while (my $seq = $fasta_total_seq->next_seq()) {
	my $id = $seq -> id;
	$id = ">".$id;
	my $seq_seq = $seq -> seq;
	if($id =~ m/$exctract_id/){$Each_seq{$id} = "\n".$seq_seq."\n"}
	else{%Each_seq = %Each_seq};
};
print %Each_seq;
