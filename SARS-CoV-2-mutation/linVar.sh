#!/usr/bin/bash
cat SARS-CoV-2_Ref.fasta Seq.fasta > SARS-CoV-2_RefvsSeq.fasta;
mafft --thread 12 SARS-CoV-2_RefvsSeq.fasta > SARS-CoV-2_RefvsSeq_maf.fasta;
perl gcgenome.pl -refname MN908947.3 -aligned SARS-CoV-2_RefvsSeq_maf.fasta;
perl ATCG_SNP.pl snp.tsv > ATCG_snp.tsv;
perl SNP2AAMut.pl SARS-CoV-2_Ref.gff3 SARS-CoV-2_Ref.fasta ATCG_snp.tsv > SNPAA.tsv;
GISAID="$(cat Seq.fasta | head -n 1 | rev | cut -d \| -f 2 | rev)";
echo ${GISAID};
LINEAGE="$(cat Seq.fasta | head -n 1 | rev | cut -d \| -f 1 | rev)";
perl strainVar.pl SARS-CoV-2_Ref.gff3 ${GISAID} ${LINEAGE};
cat Strain.var >> Total_strainVar.tsv;
rm -rf Seq.fasta && rm -rf SARS-CoV-2_RefvsSeq.fasta && rm -rf SARS-CoV-2_RefvsSeq_maf.fasta;
rm -rf snp.tsv && rm -rf insertion.tsv && rm -rf deletion.tsv && rm -rf ATCG_snp.tsv && rm -rf SNPAA.tsv && rm -rf Strain.var; 
