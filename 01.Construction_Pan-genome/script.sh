#This pipeline used to construct the duck linear pan-genome using a combination approach of Psvcp and PPsPCP pipelines.
#Five published duck genomes were used to construct the duck pan-genome, which consists of three Pekin duck genomes, one Shaoxing duck genome, and mallard duck genome.

#01.Colinear analysis between four query genomes and reference genomes

#01.1 make the bed and gff3 files
for i in CAUPekin2 duckbaserefseqv4 ASM874695v11 CAULaying ZJU
do
echo python -m jcvi.formats.gff bed --type=mRNA --key=Name ${gff_dir}/${i}.gff3 -o ${bed_dir}/${i}.bed
done

#01.2 extract the cds sequence and aa sequence
for i in CAUPekin2 duckbaserefseqv4 ASM874695v1 CAULaying ZJU
do
gffread ${gff_dir}/${i}.gff3 -g ${fa_dir}/${i}.fa -x  ${pep_dir}/${i}.cds -y ${pep_dir}/${i}.pep
done

#01.3 colinear analysis
for i in CAUPekin2 duckbaserefseqv4 ASM874695v1 CAULaying
do
mkdir ${main_dir}/${i} && cd ${main_dir}/${i}
ln -s ${pep_dir}/${i}.pep ./
ln -s ${pep_dir}/ZJU.pep ./
ln -s ${bed_dir}/${i}.bed ./
ln -s ${bed_dir}/ZJU.bed ./
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names ZJU ${i}
done

#01.4 Generate the colinear plot
for i in CAUPekin2 duckbaserefseqv4 ASM874695v1 CAULaying
do
cd ${main_dir}/${i}
python -m jcvi.compara.synteny screen --minspan=30 --simple ZJU.$i.anchors ZJU.$i.anchors.new
python -m jcvi.graphics.karyotype --font Arial --style dark --figsize 10x6 seqids layout 
mv karyotype.pdf ZJU.${i}.karyotype.pdf
done

#02. Constructing the Pan-genome.1 using Psvcp pipeline
#First step, we used Psvcp pipeline to ingrate insertions larger than 50 bp into reference genome. Briefly, duckbase.refseq.v4 genome was aligned to the initial reference genome (ZJU1.0), insertions longer than 50 bp were then identified and placed in the ZJU1.0. This process was further iterated by CAU_Pekin2.0, CAU_Laying_1.0, and ASM8764695v1 in an order of increasing phylogenetic distance to ZJU1.0 and thus the Pan-genome.1 was generated.
#According to the direction consistency of chromosome, we handled this step each chromosome by chromosome

#02.1 round1-duckbaserefseqv4
cat the Chromosome_colinear.list|while read chr1 chr2
do 
mkdir ${chr1} && cd ${chr1}
seqkit grep -p ${chr1} ZJU.revised.fa >${chr1}/${chr1}.fa
seqkit grep -p ${chr2} duckbaserefseqv4.fa >${chr1}/${chr2}.fa
bash Genome_construct_Pangenome_Psvcp.sh ${chr1}/${chr1}.fa ${chr1}/${chr2}.fa > jobs.sh
bash jobs.sh
done
#02.2 round1-CAUPekin2
scripts similar to 02.1 
#02.3 round1-CAULaying
scripts similar to 02.1
#02.4 round1-ASM874695v1
scripts similar to 02.1

#03. Constructing the Pan-genome.1 using Psvcp pipeline
#Second step, query genomes were aligned to Pan-genome.1, respectively, while novel contigs longer than 500 bp were retained after removing redundancy. Novel contigs and Pan-genome.1 were merged into the final duck pan-genome
#03.1 duckbaserefseqv4.fa
~/biosoft/ppsPCP/bin/make_pan.pl --ref pan1.fa --ref_anno pan1.gff --query duckbaserefseqv4.fa --query_anno duckbaserefseqv4.CAUwild.gff3 --thread 40 

#03.2 CAUPekin2.fa
~/biosoft/ppsPCP/bin/make_pan.pl --ref pan1.fa --ref_anno pan1.gff --query CAUPekin2.fa --query_anno CAUPekin2.CAUwild.gff3 --thread 40

#03.3 CAULaying
~/biosoft/ppsPCP/bin/make_pan.pl --ref pan1.fa --ref_anno pan1.gff --query CAULaying.fa --query_anno CAULaying.CAUwild.gff3 --thread 40

#03.4 ASM874695v1
~/biosoft/ppsPCP/bin/make_pan.pl --ref pan1.fa --ref_anno pan1.gff --query ASM874695v1.fa --query_anno ASM874695v1.gff3 --thread 40

#04 Filtering and removing redundancy
# collecting each result and filtering (absence_filtered.fa)
cat duckbaserefseqv4_absence_filtered.rename.fa CAUPekin2_absence_filtered.rename.fa CAULaying_absence_filtered.rename.fa ASM874695v1_absence_filtered.rename.fa >novel.clean.merge.fa
bioawk -c fastx 'length($seq)>500 {print ">"$name; print $seq }' novel.clean.merge.fa >novel.clean.merge.gt500.fa
~/biosoft/cdhit-4.8.1/cd-hit-est -M 0 -i novel.clean.merge.gt500.fa -c 0.9 -aS 0.8 -d 0 -sf 1 -T 100 -M 0 -o novel.clean.merge.gt500.cdhit.out


