


************** Bulk RNA-seq script


******** Fastp filtering of raw DDX43 data

nohup /home/planaria/anaconda3/bin/fastp -i C1_S1_R1_001.fastq -o C1_S1_R1_1_fastp.fastq -I C1_S1_R2_001.fastq -O /neoblast/Ankit/DDX43_KD_data/Fastp_new_results/C1_S1_R2_1_fastp.fastq -h C1_S1_1.html &> C1_S1_1_fastp_nohup.txt &

nohup /home/planaria/anaconda3/bin/fastp -i C2_S3_R1_001.fastq -o C2_S3_R1_2_fastp.fastq -I C2_S3_R2_001.fastq -O C2_S3_R2_2_fastp.fastq -h C2_S3_2.html &> C2_S3_2_fastp_nohup.txt &


nohup /home/planaria/anaconda3/bin/fastp -i KD1_S2_R1_001.fastq -o KD1_S2_R1_1_fastp.fastq -I KD1_S2_R2_001.fastq -O KD1_S2_R2_1_fastp.fastq -h KD1_S2_1.html &> KD1_S2_1_fastp_nohup.txt &

nohup /home/planaria/anaconda3/bin/fastp -i KD2_S4_R1_001.fastq -o KD2_S4_R1_2_fastp.fastq -I KD2_S4_R2_001.fastq -O KD2_S4_R2_2_fastp.fastq -h KD2_S4_2.html &> KD2_S4_2_fastp_nohup.txt &





**************** Mapping script with knockdown sample


nohup /softwares/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 10 --readFilesIn C1_S1_R1_1_fastp.fastq C1_S1_R2_1_fastp.fastq  --genomeDir /STAR_index_plana --outFileNamePrefix C1_S1_1_ --outReadsUnmapped Fastx &> C1_S1_1_STAR_nohup.txt &

nohup /softwares/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 10 --readFilesIn C2_S3_R1_2_fastp.fastq C2_S3_R2_2_fastp.fastq  --genomeDir /STAR_index_plana --outFileNamePrefix C2_S3_2_ --outReadsUnmapped Fastx &> C2_S3_2_STAR_nohup.txt &


nohup /softwares/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 10 --readFilesIn KD1_S2_R1_1_fastp.fastq KD1_S2_R2_1_fastp.fastq  --genomeDir /STAR_index_plana --outFileNamePrefix KD1_S2_1_ --outReadsUnmapped Fastx &> KD1_S2_1_STAR_nohup.txt &

nohup /softwares/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 10 --readFilesIn KD2_S4_R1_2_fastp.fastq KD2_S4_R2_2_fastp.fastq --genomeDir /STAR_index_plana --outFileNamePrefix KD2_S4_2_ --outReadsUnmapped Fastx &> KD2_S4_2_STAR_nohup.txt &




********* Counting script


nohup /softwares/subread-1.5.2-Linux-x86_64/bin/featureCounts -p -t exon -a smes_v2_repeatfilt_SMESG.gtf -o DDX43_S4_1_Aligned.txt DDX43_S4_1_Aligned.out.sam &> DDX43_S4_1_nohup.txt &

nohup /softwares/subread-1.5.2-Linux-x86_64/bin/featureCounts -p -t exon -a smes_v2_repeatfilt_SMESG.gtf -o DDX43_S5_2_Aligned.txt DDX43_S5_2_Aligned.out.sam &> DDX43_S5_2_nohup.txt &


nohup /softwares/subread-1.5.2-Linux-x86_64/bin/featureCounts -p -t exon -a smes_v2_repeatfilt_SMESG.gtf -o GFP1_S1_1_Aligned.txt GFP1_S1_1_Aligned.out.sam &> GFP1_S1_1_nohup.txt &

nohup /softwares/subread-1.5.2-Linux-x86_64/bin/featureCounts -p -t exon -a smes_v2_repeatfilt_SMESG.gtf -o GFP2_S2_2_Aligned.txt GFP2_S2_2_Aligned.out.sam &> GFP2_S2_2_nohup.txt &





********** Differential script


data <- read.table("DDX43-KD_vs_GFP-Control_2-replicates.txt",header=T,row.names=1)
groups <- factor(c(rep("TGroup",2),rep("CGroup",2)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds,independentFiltering=F)
rld <- rlogTransformation(dds)
write.table(res,file="DDX43-KD_vs_GFP-Control_two_replicates_DE_all_genes.txt",sep="\t",quote=F)