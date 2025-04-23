# Comp_Trial_Axelle_Exp
First analysis comparing Axelle's trial with the actual experiment


![Bacteria Doing Bioinformatics](https://github.com/ivandamg/Comp_Trial_Axelle_Exp/blob/main/b114de86-4949-4363-b6fc-747e29d75356.jpeg)

# 1. Change filename to discard unique names

        mv WT1_L1_R1_001_gGg9DyMWQgI6.fastq.gz WT1_L1_R1.fastq.gz
        mv WT1_L1_R2_001_GylKkOsoETPC.fastq.gz WT1_L1_R2.fastq.gz
        mv WT1_L2_R1_001_4BOXccxl6MmC.fastq.gz WT1_L2_R1.fastq.gz
        mv WT1_L2_R2_001_jiWnaH6npIIN.fastq.gz WT1_L2_R2.fastq.gz



# 2. Extract UMI

        mkdir 02_UMIExtract
        for FILE in $(ls *_R1.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_UMI --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=UMI_$(echo $FILE | cut -d'_' -f1,2).out --error=UMI_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/01_RawData;module load UMI-tools/1.0.1-foss-2021a; umi_tools extract --bc-pattern=NNNNNNNNNNNN --stdin=$FILE --stdout=../02_UMIExtract/$(echo $FILE | cut -d'_' -f1,2)_extracted_R1.fastq.gz --read2-in=$(echo $FILE | cut -d'_' -f1,2)_R2.fastq.gz --read2-out=../02_UMIExtract/$(echo $FILE | cut -d'_' -f1,2)_extracted_R2.fastq.gz"; sleep  1; done

# 3. Check quality

        for FILE in $(ls *.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2,3)fastQC --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=$(echo $FILE | cut -d'_' -f1,2_3)_fastQC.out --error=$(echo $FILE | cut -d'_' -f1,2_3)_fastQC.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/01* ; module load FastQC; fastqc -t 4 $FILE"; sleep  1; done

# 4. Concatenate L1 and L2

        mkdir ../03_Concatenated
        cat WT1_L1_extracted_R1.fastq.gz WT1_L2_extracted_R1.fastq.gz > ../03_Concatenated/WT1_extracted_R1.fastq.gz
        cat WT1_L1_extracted_R2.fastq.gz WT1_L2_extracted_R2.fastq.gz > ../03_Concatenated/WT1_extracted_R2.fastq.gz

# 4. fastp trimming and filter reads

        mkdir ../04_TrimmedData/
    for FILE in $(ls *extracted_R1.fastq.gz); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-08:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=4 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/03_Concatenated ; module load FastQC; module load fastp/0.23.4-GCC-10.3.0; fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_extracted_R2.fastq.gz --out1 ../04_TrimmedData/$(echo $FILE | cut -d'_' -f1)_R1_trimmed.fastq.gz --out2 ../04_TrimmedData/$(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz -h ../04_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../04_TrimmedData/$(echo $FILE | cut -d'_' -f1)_R1_trimmed.fastq.gz; fastqc -t 4 ../04_TrimmedData/$(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz"; sleep  1; done


# 5. index genome medicago

    sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/test3/00_ReferenceGenomes; module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/test3/00_ReferenceGenomes --genomeFastaFiles GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna --sjdbGTFfile GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff --sjdbOverhang 99 --genomeSAindexNbases 10"


# 6. Map reads to medicago

        mkdir ../05_MappedMedicago
        mkdir ../06_UnmappedReads
       for FILE in $(ls *_R1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1)_1STAR --time=1-12:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_STAR.out --error=$(echo $FILE | cut -d'_' -f1)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/04_TrimmedData; STAR --runThreadN 8 --genomeDir /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/00_References --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1)_R2_trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix $(echo $FILE | cut -d'_' -f1)_MappedMedicago --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --limitBAMsortRAM 5919206202 --limitOutSJcollapsed 5000000; mv $(echo $FILE | cut -d'_' -f1)_MappedMedicago.bam ../05_MappedMedicago/$(echo $FILE | cut -d'_' -f1)_MappedMedicago.bam ; mv $(echo $FILE | cut -d'_' -f1)_MappedMedicagoUnmapped.out.mate1 ../06_UnmappedReads/$(echo $FILE | cut -d'_' -f1)_UnmappedMedicago_R1.fastq; mv $(echo $FILE | cut -d'_' -f1)_MappedMedicagoUnmapped.out.mate2 ../06_UnmappedReads/$(echo $FILE | cut -d'_' -f1)_UnmappedMedicago_R2.fastq "; sleep  1; done

# 7. gzip unmapped reads

    for FILE in $(ls *.fastq); do echo $FILE; sbatch --partition=pshort_el8 --job-name=gzip$(echo $FILE | cut -d'_' -f1) --time=0-02:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_gzip.out --error=$(echo $FILE | cut -d'_' -f1)_gzip.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/12_dualRNAseqv3/06_UnmappedReads/ ; gzip $FILE"; done


# 5. divide gff into proteins and others

    cat GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff |grep "gene_biotype=p" > Medicago_proteins.gff

    cat GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gff |grep -v "gene_biotype=p" > Medicago_others.gff


# 6. Counts reads with feature counts -M --primary (count also multimapped reads but only once..)

### 6a. count reads to others tRNA and rRNA,etc

    for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/Medicago_others.gff  -o CountsTableMedicago_UniqueMultiple_Others_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done

### 6b.   count reads to proteins

       for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/Medicago_proteins.gff  -o CountsTableMedicago_UniqueMultiple_Proteins_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done


# 7. index genome Rhizobia

        sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia; module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia --genomeFastaFiles GCF_037023865.1_ASM3702386v1_genomic.fna --sjdbGTFfile GCF_037023865.1_ASM3702386v1_genomic.gff --sjdbOverhang 99 --genomeSAindexNbases 10"




# 8. Map unmapped reads to rhizobia

           for FILE in $(ls WT*_UnmappedMedicago_R1.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=1-12:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/04_UnmappedReads; STAR --runThreadN 8 --genomeDir /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2)_UnmappedMedicago_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $(echo $FILE | cut -d'_' -f1,2)_MappedRhizobia --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5919206202 --limitOutSJcollapsed 5000000"; sleep  1; done

# 9. divide gff into proteins and others

        cat GCF_037023865.1_ASM3702386v1_genomic.gff |grep "gene_biotype=p" > Rhizobium_proteins.gff

        cat GCF_037023865.1_ASM3702386v1_genomic.gff |grep -v "gene_biotype=p" > Rhizobium_others.gff



# 10. Counts also multiple mapped reads once
###  10a. count  reads to others tRNA and rRNA

        for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia/Rhizobium_others.gff  -o CountsTableRhizobia_UniqueMultiple_Others_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done

###  10b.  count  reads to proteins

       for FILE in $(ls *.bam ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=FC_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=64G --ntasks=1 --cpus-per-task=1 --output=FC_$(echo $FILE | cut -d'_' -f1,2).out --error=FC_$(echo $FILE | cut -d'_' -f1,2).error --mail-type=END,FAIL --wrap "module load Subread; featureCounts -p -M --primary --countReadPairs -t gene -g ID -a /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/00_ReferenceGenomes/01_Rhizobia/Rhizobium_proteins.gff  -o CountsTableRhizobia_UniqueMultiple_Proteins_$(echo $FILE | cut -d'_' -f1,2).txt $FILE -T 8"; sleep  1; done

# 11. IGView viz create bai file for vizualisation

       for FILE in $(ls *.bam); do echo $FILE; sbatch --partition=pshort_el8 --job-name=index_$(echo $FILE | cut -d'_' -f1,2) --time=0-02:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_index.out --error=$(echo $FILE | cut -d'_' -f1,2)_index.error --mail-type=END,FAIL --wrap "module load SAMtools/1.13-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/03_TrimmedData; samtools index $FILE"; done 

