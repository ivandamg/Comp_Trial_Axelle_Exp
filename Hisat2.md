2. align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html

           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H2Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap2_HiSat2index.log --error=Hap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/TestHiSat; module load HISAT2; hisat2-build GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna Medicago_hisat_index -p 12"

   2. Map reads to genome, sort and compress

              # run from /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/03_TrimmedData
            for FILE in $(ls *_*_1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=3-08:00:00 --mem-per-cpu=128G --ntasks=16 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/comp_trial_Axelle/03_TrimmedData; module load HISAT2/2.2.1-gompi-2021a; hisat2 --phred33 --dta -x /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/TestHiSat/Medicago_hisat_index -1 $FILE -2 $(echo $FILE | cut -d'_' -f1,2)_2_trimmed.fastq.gz -S /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/TestHiSat/$(echo $FILE | cut -d'_' -f1,2)_MappedMedicago.sam --un-conc /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/TestHiSat/$(echo $FILE | cut -d'_' -f1,2)_UnmapMedicago -p 48"; done

3. Convert sam to bam

           for FILE in $(ls *_*_1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pibu_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)_1STAR --time=1-08:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR; cd /data/projects/p495_SinorhizobiumMeliloti/11_dualRNAseqv2/TestHiSat/;module load SAMtools/1.13-GCC-10.3.0; samtools view --threads 48 -b -o $(echo $FILE | cut -d'_' -f1,2)_MappedMedicago.bam $(echo $FILE | cut -d'_' -f1,2)_MappedMedicago.sam; samtools sort -o $(echo $FILE | cut -d'_' -f1,2)_MappedMedicago_sorted.bam -T temp_$(echo $FILE | cut -d'_' -f1,2) --threads 48 $(echo $FILE | cut -d'_' -f1,2)_MappedMedicago.bam"; done 

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/40_S_spinosum_FinalFinal/03_BRAKER/Ref_RnaSeq/02_PublishedData/; module load SAMtools/1.13-GCC-10.3.0; samtools view --threads 48 -b -o PRJNA863910_30samples_Hap1.bam PRJNA863910_30samples_Hap1.sam; samtools sort -o PRJNA863910_30samples_Hap1_sorted.bam -T PRJNA863910_30samples_Hap1_temp --threads 48 PRJNA863910_30samples_Hap1.bam"
