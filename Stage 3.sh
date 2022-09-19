#make a new directory for project work
mkdir Oluwaseun

mkdir WGS_datasets
cd WGS_datasets

#download datasets (from SRA Archive):
#China-Myanmar border (Asia)
mkdir China-Myanmar
cd China-Myanmar
nano China-Myanmar.sh
#enter the download links and save 
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/002/ERR3325132/ERR3325132_1.fastq.gz -o ERR3325132_Illumina_HiSeq_2000_paired_end_sequencing_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/002/ERR3325132/ERR3325132_2.fastq.gz -o ERR3325132_Illumina_HiSeq_2000_paired_end_sequencing_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/001/ERR3325131/ERR3325131_1.fastq.gz -o ERR3325131_Illumina_HiSeq_2000_paired_end_sequencing_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/001/ERR3325131/ERR3325131_2.fastq.gz -o ERR3325131_Illumina_HiSeq_2000_paired_end_sequencing_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/003/ERR3325133/ERR3325133_1.fastq.gz -o ERR3325133_Illumina_HiSeq_2000_paired_end_sequencing_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/003/ERR3325133/ERR3325133_2.fastq.gz -o ERR3325133_Illumina_HiSeq_2000_paired_end_sequencing_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/004/ERR3325134/ERR3325134_1.fastq.gz -o ERR3325134_Illumina_HiSeq_2000_paired_end_sequencing_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/004/ERR3325134/ERR3325134_2.fastq.gz -o ERR3325134_Illumina_HiSeq_2000_paired_end_sequencing_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/005/ERR3325135/ERR3325135_1.fastq.gz -o ERR3325135_Illumina_HiSeq_2000_paired_end_sequencing_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR332/005/ERR3325135/ERR3325135_2.fastq.gz -o ERR3325135_Illumina_HiSeq_2000_paired_end_sequencing_2.fastq.gz

bash China-Myanmar.sh


#Malawi (Africa):
mkdir Malawi
cd Malawi
nano Malawi.sh
#enter the download links and save 
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/002/SRR7465152/SRR7465152_1.fastq.gz -o SRR7465152_Other_Sequencing_of_malaria_parasite_P._falciparum_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/002/SRR7465152/SRR7465152_2.fastq.gz -o SRR7465152_Other_Sequencing_of_malaria_parasite_P._falciparum_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/004/SRR7464834/SRR7464834_1.fastq.gz -o SRR7464834_Other_Sequencing_of_malaria_parasite_P._falciparum_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/004/SRR7464834/SRR7464834_2.fastq.gz -o SRR7464834_Other_Sequencing_of_malaria_parasite_P._falciparum_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/001/SRR7465161/SRR7465161_1.fastq.gz -o SRR7465161_Other_Sequencing_of_malaria_parasite_P._falciparum_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/001/SRR7465161/SRR7465161_2.fastq.gz -o SRR7465161_Other_Sequencing_of_malaria_parasite_P._falciparum_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/000/SRR7465160/SRR7465160_1.fastq.gz -o SRR7465160_Other_Sequencing_of_malaria_parasite_P._falciparum_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/000/SRR7465160/SRR7465160_2.fastq.gz -o SRR7465160_Other_Sequencing_of_malaria_parasite_P._falciparum_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/002/SRR7465162/SRR7465162_1.fastq.gz -o SRR7465162_Other_Sequencing_of_malaria_parasite_P._falciparum_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR746/002/SRR7465162/SRR7465162_2.fastq.gz -o SRR7465162_Other_Sequencing_of_malaria_parasite_P._falciparum_2.fastq.gz

bash Malawi.sh


#Vietnam (Asia):
mkdir Vietnam
cd Vietnam
nano Vietnam.sh
#enter the download links and save 
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629073/SRR629073_1.fastq.gz -o SRR629073_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629073/SRR629073_2.fastq.gz -o SRR629073_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530502/SRR530502_1.fastq.gz -o SRR530502_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-69845_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR530/SRR530502/SRR530502_2.fastq.gz -o SRR530502_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-69845_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629066/SRR629066_1.fastq.gz -o SRR629066_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629066/SRR629066_2.fastq.gz -o SRR629066_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629042/SRR629042_1.fastq.gz -o SRR629042_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629042/SRR629042_2.fastq.gz -o SRR629042_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629067/SRR629067_1.fastq.gz -o SRR629067_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR629/SRR629067/SRR629067_2.fastq.gz -o SRR629067_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO__2.fastq.gz

bash Vietnam.sh


#Trimming with FastP:
cd China-Myanmar
nano trim.sh

#enter the command codes below and paste in the trim.sh file
#!/bin/bash
mkdir qc_reads
SAMPLES=(
  "ERR3325131_Illumina_HiSeq_2000_paired_end_sequencing”
  “ERR3325132_Illumina_HiSeq_2000_paired_end_sequencing
  “ERR3325133_Illumina_HiSeq_2000_paired_end_sequencing"
  “ERR3325134_Illumina_HiSeq_2000_paired_end_sequencing"
  “ERR3325135_Illumina_HiSeq_2000_paired_end_sequencing”
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_1.fastq.gz" \
    -I "$PWD/${SAMPLE}_2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done

bash trim.sh
mv qc_reads trimmed_reads


cd Malawi
nano trim.sh

#enter the command codes below and paste in the trim.sh file
#!/bin/bash
mkdir qc_reads
SAMPLES=(
  "SRR7464834_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465152_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465160_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465161_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465162_Other_Sequencing_of_malaria_parasite_P._falciparum"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_1.fastq.gz" \
    -I "$PWD/${SAMPLE}_2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done

bash trim.sh
mv qc_reads trimmed_reads


cd Vietnam
nano trim.sh

#enter the command codes below and paste in the trim.sh file
#!/bin/bash
mkdir qc_reads
SAMPLES=(
  "SRR530502_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-69845_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629042_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629066_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629067_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629073_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
)

for SAMPLE in "${SAMPLES[@]}"; do

  fastp \
    -i "$PWD/${SAMPLE}_1.fastq.gz" \
    -I "$PWD/${SAMPLE}_2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
done

bash trim.sh
mv qc_reads trimmed_reads



#Genome Assembly with Spades.py:
cd China-Myanmar
cd trimmed_reads
nano CM.sh
#enter the command codes below and paste in the .sh file
#!/bin/bash
samples=(
  "ERR3325131_Illumina_HiSeq_2000_paired_end_sequencing"
  "ERR3325132_Illumina_HiSeq_2000_paired_end_sequencing"
  "ERR3325133_Illumina_HiSeq_2000_paired_end_sequencing"
  "ERR3325134_Illumina_HiSeq_2000_paired_end_sequencing"
  "ERR3325135_Illumina_HiSeq_2000_paired_end_sequencing"
)
 
for sample in ${samples[@]}
 do
    spades.py -1 "$PWD/$sample"_1.fastq.gz -2 "$PWD/$sample"_2.fastq.gz -o "CMassembly"
 done

bash CM.sh


cd Malawi
cd trimmed_reads
nano Mal.sh
#enter the command codes below and paste in the .sh file
#!/bin/bash
samples=(
  "SRR7464834_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465152_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465160_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465161_Other_Sequencing_of_malaria_parasite_P._falciparum"
  "SRR7465162_Other_Sequencing_of_malaria_parasite_P._falciparum"
)
 
for sample in ${samples[@]}
 do
    spades.py -1 "$PWD/$sample"_1.fastq.gz -2 "$PWD/$sample"_2.fastq.gz -o "Malassembly"
 done

bash Mal.sh


cd Vietnam
cd trimmed_reads
nano Viet.sh
#enter the command codes below and paste in the .sh file
#!/bin/bash
samples=(
  "SRR530502_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-69845_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629042_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629066_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Solexa-120570_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629067_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
  "SRR629073_Illumina_whole_genome_shotgun_sequencing_of_genomic_DNA_paired-end_library_Pond-108958_containing_sample_Plasmodium_falciparum_Vietnam_Oak-Knoll_FVO_"
)
 
for sample in ${samples[@]}
 do
    spades.py -1 "$PWD/$sample"_1.fastq.gz -2 "$PWD/$sample"_2.fastq.gz -o "Vietassembly"
 done

bash Viet.sh
