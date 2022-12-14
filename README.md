**Stage 3 Project**

**#Project- Genotypying Malaria Drug Resistance
**

INTRODUCTION: The global burden of malaria has increased during the last decade and it continues to cause signifcant morbidity and mortality globally, despite the control interventions. Current malaria control and elimination strategies rely mainly on effective antimalarial drugs. However, drug resistance has become a major threat facing malaria control programs because Plasmodium falciparum as a parasite has adapted to anti-malarials and this can be attributed partly to its enormous genetic diveristy. Identifying, monitoring and understanding the molecular markers of adaptive mechanisms against anti-malarials are of importance as anti-malarial therapy and prophylaxis have been key in reducing mortality and morbidity from malaria. Although, mutations in genes determining resistance to drugs such as chloroquine have been identified, there is still a need to have full understanding of the resistance mechanisms, and genes that contribute to resistance to many other drugs.

AIM: This project aims to identify the presence of one or more of these variants (Pfcrt, Pfmdr1, Pfdhfr, and Pfdhps, Pfarps10, Pfferredoxin, Pfexonuclease and Pfmdr2) in Plasmodium falciparum genome sequences from different countries in the world and make drug prescriptions that fit each.

METHODOLOGY: Five (5) Whole Genome Sequences each of Plasmodium falciparum from three countries – China-Myanmar border, Malawi, and Vietnam were retrieved from the SRA acrchive of NCBI (https://www.ncbi.nlm.nih.gov/sra) a public biological database. The Accession IDs were used to retrieve the datasets and their metadata from the SRA (https://sra-explorer.info/), an open and free access repository of high throughput sequencing data. The link to the metadata of the datasets is https://github.com/ruthobaado/Project-Malaria-Drug-Resistance/blob/main/Metadata_samples.tsv The script file for the codes used for this project can be accessed via the link: https://github.com/ruthobaado/Project-Malaria-Drug-Resistance/blob/main/Stage%203.sh

Software packages used for the analysis: FastP: This was used to trim the adapters in the datasets. SPAdes 3.15.5: This was used for the genome assembly of the trimmed data. #Bash scripts were used for trimming and the genome assembly ResFinder 4.1: The contigs.fasta files of the assembled genomes were uploaded to ResFiinder (https://cge.food.dtu.dk/services/ResFinder-4.1/), an open online resource for identification of antimicrobial resistance genes in high-throughput sequencing data and prediction of phenotypes from genotypes. Bioinformatics webtools was used to visualize the compiled result from ResFinder (https://bioinformatics.psb.ugent.be/webtools/Venn/)

The image and the summary table for the visualized result can be accessed via the links: 
https://github.com/ruthobaado/Project-Malaria-Drug-Resistance/blob/main/venn_result%20(3%20Countries).png
https://github.com/ruthobaado/Project-Malaria-Drug-Resistance/blob/main/Summary-table.txt
