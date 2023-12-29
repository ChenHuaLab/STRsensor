STRsensor Usage Guide
======================

**Note**
* The version (e.g. hg19 ) of reference genome must be exactly the same as the coordinate defined in the STR locus file. <br>
* The chromosome (e.g. chr12) in reference genome must be exactly the same as defined in the STR locus file. <br>


**Example**

    Example/
    ├── Data
    │   ├── Fastq     // simulated single-end fastq files
    │   │   ├── P1_R1.fq.gz
    │   │   ├── P2_R1.fq.gz
    │   │   ├── P3_R1.fq.gz
    │   │   └── ...
    │   └── Mapping   // align the fastq read to GRCH37 reference genome
    │       ├── P1.bam
    │       ├── P1.bam.bai
    │       └── ...
    ├── Materials
    │   ├── BamList.txt     // the BAM file path of each sample
    │   ├── InFreq.txt      // input frequency parameters of each STR locus
    │   ├── InStutter.txt   // input stutter parameters of each STR locus
    │   └── STRLocus.txt    // STR locus information (GRCH37 coordinate)
    ├── Prog
    │   └── Mapping.sh      // bash script used to align the fastq to reference
    ├── README.md
    └── Result        // output directory contains allele information of each locus
        ├── CSF1PO.txt
        ├── D12S391.txt
        ├── D13S317.txt
        └── ...

**Processing**

1. Dependence

       (1) BWA (v0.7.17)
       (2) samtools (v1.10)

2. Preparation

        # download reference genome
        ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

        # build reference index
        samtools faidx GRCh37.p13.genome.fa  --> GRCh37.p13.genome.fa.fai

        # build mapping index
        bwa index GRCh37.p13.genome.fa

3. Mapping

        # mapping (BWA is recommended)
        bwa mem -M <GRCH37.p13.genome.fa> <In.fq.gz> | samtools view -bS | samtools sort -o <Out.bam>

4. Running STRsensor

        # Go to the Example folder
        cd STRsensor/Example

        # Run the following command
        ../STRsensor-1.X
            -i Materials/BamList.txt 
            -r Materials/STRLocus.txt 
            -f Downloads/Reference/GRCh37.p13.genome.fa 
            -o Result 
            -s Materials/InStutter.txt 
            -q Materials/InFreq.txt

