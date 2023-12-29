STRsensor
=========================
STR allele-typing with both WGS dataset and multiplex deep sequencing dataset.


__PROGRAM: STRsensor__<br>
__VERSION: 1.2.1__<br>
__PLATFORM: Linux__<br>
__COMPILER: gcc-4.8.5__<br>
__AUTHOR: xiaolong zhang__<br>
__EMAIL: xiaolongzhang2015@163.com__<br>
__DATE:   2020-10-28__<br>
__UPDATE: 2023-12-29__<br>
__DEPENDENCE__<br>
* __GNU make and gcc__<br>



Description
=========================
* STRsensor is designed for typing of STR loci, which including forensic CODIS STR locus and user-specified STR locus. <br>
* It is developed to work with both low-coverage WGS data and high-coverage target/amplicon sequencing data.
* 'Cigar-based' and 'Kmer-based' matching algorithm are dperformed to extract enough alleles from spanning reads for final-allele determination.<br>
* The Maximum Likelihood Estimation (MLE) is performed to quantify the effect of PCR stutter induced by polymerase slippages on STR typing.<br>
* The Maximum A Posteriori estimation (MAP) is applied to obtain the most likely STR allele from the candidates.<br>



Building
=========================

See INSTALL for complete details.


Usage
========================

      Options:
            -h|--help             print help infomation

      [Required]
            -i|--infile     FILE   bam file list, one sample per line [.txt]
            -r|--region     FILE   region file of the STR locus [.txt]
            -f|--fasta      FILE   reference genome (must be same as the mapping process used) [.fa]
            -o|--outpath    PATH   the path for output file

      [Optional input parameters]
            -s|--stutter    FILE   the stutter parameters file for each STR locus [.txt]
            -q|--frequency  FILE   the frequency file for each STR locus on each potential allele [.txt]

      [Optional output parameters]
            -p|--params_out        output parameters learned by given samples (including
                                   stutter and allele frequency)

      [Optional filter rules]
            -c|--min_prob   FLOAT  the minimal probability allowed to call an allele [0.00]
            -n|--min_reads  INT    the minimum reads needed to genotype a STR locus for an individual [10]
            -m|--mis_match  INT    the maximum mismatch(mismatch/insert/delete) bases allowed in both
                                   3' and 5' flanking sequence [2]

      [Optional]
            -a|--allow_dup         duplication is allowed in target and amplicon sequencing, but not in WGS
            -t|--threads    INT    the number of threads used [1]

                          

Option
========================

#### \[-h|--help]
      Print the help infomation

#### \[-i|--infile]
      A text file (without header) contains the full path of each samples, one sample per line. 
      Each BAM file (*.bam) should have a corresponding index file (*.bam.bai)
      e.g. 
            /home/xlzh/data/sample1.bam
            /home/xlzh/data/sample2.bam
            ...

#### \[-r|--region]
      A text file (with header) contains the details of each STR locus.
      e.g.
            STRLocus  Chrom  Start     End       MotifLen  NotCountedBase  AsHaplotype
            D19S433   chr19  30417142  30417205  4         8               No

      [FIELDS]
            1. STRLocus: the locus name of the STR
            2. Chrom: the chromosome of the STR locus
            3. Start: the start position in the reference genome (1-based)
            4. End: the end position in the reference genome (1-based)
            5. MotifLen: the motif length of the locus.
            6. NotCountedBases: the number of bases should be excluded from the STR region
            7. AsHaplotype: regard the STR locus as haplotype or not ('Yes' | 'No')

      [NOTE]
            For STR locus with two fragments, such as DYS385ab, the two fragments need to be 
            defined in two lines and ended with the identifier of lowercase character of "a" 
            and "b" (DYS385a and DYS385b).
            e.g.
                DYS385a  chrY  20801599  20801642  4  0  No
                DYS385b  chrY  20842518  20842573  4  0  No

#### \[-f|--fasta]
      Reference genome (FASTA), which should be the one used in sequence alignment.
      e.g. GRCH37_hg19.fa

#### \[-o|--outpath]
      The OUT_PATH for output files.
      eg. /home/xlzh/result

#### \[-s|--stutter]
      A text file (without header) contains the stutter parameters of each STR locus.
      e.g.
            DYS447  0.964309  0.007595  0.001781  0.002708  0.023607
      
      [FIELDS]
            1. STRLocus:  the locus name of the STR
            2. p_normal:  probability of NO STUTTER occured
            3. p_add1:    probability of one repeat unit adds
            4. p_add2:    probability of two repeat unit adds
            5. p_remove2: probability of two repeat unit removes
            6. p_remove1: probability of one repeat unit removes

#### \[-q|--frequency]
      A text file (without header) contains the frequency distribution of each STR locus
      e.g.
            DYS437  14.0:0.645401  15.0:0.339763  16.0:0.005935  13.0:0.007418  17.0:0.001484

      [FIELDS]
            1. STRLocus: the locus name of the STR
            2. allele1:frequency1
            3. allele2:frequency2
            ...

#### \[-p|--params_out]
      Output the parameters (stutter and frequeny) learned from the given samples.

      [Output]
          File1: OUT_PATH/Stutter.txt
          File2: OUT_PATH/Frequency.txt

      [Rules]
          1. parameters of the locus is not defined in stutter/frequency file
             (set by '--stutter' and '--frequency' options)
          2. the number of extracted alleles is greater than 100 for stutter evaluation
          3. the number of given samples is greater than 20 for frequency evaluation
      
#### \[-c|--min_prob]
      The minimal probability allowed to call an allele.
      Default: 0.0

#### \[-n|--min_reads]
      The minimum number of spanning reads needed to genotype a STR locus for an individual.
      Default: 10

#### \[-m|--mis_match]
      The maximum number of mismatch (mismatch/insert/delete) bases allow in both 3' 
      and 5' flanking sequence.
      Default: 2
      
#### \[-a|--allow_dup]
      Duplication is only recommended for amplicon/target sequencing dataset.
      For WGS dataset, NOT recommend to set this option.

#### \[-t|--threads]
      The number of threads used by STRsensor at run time.
      Default: 1


STR Locus
=========================
**One Fake STR Locus**

            5' flanks                       STR Region                          3' flanks
        ++++++++++++++++++TCTATCTATCTATCTATCTA AACC TCTATCTATCTATCTATCTATCTA++++++++++++++++++
                          |                    \  /                        |
                      STR_Start          Not_Counted_Bases              STR_End
                      (149347)                                          (149394)


**Locus.txt**

    STRLocus<TAB>Chrom<TAB>Start<TAB>End<TAB>MotifLen<TAB>NotCountedBase<TAB>AsHaplotype
    FakeSTR<TAB>chr2<TAB>149347<TAB>149394<TAB>4<TAB>4<TAB>No

**Description**

* The name of STRLocus should **NOT EXCEED 64** characters!
* The 'STR_Start' is the start position (1-based) of STR region in the reference genome
* The 'STR_End' is the end position (1-based) of STR region in the reference genome
* The motif is 'TCTA', therefore, MotifLen = 4
* The sequence of 'AACC' should be excluded from allele determination, therefore, NotCountedBases = 4
* The FakeSTR belong to chromosome 2, therefore, AsHaplotype = 'No'



Output
=========================
All the files generated by STRsensor will be output to **OUT_PATH** folder<br>
e.g.<br>

    --outpath /home/xlzh/result

(1) Stutter.txt & Frequency.txt (if '--params_out' is given)
* The **FILE FORMAT** is exactly the same as described in the options of '--stutter' and '--frequency'. <br><br>

(2) LOCUS_NAME.txt <br> 
* Example: DYS19.txt

      #Locus  DYS19
      #Position       chrY:9521989-9522052
      #MotifLength    4
      #IsHaplotype    Yes
      #MinimalReads   30
      #MinimalProbability     0.99
      #CandidateAllele        13.0:0.033794,14.0:0.205837,15.0:0.491551,16.0:0.205837,17.0:0.061444,18.0:0.001536
      #StutterParams  d_2:0.004609,d_1:0.068775,n:0.920281,u_1:0.005676,u_2:0.000660
      #Sample Allele  Probability     TotalReads      ValidReads      SpanReads       AlleleList
      Sample_25.bam   15.0    1.000000        2548    2506    2290    15.0:2108,14.0:151,16.0:21,13.0:6,17.0:3,10.0:1
      Sample_65.bam   15.0    1.000000        3274    3253    3154    15.0:2912,14.0:212,16.0:18,13.0:12
      Sample_82.bam   17.0    1.000000        2776    2766    2679    17.0:2424,16.0:210,18.0:26,15.0:16,14.0:2,13.0:1

* Description and Fields

      [Description]
            Allele information for each individual at the locus, one sample per line

      [StutterParams]
            d_2: probability of two repeat unit removes
            d_1: probability of one repeat unit removes
              n: probability of NO STUTTER occured
            u_1: probability of one repeat unit adds
            u_2: probability of two repeat unit adds

      [Fields]
            1. Sample: sample name that extracted from sample's path
            2. Allele: allele that determined by STRsensor
            3. Probability: the probability of the allele [0.0 ~ 1.0]
            4. TotalReads: the total number of reads that have overlap with STR region
            5. ValidReads: the number of reads that passed the filter rules
            6. SpanReads: the number of reads that fully span the entire STR region
            7. AlleleList: extracted allele and its corresponding number 

