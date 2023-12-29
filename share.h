/*************************************************************************
    > File Name: share.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2020年10月27日 星期一 16时24分32秒
 ************************************************************************/

#ifndef SHARE_H
#define SHARE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "index.h"

#define __STRSENSOR_VERSION__ "1.2.1"
#define __STRSENSOR_CREATE_DATE__ "2020-10-28"
#define __STRSENSOR_UPDATE_DATE__ "2023-12-29"


/*! @typedef arg_t
 @abstract structure for the comand line args
 @field help            [0|1] 1: print the help infomation
 @field allow_dup       [0|1] 1: keep the read that is duplicate
 @field params_out      [0|1] 1: output the parameters (stutter and frequency)
 @field threads         the number of threads used
 @field kmer            kmer=12(Amplicon; allow_dup = 1), kmer=8(WGS; allow_dup = 0)
 @field mismatch        maximum mismatch bases (mismatch/insert/delete) allowed
 @field min_reads       configuration file with level information
 @field cutoff          minimum probability to call an allele
 @field infile          bam file list (eg. samples.list)
 @field outpath         output file path (eg. /home/xlzh/Result/)
 @field region          region file of the STR Locus
 @field fasta           reference genome (eg. hg19.fa)
 @field stutter         stutter parameter for each STR locus
 @field frequency       frequency file for each locus on each potential allele
*/
typedef struct arg_t {
    int help;
    int allow_dup, params_out;
    int threads, kmer;
    int mis_match, min_reads;
    float min_prob;
    char *infile, *outpath;
    char *region, *fasta;
    char *stutter, *frequency;
} arg_t;


/*! @typedef locus_t
 @abstract structure for the STR locus
 @field name              STR locus name, eg. DYS390
 @field chr               chromosome of the STR locus
 @field start, end        start and end at the reference genome
 @field motif_len         length of STR motif
 @field exclude_base      the number of bases should excluded from the core STR repeat region
 @field index_kmer        minimum kmer length for indexing
 @field is_haplotype      the status of STR locus (DYS389I->is haplotype; D2S1338->not haplotype)
 @field fa_obj;           the pointer to the fasta object [fasta_t *]
 @field index5, index3    index structure for flank5 and flank3 (index->{kmer, flank, idx})
*/
typedef struct locus_t {
    char name[0x40];
    char chr[0x40];
    uint32_t start, end;
    uint32_t motif_len;
    uint32_t exclude_base;
    uint32_t index_kmer;
    uint32_t is_haplotype;
    fasta_t *fa_obj;
    index_t *index5, *index3;
} locus_t;


/*! @typedef region_t
 @abstract structure for the STR region
 @field n                 number of STR locus [int]
 @field m                 recorde the memory allocated locus_t [int]
 @field loci              locus name for each STR region [char **]
 @field locus             structure for each locus [locus_t]
*/
typedef struct region_t {
    int n, m;
    char **loci;
    locus_t *locus;
} region_t;


/* ------------------------ prototype function ------------------------ */

/*! @function
  * @abstract     command line parameters parsing.
  @param  argc    the number of args [int].
  @param  argv    the parameters [char **].
  @return         args structure [arg_t *].
 */
arg_t *args_parse(int argc, char **argv);


#endif // SHARE_H
