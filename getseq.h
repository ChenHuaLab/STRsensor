/*************************************************************************
    > File Name: getseq.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月28日 星期三 15时53分18秒
 ************************************************************************/

#ifndef GETSEQ_H
#define GETSEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>


#ifdef _WIN32
  #define ftell ftello64
  #define fseek fseeko64
#endif


typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;

typedef struct {
    FILE *fai_fp;
    int32_t n;
    char **name;
    faidx1_t *line_fai;
} faidx_t;

typedef struct {
    char *fa_file;
    FILE *fa_fp;
    faidx_t *fai;
} fasta_t;


/*! @funciton: fasta (.fa) initiate
 *   @parmeters:
 *   @    argc        the file name of reference genome [char *]
 *   @return:
 *   @    fasta_obj   the structure of fasta [fasta_t *]
*/
fasta_t *fa_init(char *fa_file);


/*! @funciton: get sequene from reference genome by given position
 *   @parmeters:
 *   @    chr            the chromosome [char *]
 *   @    start, end     start and end of the STR locus [uint32_t]
 *   @    fa_obj         the fasta struct [fasta_t *]
 *   @return:
 *   @    sequence       sequence from reference genome [char *]
*/
char *get_seq(char *chr, int32_t start, int32_t end, fasta_t *fa_obj);


#endif // GETSEQ_H