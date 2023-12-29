/*************************************************************************
    > File Name: cigar.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年11月19日 星期四 09时28分52秒
    > Updated Time: 2021年01月03日 星期六 20时51分18秒
 ************************************************************************/

#ifndef CIGAR_H
#define CIGAR_H

#include "bamio.h"

#ifndef MIN_MATCH
  #define MIN_MATCH 8
#endif

/*! @typedef locate_t
 @abstract structure for flank location on the read sequence
 @field status       flank query status (0->failed; 1->succeed)
 @field loc5         located start of flanks5' [int64_t]
 @field loc3         located end of flanks3' [int64_t]
*/
typedef struct {
    int status;
    int64_t loc5;
    int64_t loc3;
} locate_t;


/*! @typedef indel1_t
 @abstract structure for each processed cigar 
 @field ref_pos       reference position (1-based)
 @field c_clen        the relative postion of cigar in the read sequence
 @field c_value       the cigar length. eg. 105M -> 105
 @field c_type        the cigar type. eg. 105M -> 'M'
*/
typedef struct {
    int64_t ref_pos;
    int c_clen;
    int c_value;
    char c_type;
} indel1_t;


/*! @typedef indel_t
 @abstract structure for processed cigar list
 @field n            number of cigars [int]
 @field m            max number of memory allocated [int]
 @field core         pointer to the new processed cigar [indel1_t *]
*/
typedef struct {
    int n, m;
    indel1_t *core;
} indel_t;


/*! @function: get allele by cigar value provided by alignment
 * @param read_seq          the read sequence used to check the flank sequencing
 * @param r_cigar           the raw cigar list parsed from alignment
 * @param ref_pos           the reference position (1-based)
 * @param locus             STR locus object [locus_t *]
 * @return                  detected allele [float]
*/
float cigar_allele(char *read_seq, xcigar_t *r_cigar, int64_t ref_pos, locus_t *locus, int mis_match);


#endif /* CIGAR_H*/