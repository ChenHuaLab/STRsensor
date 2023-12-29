/*************************************************************************
    > File Name: bamio.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2020年11月04日 星期三 22时06分30秒
 ************************************************************************/

#ifndef BAMIO_H
#define BAMIO_H


#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "htslib/sam.h"
#include "htslib/kstring.h"

/*! @function 
   @abstract Get the next read from a BAM/CRAM multi-iterator
   @param htsfp       Htsfile pointer for the input file
   @param itr         Iterator
   @param r           Pointer to a bam1_t struct
   @return >= 0 on success; -1 when there is no more data; < -1 on error
*/
#define xbam_next(htsfp, itr, aln) sam_itr_next(htsfp, itr, aln)

/*! @function
 @abstract  Get whether the query is secondary
 @param  b  pointer to an alignment
 @return    boolean true if query is secondary
 */
#define xbam_is_secondary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)

/*! @function
 @abstract  Get whether the query is supplementary
 @param  b  pointer to an alignment
 @return    boolean true if query is supplementary
 */
#define xbam_is_supplementary(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)


/* @fucntion
 Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T 
 and 15 for N. Two bases are packed in one byte with the base
 at the higher 4 bits having smaller coordinate on the read.
*/
static const char HTSLIB_BASE_TABLE[0x10] = {
    ' ', 'A', 'C', ' ',  //  0,  (1), (2), 3
    'G', ' ', ' ', ' ',  // (4),  5 ,  6,  7
    'T', ' ', ' ', ' ',  // (8),  9,   10, 11
    ' ', ' ', ' ', 'N'   // 12,   13,  14, (15)
};


/*! @typedef xbam_t
 @abstract structure of bam info structure
 @field bam_hd           bam file handle [htsFile *]
 @field header           bam file header [sam_hdr_t *]
 @field bam_idx          bam file index [hts_idx_t *]
 @field aln              alignment [bam1_t *]
 @field iter             iterator for bam file (DYS385ab: iter[0] and iter[1])[hts_itr_t *]
*/
typedef struct xbam_t {
    htsFile *bam_hd;
    bam_hdr_t *header;
    hts_idx_t *bam_idx;
    bam1_t *aln;
    hts_itr_t *iter[2];
    int n_iter;
} xbam_t;


/*! @typedef xcigar_t
 @abstract structure of bam cigar structure
 @field n                number of cigars [int]
 @field m                max number of memory allocated [int]
 @field cigar_op         cigar operator [char *]
 @field bam_idx          cigar length [int *]
*/
typedef struct xcigar_t {
    int n, m;
    char *cigar_op;
    int *cigar_len;
} xcigar_t;


/* prototype function */
void xbam_seq(bam1_t *aln, kstring_t *ks);
void xbam_cigar(bam1_t *aln, xcigar_t *cigar);
xcigar_t *xbam_cigar_init(int n_cigar);
xbam_t *xbam_init(char *bam_file);
int xbam_fetch(char *chr, int start, int end, xbam_t *bam);
void xbam_destroy(xbam_t *bam);
void xbam_cigar_destroy(xcigar_t *cigar);


#endif