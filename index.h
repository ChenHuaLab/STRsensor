/*************************************************************************
    > File Name: index.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月29日 星期四 10时15分20秒
 ************************************************************************/

#ifndef INDEX_H
#define INDEX_H

#include "getseq.h"
#include "utils.h"
#include "htslib/khash.h"


KHASH_MAP_INIT_STR(index, uint32_t) /* hash table: str -> uint32_t */

#define ShiftScope 200 /* defined search length of flanks during indexing */
#define MaxFlank 25 /* defined flanks length*/

/* type of flanking sequence */
enum FLANKTYPE { FLANK5 = 0, FLANK3 = 1 };


/*! @typedef index_t
 @abstract  structure for the STR flank index
 @field kmer        actual kmer length used in flanks indexing
 @field flank       the flanking sequence (length is defined as "ShiftScope")
 @field idx         hash table that contain the index of flanks
*/
typedef struct {
    uint32_t kmer;
    char *flank;
    khash_t(index) *idx;
} index_t;


/*! @function
  * @abstract           get allele size from read sequence by 'k-mers' algorithm.
  @param  chr           the chromosome [char *].
  @param  start         start of the STR locus [uint32_t].
  @param  end           end of the STR locus [uint32_t].
  @param  fa_obj        the fasta struct [fasta_t *].
  @param  motif_len     the motif length [uint32_t].
  @param  index_kmer    the minimum length of kmer indexing [uint32_t].
  @param  flank_type    flanks type (FLANK5 and FLANK3) [int].
  @return               index structure [index_t *].
 */
index_t *flanks_index(char *chr, uint32_t start, uint32_t end, fasta_t *fa_obj, 
            uint32_t motif_len, uint32_t index_kmer, int flank_type);


#endif // INDEX_H