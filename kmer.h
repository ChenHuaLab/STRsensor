/*************************************************************************
    > File Name: kmer.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月30日 星期五 12时40分53秒
 ************************************************************************/

#ifndef KMER_H
#define KMER_H

#include "share.h"
#include "utils.h"


/*! @typedef loc_t
 @abstract structure for kmer location on the read sequence
 @field status       flank query status (0->failed; 1->succeed)
 @field loc5         located start and end of 5' flanks [uint32_t]
 @field loc3         located start and end of 3' flanks [uint32_t]
*/
typedef struct {
    int status;
    uint32_t loc5[2];
    uint32_t loc3[2];
} loc_t;


/*! @typedef hit_t
 @abstract structure for index of kmer hit on the read sequence
 @field k_idx        index on the kmer sequence (0-based) [uint32_t]
 @field r_idx        index on the read sequence (0-based) [uint32_t]
*/
typedef struct {
    uint32_t k_idx;
    uint32_t r_idx;
} hit_t;


/*! @typedef list_t
 @abstract structure for index of kmer hit on the read sequence
 @field n        number of elements [size_t]
 @field m        maximum allowed number of elements [size_t]
 @field a        hit index for each kmer [loc_t *]
*/
typedef struct {
    size_t n, m;
    hit_t *a;
} list_t;


/* core sort function for hits of slice_list */

static int _compare_ascend(const void *a, const void *b) {
    /* sorted by ascending order */
    return ((hit_t*)a)->k_idx - ((hit_t*)b)->k_idx;
}

static int _compare_descend(const void *a, const void *b) {
    /* sorted by descending order */
    return ((hit_t*)b)->k_idx - ((hit_t*)a)->k_idx;
}


/* list operation function */
#define x_destroy(v) free((v)->a)
#define x_size(v) ((v)->n)
#define x_get(v) ((v)->a)

#define x_push(v, x) do {                                               \
		if ((v)->n == (v)->m) {                                         \
			(v)->m = (v)->m? (v)->m<<1 : 16;                            \
			(v)->a = (hit_t*)realloc((v)->a, sizeof(hit_t) * (v)->m);   \
		}                                                               \
		(v)->a[(v)->n++] = (x);                                         \
	} while (0)


/*! @function: get allele size from read sequence by 'k-mers' algorithm.
  @param  read_seq    Name of the hash table [char *].
  @param  locus       Pointer to the hash table [locus_t *].
  @param  mis_match   Key [int].
  @return             size of allele deteted [float]
 */
float kmer_allele(char *read_seq, locus_t *locus, int mis_match);

#endif // KMER_H