/*************************************************************************
    > File Name: kmer.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月30日 星期五 12时41分22秒
 ************************************************************************/

#include "kmer.h"
#include "index.h"


/*! @function: slice the read sequence to get all potential kmer hits
 * @param read_seq          the read sequence of the alignment [char *].
 * @param flank_idx         flanking sequence index [khash_t(index) *].
 * @param kmer              the kmer length used in the flanks indexing [uint32_t].
 * @return                  the hit list for the read sequence kmer [list_t *].
*/
static list_t *get_hit(char *read_seq, khash_t(index) *flank_idx, uint32_t kmer)
{
    list_t *slice_list;
    hit_t hit;

    err_calloc(slice_list, 1, list_t);
    int shift_len = strlen(read_seq) - kmer;

    khint_t k;
    char key[0xff];

    for (int i=0; i <= shift_len; ++i) {
        strncpy(key, read_seq+i, kmer); key[kmer] = '\0';

        k = kh_get(index, flank_idx, key);
        if (k == kh_end(flank_idx)) continue; /* the key is not existed in the hash-table */
        hit.k_idx = kh_value(flank_idx, k); hit.r_idx = i;
        x_push(slice_list, hit);
    }

    return slice_list;
}


/*! @function: check the hamming distance between the two given sequence.
 * @param str1          the first sequence used to check the hamming distance.
 * @param str2          the second sequence used to check the hamming distance.
 * @param min_match     minimum matched bases need [int].
 * @param mis_match     number of mismath allowed [int].
 * @return              pass the check (1) and failed to check (0)
*/
static int hamming_check(char *str1, char *str2, int min_match, int mis_match)
{
    int len1, len2, min_len;
    int match_num=0, mis_num=0;

    #define ACCEPT 1
    #define DISCARD 0

    len1 = strlen(str1); len2 = strlen(str2);
    min_len = len1 < len2 ? len1 : len2;

    for (int i=0; i < min_len; ++i) {
        if (str1[i] != str2[i]) mis_num += 1;
        else match_num += 1;
        if (mis_num > mis_match) return DISCARD;
    }
    if (match_num < min_match) return DISCARD;

    return ACCEPT;
}


/*! @function: get the loc index of flank_seq on the read_seq.
  @param  read_seq      the read sequence [char *].
  @param  flank_seq     the flanking sequence [char *].
  @param  loc_list      hit list of read sequence kmer [list_t *].
  @param  min_match     minimum matched bases need [int].
  @param  mis_match     number of mismath allowed [int].
  @return               location index of the flanks on the read sequence [loc_t]
 */
static loc_t flank_check(char *read_seq, char *flank_seq, list_t *loc_list, int min_match, int mis_match)
{
    loc_t loc = {0};
    
    hit_t *hit_list = x_get(loc_list); /* list for each kmer hit */
    int read_len = strlen(read_seq);
    int flank_len = strlen(flank_seq);

    int shift, pstart, pend, accept;
    char *r_seq, *f_seq; /* part of read_seq or flank_seq used */

    for (int i=0; i < x_size(loc_list); ++i) {
        shift = hit_list[i].r_idx - hit_list[i].k_idx;

        if (shift >= 0) {
            r_seq = read_seq + shift;
            f_seq = flank_seq;
        }
        else { /* shift < 0 */
            r_seq = read_seq;
            f_seq = flank_seq + (-1 * shift);
        }

        pstart = shift < 0 ? 0 : shift;
        accept = hamming_check(r_seq, f_seq, min_match, mis_match);
        if (accept != 1) return loc; /* it is not a proper kmer match */

        if (shift < 0) pend = flank_len + shift - 1;
        else if (pstart+flank_len >= read_len) pend = read_len - 1;
        else pend = pstart + flank_len - 1;

        loc.status = 1; 
        loc.loc5[0] = pstart; loc.loc5[1] = pend;
        return loc;
    } 
    /* failed to get proper kmer match finally */
    loc.status = 0; return loc;
}


/*! @function: get the location of both 5' and 3' flanks.
    loc_idx.status: 
        0 -> can't find any kmer match
        1 -> only flank5 find proper kmer match
        3 -> both flank5 and flank3 find proper kmer match
        4 -> the kmer position is opposit, which indicate the wrong kmer-matching

  @param  read_seq    the read sequence [char *].
  @param  locus       Pointer to the hash table [locus_t *].
  @param  mis_match   number of mismath allowed [int].
  @return             location index of 5' and 3' flanks on the read sequence [loc_t *]
 */
static loc_t flank_query(char *read_seq, locus_t *locus, int mis_match)
{
    loc_t loc_idx = {0}, loc = {0};
    list_t *slice_list;
    index_t *fidx;

    /* get the index of flank5 */
    fidx = locus->index5;
    slice_list = get_hit(read_seq, fidx->idx, fidx->kmer);
    if (x_size(slice_list) == 0) return loc_idx; /* no kmer match for flank5 */
 
    qsort(x_get(slice_list), x_size(slice_list), sizeof(hit_t), _compare_descend);
    loc = flank_check(read_seq, fidx->flank, slice_list, fidx->kmer+4, mis_match);
    if (!loc.status)
        return loc_idx; /* failed to located the index at the flank5 */
    else {
        loc_idx.status += 1;
        loc_idx.loc5[0] = loc.loc5[0]; 
        loc_idx.loc5[1] = loc.loc5[1];
    }
    x_destroy(slice_list);

    /* get the index of flank3 */
    fidx = locus->index3;
    slice_list = get_hit(read_seq, fidx->idx, fidx->kmer);
    if (x_size(slice_list) == 0) return loc_idx; /* no kmer match for flank3 */
 
    qsort(x_get(slice_list), x_size(slice_list), sizeof(hit_t), _compare_ascend);
    loc = flank_check(read_seq, fidx->flank, slice_list, fidx->kmer+4, mis_match);
    if (!loc.status)
        return loc_idx; /* failed to located the index at the flank3 */
    else { /* only the memory of loc.loc5 is used */
        loc_idx.status += 2;
        loc_idx.loc3[0] = loc.loc5[0];
        loc_idx.loc3[1] = loc.loc5[1];
    }
    x_destroy(slice_list);

    /* check whether the match positon is right */
    if (loc_idx.loc5[1] >= loc_idx.loc3[0]) loc_idx.status = 4;

    return loc_idx;
}


/*! @function: get allele size from read sequence by 'k-mers' algorithm.
  @param  read_seq    read sequence of the alignment [char *].
  @param  locus       Pointer to the hash table [locus_t *].
  @param  mis_match   number of mismath allowed [int].
  @return             size of allele deteted [float]
 */
float kmer_allele(char *read_seq, locus_t *locus, int mis_match)
{
    /*                    STR region
     * +++++++++++++===========================++++++++++++++
     *   flank5'    |                         |    flank3'
     *           ins_start                 ins_end
    */
    loc_t loc_idx;
    uint32_t ins_start, ins_end, allele_len;
    uint32_t allele_pre, allele_rest;

    loc_idx = flank_query(read_seq, locus, mis_match);
    if (loc_idx.status != 3) return -1.0; /* failed to located the index at flank5 or flank3 */

    ins_start = loc_idx.loc5[1] - locus->motif_len + 1;
    ins_end = loc_idx.loc3[0] + locus->motif_len;

    /* calculate the size of STR allele */
    allele_len = ins_end - ins_start - locus->exclude_base;
    allele_pre = allele_len / locus->motif_len;
    allele_rest = allele_len % locus->motif_len;

    return allele_pre + 0.1 * allele_rest;
}