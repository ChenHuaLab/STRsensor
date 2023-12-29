/*************************************************************************
    > File Name: index.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月28日 星期三 17时00分18秒
 ************************************************************************/

#include "index.h"


static uint32_t get_unit_count(char *full_seq, char *sub_seq)
{
    uint32_t unit_num = 0;
    char buf[MaxFlank<<1];
    
    uint32_t sub_len = strlen(sub_seq);
    uint32_t scope = strlen(full_seq) - sub_len;

    for (int i=0; i <= scope; ++i) {
        strncpy(buf, full_seq+i, sub_len); buf[sub_len] = '\0';
        if (!strcmp(buf, sub_seq)) unit_num += 1;
    }
    return unit_num;
}


static uint32_t get_min_kmer(char *flank_seq, char *full_seq, int flank_type, uint32_t index_kmer)
{
    char subseq[MaxFlank<<1];
    uint32_t max_len = strlen(flank_seq);
    uint32_t count;

    if (flank_type == FLANK5) {
        for (int i=index_kmer; i <= max_len; ++i) {
            strncpy(subseq, flank_seq + (max_len-i), i); subseq[i] = '\0';
            count = get_unit_count(full_seq, subseq);
            if (count == 1) return i;
        }
    }
    else { /* FLANK3 */
        for (int i=0; i < max_len-index_kmer; ++i) {
            strncpy(subseq, flank_seq, i+index_kmer); subseq[i+index_kmer] = '\0';
            count = get_unit_count(full_seq, subseq);
            if (count == 1) return index_kmer + i;
        }
    }
    return max_len;
}


static khash_t(index) *get_unique_index(char *flank_seq, char *full_seq, uint32_t min_kmer)
{
    uint32_t key_num;
    char key[MaxFlank<<1];
    
    khint_t k;
    int absent;
    khash_t(index) *hash = kh_init(index);
    uint32_t scope = strlen(flank_seq) - min_kmer;

    for (int i=0; i <= scope; ++i) {
        strncpy(key, flank_seq+i, min_kmer); key[min_kmer] = '\0';
        key_num = get_unit_count(full_seq, key);

        /* "key_num = 1" is defined as local unique kmer sequence */
        if (key_num != 1) continue;

        k = kh_put(index, hash, key, &absent);
        if (absent == 1) { /* the key is not existed in the hash-table */
            kh_key(hash, k) = x_strncopy(key, strlen(key)); kh_value(hash, k) = i;
        }
        else {
            /* all the key are unique for key_num = 1 */
            fprintf(stderr, "[index/get_unique_index] Error: No possible to be here!\n");
            exit(-2);
        }
    }
    return hash;
}


index_t *flanks_index(char *chr, uint32_t start, uint32_t end, fasta_t *fa_obj, 
            uint32_t motif_len, uint32_t index_kmer, int flank_type)
{
    index_t *index;
    uint32_t pos_s, pos_e;

    err_calloc(index, 1, index_t);
    char *full_seq = get_seq(chr, start-ShiftScope, end+ShiftScope, fa_obj);

    if (flank_type == FLANK5) {
        pos_s = start - MaxFlank;
        pos_e = start + motif_len - 1;
    }
    else { /* FLANK3 */
        pos_s = end - motif_len + 1;
        pos_e = end + MaxFlank;
    }

    index->flank = get_seq(chr, pos_s, pos_e, fa_obj);
    index->kmer = get_min_kmer(index->flank, full_seq, flank_type, index_kmer);
    index->idx = get_unique_index(index->flank, full_seq, index->kmer);
    free(full_seq);

    return index;
}


