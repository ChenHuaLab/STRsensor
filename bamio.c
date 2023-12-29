/*************************************************************************
    > File Name: bamio.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2020年11月04日 星期三 22时07分12秒
 ************************************************************************/

#include "bamio.h"
#include "utils.h"


void xbam_seq(bam1_t *aln, kstring_t *ks)
{
    int base_len;
    uint8_t *raw_seq;

    base_len = aln->core.l_qseq;
    raw_seq = bam_get_seq(aln);

    ks_resize(ks, base_len + 1);
    for (int i=0; i < base_len; ++i) {
        ks->s[i] = HTSLIB_BASE_TABLE[bam_seqi(raw_seq, i)];
    }
    ks->s[base_len] = '\0';
}


void xbam_cigar(bam1_t *aln, xcigar_t *cigar)
{
    int new_size;
    uint32_t *data;

    new_size = aln->core.n_cigar;
    
    if (cigar->m < new_size) {
        kroundup32(new_size);
        err_realloc(cigar->cigar_op, new_size, char);
        err_realloc(cigar->cigar_len, new_size, int);
        cigar->m = new_size;
    }
    cigar->n = aln->core.n_cigar;
    data = bam_get_cigar(aln);

    for (int i=0; i < cigar->n; ++i) {
        cigar->cigar_op[i] = bam_cigar_opchr(data[i]);
        cigar->cigar_len[i] = bam_cigar_oplen(data[i]);
    }
}


xcigar_t *xbam_cigar_init(int n_cigar)
{
    xcigar_t *cigar;

    err_calloc(cigar, 1, xcigar_t);
    err_calloc(cigar->cigar_op, n_cigar, char); 
    err_calloc(cigar->cigar_len, n_cigar, int);
    cigar->n = cigar->m = n_cigar;

    return cigar;
}


xbam_t *xbam_init(char *bam_file)
{
    xbam_t *bam;

    err_calloc(bam, 1, xbam_t);
    bam->aln = bam_init1(); /* bam alignment initiate */

    /* open bam file */
    bam->bam_hd = hts_open(bam_file, "rb");
    if (bam->bam_hd == NULL) {
        fprintf(stderr, "[bamio/xbam_init] Error: failed to open file of %s\n", bam_file);
        return NULL;
    }

    /* load bam header */
    bam->header = sam_hdr_read(bam->bam_hd);
    if (bam->header == NULL) {
        fprintf(stderr, "[bamio/xbam_init] Error: failed to load bam header of %s\n", bam_file);
        return NULL;
    }

    /* load index of bam file */
    bam->bam_idx = sam_index_load(bam->bam_hd, bam_file);
    if (bam->bam_idx == 0) {
        fprintf(stderr, "[bamio/xbam_init] Error: failed to load bam index of %s.bai\n", bam_file);
        return NULL;
    }

    return bam;
}


int xbam_fetch(char *chr, int start, int end, xbam_t *bam)
{
    char position[0x40];

    sprintf(position, "%s:%d-%d", chr, start, end);
    bam->iter[bam->n_iter++] = sam_itr_querys(bam->bam_idx, bam->header, position);

    if (bam->iter[bam->n_iter-1] == NULL) {
        bam->n_iter--;
        fprintf(stderr, "[bamio/xbam_fetch] Error: failed to fetch position of %s\n", position);
        return -1;
    }
    return 0;
}


void xbam_destroy(xbam_t *bam)
{
    bam_hdr_destroy(bam->header); /* destroy the header structure */
    bam_destroy1(bam->aln); /* destroy the alignment structure */
    hts_idx_destroy(bam->bam_idx); /* destroy the bam index structure */
    hts_close(bam->bam_hd); /* close the bam handle */
    for (int i=0; i < bam->n_iter; ++i) hts_itr_destroy(bam->iter[i]); /* destroy iters */
}


void xbam_cigar_destroy(xcigar_t *cigar)
{
    free(cigar->cigar_op); /* destroy the cigar operator */
    free(cigar->cigar_len); /* destroy the cigar length */
    cigar->n = cigar->m = 0; /* I just want to do it, HaHaHa */
}


