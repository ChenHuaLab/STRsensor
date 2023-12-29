/*************************************************************************
    > File Name: getseq.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2018年05月21日 星期一 17时12分12秒
 ************************************************************************/

#include "getseq.h"
#include "utils.h"


static faidx_t *fai_load(char *faifile)
{
    faidx_t *fai;
    char buf[0x1000], chr[0x100];

    err_calloc(fai, 1, faidx_t);
    err_open(fai->fai_fp, faifile, "r");

    while (fgets(buf, 0x1000, fai->fai_fp)) {
        if (fai->n % 0x20 == 0) {
            err_realloc(fai->line_fai, fai->n+0x20, faidx1_t);
            err_realloc(fai->name, fai->n+0x20, char *);
        }
        faidx1_t *f = &fai->line_fai[fai->n];
        sscanf(buf, "%s%"PRId64"%"PRIu64"%"PRId32"%"PRId32,
            chr, &f->len, &f->offset, &f->line_blen, &f->line_len);
        fai->name[fai->n++] = x_strncopy(chr, strlen(chr));
    } fclose(fai->fai_fp);

    return fai;
}


static int32_t get_chr(char **name, int32_t n, char *chr)
{
    int32_t i;

    for (i=0; i < n; ++i) {
        if (!strcmp(name[i], chr))
            return i;
    }
    return -1;
}


fasta_t *fa_init(char *fa_file)
{
    fasta_t *fa_obj;
    char buf[0x1000];

    err_calloc(fa_obj, 1, fasta_t);
    sprintf(buf, "%s%s", fa_file, ".fai");

    if (!fopen(buf, "r")) {
        fprintf(stderr, "[get_seq:fasta_init] Error: failed to open index of %s\n", fa_file);
        exit(-1);
    }

    fa_obj->fai = fai_load(buf);
    err_open(fa_obj->fa_fp, fa_file, "r");
    fa_obj->fa_file = x_strncopy(fa_file, strlen(fa_file));

    return fa_obj;
}


char *get_seq(char *chr, int32_t start, int32_t end, fasta_t *fa_obj)
{
    char *seq;
    uint64_t d_offset;

    faidx_t *fai = fa_obj->fai;
    FILE *fa = fa_obj->fa_fp;
    int32_t index = get_chr(fai->name, fai->n, chr);

    if (index == -1 || end < start) {
        fprintf(stderr, "Err: Invaild bed line %s:%d:%d\n", chr,start,end);
        return NULL;
    }
    faidx1_t *f = &fai->line_fai[index];
    int32_t s = f->line_len - f->line_blen;
    d_offset = f->offset + (start-1) + (start-1)/f->line_blen*s;

    /* allocate necessary memeory and move the pointer handle */
    err_calloc(seq, end-start+2, char);
    fseek(fa, d_offset, SEEK_SET);

    int num =0;
    for (int ch=0; (ch = getc(fa)) != EOF; ) {
        if (num >= end-start+1) break;
        if (ch != '\r' && ch != '\n') seq[num++] = ch;
    } seq[num] = '\0';

    return seq;
}

