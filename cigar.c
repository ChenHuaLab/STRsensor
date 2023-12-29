/*************************************************************************
    > File Name: cigar.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年11月19日 星期四 09时29分20秒
    > Updated Time: 2021年01月05日 星期二 16时45分11秒
 ************************************************************************/

#include "share.h"
#include "cigar.h"
#include "htslib/kstring.h"


/*! @function: check whether the indel is an allele-near indel
 * @param indel         one indel1_t structure, like (ref_pos, cigar_clen, cigar_type, cigar_value)
 * @param read_seq      the read sequence
 * @param locus         locus object [locus_t *]
 * @return              normal indel -> 0;  shift indel -> shift_bp
*/
static int cigar_shift(indel1_t *indel, char *read_seq, locus_t *locus)
{
    
    int ind_start, ind_end;

    ind_start = indel->ref_pos;
    ind_end = indel->ref_pos + indel->c_value;

    /* indel occured within the allele region */
    if (ind_start >= locus->start && ind_end <= locus->end)
        return 0;

    /* treat as normal indel */
    if (indel->c_value % locus->motif_len != 0)
        return 0;

    int shift_bp = 0;
    int safe_dist = 2 * locus->motif_len;

    /* indel before STR region (major condition) */
    if (ind_start < locus->start && ind_start > locus->start-safe_dist) {
        shift_bp = locus->start - ind_start;
        return shift_bp;
    }
    /* indel before STR region (minor condition) */
    if (ind_end > locus->end && ind_end < locus->end+safe_dist) {
        shift_bp = locus->end - ind_end;
        return shift_bp;
    }
    
    return shift_bp;
}


/*! @function: move the str unit (diff from ref) to the allele region
 * @param r_cigar       the list of cigar value, e.g. 40M4I106M -> [('M',40), ('I',4), ('M',106)]
 * @param cigar_start   the reference start of the alignment (convert to 1-based)
 * @param read_seq      the read sequence
 * @param locus         locus object [locus_t *]
 * @return              treated cigar list [xcigar_t *]
*/
static xcigar_t *cigar_process(xcigar_t *r_cigar, int64_t cigar_start, char *read_seq, locus_t *locus)
{
    xcigar_t *t_cigar; /* treated cigar */
    indel_t w_cigar; /* new cigar */
    int cigar_rlen=0, cigar_clen=0, indel_num=0;
    
    w_cigar.n = r_cigar->n;
    err_calloc(w_cigar.core, r_cigar->n, indel1_t);
    t_cigar = xbam_cigar_init(r_cigar->n);
    
    indel1_t *cc;
    for (int i=0; i < r_cigar->n; ++i) { /* construct the new cigar list */
        cc = &w_cigar.core[i];
        cc->ref_pos = cigar_start + cigar_rlen; 
        cc->c_clen = cigar_clen;
        cc->c_type = r_cigar->cigar_op[i]; 
        cc->c_value = r_cigar->cigar_len[i];

        if (r_cigar->cigar_op[i] != 'D') cigar_clen += r_cigar->cigar_len[i];
        if (strchr("MD=X", r_cigar->cigar_op[i])) cigar_rlen += r_cigar->cigar_len[i];
        if (strchr("ID", r_cigar->cigar_op[i])) indel_num++;
    }
    if (indel_num == 0) { /* NO indel occured */
        t_cigar->n = r_cigar->n;
        memcpy(t_cigar->cigar_op, r_cigar->cigar_op, t_cigar->n * sizeof(char));
        memcpy(t_cigar->cigar_len, r_cigar->cigar_len, t_cigar->n * sizeof(int));
        free(w_cigar.core); return t_cigar;
    }
    /* processing the indel(I/D) in the cigar value */
    int shift_bp;
    int left_len, right_len;

    for (int i=0; i < w_cigar.n; ++i) {
        cc = &w_cigar.core[i];
        if (!strchr("ID", cc->c_type)) { /* not Insert(I) or Delete(D) */
            t_cigar->cigar_op[i] = cc->c_type; t_cigar->cigar_len[i] = cc->c_value;
            continue;
        }
        shift_bp = cigar_shift(cc, read_seq, locus);
        if (shift_bp == 0) { /* no need to shift the indel */
            t_cigar->cigar_op[i] = cc->c_type; t_cigar->cigar_len[i] = cc->c_value;
        }
        else { /* only one real 'allele-shift' is possible */
            for (int j=i; j < w_cigar.n; ++j) {
                t_cigar->cigar_op[j] = w_cigar.core[j].c_type; 
                t_cigar->cigar_len[j] = w_cigar.core[j].c_value;
            }
            left_len = t_cigar->cigar_len[i-1]+shift_bp;
            right_len = t_cigar->cigar_len[i+1]-shift_bp;
            if (left_len > 0 && right_len > 0) {
                t_cigar->cigar_len[i-1] = left_len; t_cigar->cigar_len[i+1] = right_len;
            }
            break; /* the remaining cigar has been copied to t_cigar */
        }
    } free(w_cigar.core); return t_cigar;
}


/*! @function: locate the allele index at the read sequence
 * @param t_cigar          processed treated cigar list
 * @param cigar_start      the reference start of the alignment (convert to 1-based)
 * @param l_start          STR locus start (1-based)
 * @param l_end            STR locus end (1-based)
 * @return                 allele start/end index on the read (0-based)
*/
static locate_t locate_index(xcigar_t *t_cigar, int64_t cigar_start, int64_t l_start, int64_t l_end)
{
    locate_t locate={0, -1, -1};
    char *cigar_string;
    int cigar_rlen=0, cigar_tlen=0, tlen=0;

    /* prepare the cigar string */
    for (int i=0; i < t_cigar->n; ++i) tlen += t_cigar->cigar_len[i];
    err_malloc(cigar_string, tlen+1, char); cigar_string[tlen] = '\0';

    for (int i=0; i < t_cigar->n; ++i) {
        if (strchr("MD=X", t_cigar->cigar_op[i])) cigar_rlen += t_cigar->cigar_len[i];
        memset(cigar_string+cigar_tlen, t_cigar->cigar_op[i], t_cigar->cigar_len[i]);
        cigar_tlen += t_cigar->cigar_len[i];
    }
    /* check the allele start and end */
    if (l_start < cigar_start || l_end >= cigar_start+cigar_rlen) {
        free(cigar_string); return locate;
    }

    /* get the start index allele in the read */
    int64_t pos=cigar_start;
    int cigar_sidx=0, last_midx=0; /* cigar start index & last match index */

    while (pos < l_start) {
        if (strchr("MD=X", cigar_string[cigar_sidx])) pos++;
        if (strchr("M=X", cigar_string[cigar_sidx])) last_midx = cigar_sidx;
        cigar_sidx++;
    }
    cigar_sidx = last_midx;
    if (cigar_sidx == 0 && (!strchr("M=X", cigar_string[cigar_sidx]))) { /* not in 'M=X' */
        free(cigar_string); return locate;
    }

    /* get the end index allele in the read */
    int cigar_eidx;
    pos = cigar_start + cigar_rlen;
    cigar_eidx = last_midx = cigar_tlen - 1;

    while (pos > l_end) {
        if (strchr("MD=X", cigar_string[cigar_eidx])) pos--;
        if (strchr("M=X", cigar_string[cigar_eidx])) last_midx = cigar_eidx;
        if (cigar_eidx == 0) break;
        cigar_eidx--;
    }
    cigar_eidx = last_midx;
    if (cigar_eidx == cigar_tlen-1 && (!strchr("M=X", cigar_string[cigar_eidx]))) {
        free(cigar_string); return locate;
    }

    /* calculate the delete bases on cigarstring if D existed */
    int d_left = 0; // normal deletion before STR region
    int d_core = 0; // usually motif (deletion)
    if (cigar_sidx >= cigar_eidx) { free(cigar_string); return locate; }
    for (int i=0; i <= cigar_sidx; ++i) d_left += (cigar_string[i] == 'D') ? 1 : 0;
    for (int i=cigar_sidx+1; i <= cigar_eidx; ++i) d_core += (cigar_string[i] == 'D') ? 1 : 0;

    /* given the final STR index on the read sequence */
    locate.status = 1;
    locate.loc5 = cigar_sidx-d_left+1; locate.loc3 = cigar_eidx-d_left-d_core;
    free(cigar_string); return locate;
}


/*! @function: check whether there have MIN_MATCH flank bases on both side
 * @param str1          the first sequence used to check the flank sequencing
 * @param str2          the second sequence used to check the flank sequencing
 * @param min_match     equal to MIN_MATCH + motif_len
 * @param mis_match     maximum number of mismatch bases allowed
 * @return              pass the check (1) and failed to check (0)
*/
static int match_check(char *str1, char *str2, int min_match, int mis_match)
{
    int len1, len2;
    int misnum = 0;
    int matchnum = 0;

    #ifndef ACCEPT
      #define ACCEPT 1
      #define DISCARD 0
    #endif

    len1 = strlen(str1); len2 = strlen(str2);
    if (len1 < min_match || len2 < min_match)
        return DISCARD; /* the minimun match bases need */

    for (int i=0; i < min_match; ++i) {
        if (str1[i] != str2[i]) misnum++;
        else matchnum++;
        if (misnum > mis_match) return DISCARD;
    }
    return ACCEPT;
}


/*! @function: check whether there have MIN_MATCH flank bases on both side
 * @param read_seq          the read sequence used to check the flank sequencing
 * @param flank5            the sequence of flank5 (length=MaxFlank)
 * @param flank3            the sequence of flank3 (length=MaxFlank)
 * @param motif_len         STR motif length. eg. ATTC -> 4
 * @param loc_idx           index on the read sequence (0-based)
 * @param mis_match         maximum number of mismatch bases allowed
 * @return                  passed:1  failed:0
*/
static int64_t flank_check(char *read_seq, char *flank5, char *flank3, int motif_len, locate_t *loc_idx, int mis_match)
{
    int64_t keep_flag = 0;
    char *f5_seq, *f5_flank, *f3_seq, *f3_flank;

    if (loc_idx->loc5 < MIN_MATCH || strlen(read_seq+loc_idx->loc3) < MIN_MATCH)
        return keep_flag; /* not enough bases */

    f5_seq = read_seq + (loc_idx->loc5 - MIN_MATCH);
    f5_flank = flank5 + (strlen(flank5) - motif_len - MIN_MATCH);
    keep_flag = match_check(f5_seq, f5_flank, motif_len+MIN_MATCH, mis_match);
    if (keep_flag == 0) return keep_flag;

    f3_seq = read_seq + (loc_idx->loc3 - motif_len + 1);
    f3_flank = flank3;
    keep_flag = match_check(f3_seq, f3_flank, motif_len+MIN_MATCH, mis_match);
    if (keep_flag == 0) return keep_flag;

    return 1; /* Successfully PASSED */
}


/*! @function: get allele by cigar value provided by alignment
 * @param read_seq          the read sequence used to check the flank sequencing
 * @param r_cigar           the raw cigar list parsed from alignment
 * @param ref_pos           the reference position (1-based)
 * @param locus             STR locus object [locus_t *]
 * @return                  detected allele [float]
*/
float cigar_allele(char *read_seq, xcigar_t *r_cigar, int64_t ref_pos, locus_t *locus, int mis_match)
{
    int64_t keep_status;
    xcigar_t *t_cigar; /* t_cigar: treated cigar; r_cigar: raw cigar */
    locate_t loc_idx;

    t_cigar = cigar_process(r_cigar, ref_pos, read_seq, locus);
    loc_idx = locate_index(t_cigar, ref_pos, locus->start, locus->end);
    if (loc_idx.status == 0)
        return -1.0; /* failed to get 5' and 3' flank index */

    keep_status = flank_check(read_seq, locus->index5->flank, locus->index3->flank, locus->motif_len, &loc_idx, mis_match);
    if (keep_status == 0)
        return -1.0; /* can't make proper mathing with flank5 or flank3 */

    uint32_t allele_len;
    uint32_t allele_pre, allele_rest;

    allele_len = loc_idx.loc3 - loc_idx.loc5 - locus->exclude_base + 1;
    allele_pre = allele_len / locus->motif_len;
    allele_rest = allele_len % locus->motif_len;
    xbam_cigar_destroy(t_cigar); /* destroy the memory allocated */

    return allele_pre + 0.1 * allele_rest;
}