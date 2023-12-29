/*************************************************************************
    > File Name: model.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月31日 星期六 20时02分18秒
 ************************************************************************/

#include "model.h"
#include "utils.h"


static int allele_idx_search(cal1_t *cal, int n, float allele)
{
    int idx = -1;

    for (int i=0; i < n; ++i) {
        if (fabs(cal[i].allele - allele) < 0.01) { idx = i; break; }
    }
    return idx;
}


void allele_count(cal_t *cal, float allele)
{
    int idx;

    idx = allele_idx_search(cal->cal, cal->n, allele);
    if (idx != -1) { /* the allele existed in cal structure */
        cal->cal[idx].nf++; 
        return ;
    }
    if (cal->n == cal->m) { /* realloc memory if necessary */
        cal->m = cal->m ? cal->m + 8 : 8;
        cal->cal = (cal1_t*)realloc(cal->cal, sizeof(cal1_t) * cal->m);
    }
    cal->cal[cal->n].allele = allele;
    cal->cal[cal->n++].nf = 1;
}


static void freq_push(cal_t *cal, float allele, float freq)
{
    int idx;

    idx = allele_idx_search(cal->cal, cal->n, allele);
    if (idx != -1) return ; /* the allele existed in cal structure */

    if (cal->n == cal->m) { /* realloc memory if necessary */
        cal->m = cal->m ? cal->m + 8 : 8;
        cal->cal = (cal1_t*)realloc(cal->cal, sizeof(cal1_t) * cal->m);
    }
    cal->cal[cal->n].allele = allele;
    cal->cal[cal->n++].nf = freq;
}


static void haplotype_stutter(float *a_list, int n_allele, float allele, int *stutter)
{
    float diff;
    int idx;

    for (int i=0; i < n_allele; ++i) {
        diff = fabs(a_list[i] - allele);
        if (diff == 0.0||diff==1.0||diff==2.0) {
            idx = a_list[i] < allele ? (int)(-1*diff) + 5 : (int)diff;
            stutter[idx]++;
        }
    } return ;
}


static int get_close_diff(float a, float *a_order)
{
    float diff1, diff2;

    diff1 = a - a_order[0]; /* left length*/
    diff2 = a_order[1] - a; /* right length */

    if (diff2==diff1 && (diff2==1.0 || diff2==2.0))
        return (int)(-1*diff2) + 5;

    if (diff1 < diff2 && (diff1==1.0 || diff1==2.0))
        return (int)diff1;

    if (diff1 > diff2 && (diff2==1.0 || diff2==2.0))
        return (int)(-1*diff2) + 5;

    return -1; /* other confused condition */
}


static void autosome_stutter(float *a_list, int n_allele, float *alleles, int *stutter)
{
    float a_order[2];
    float diff;
    int idx;

     /* sort the allele by ascending order */
    if (alleles[0] > alleles[1]) { 
        a_order[0] = alleles[1]; a_order[1] = alleles[0]; 
    }
    else memcpy(a_order, alleles, 2 * sizeof(float));
    
    /* determine the stutter by the given alleles */
    for (int i=0; i < n_allele; ++i) {
        if (a_list[i]==a_order[0] || a_list[i]==a_order[1]) stutter[0]++; /* no stutter */
        else if (a_list[i] < a_order[0]) { /* stutter down */
            diff = a_list[i] - a_order[0];
            if (diff == -1.0 || diff == -2.0) stutter[(int)diff+5]++;
        }
        else if (a_list[i] > a_order[1]) { /* stutter up */
            diff = a_list[i] - a_order[1];
            if (diff == 1.0 || diff == 2.0) stutter[(int)diff]++;
        }
        else { /* within the allele region */
            idx = get_close_diff(a_list[i], a_order);
            if (idx != -1) stutter[idx]++;
        }
    } return ;
}


static void person_process(float *a_list, int n_allele, evaluate_t *eva, int min_reads, int is_haplotype)
{
    cal_t cal = {0}; /* allele (from spanning reads) distribution */

    memset(eva, 0, sizeof(evaluate_t));
    for (int i=0; i < n_allele; ++i) allele_count(&cal, a_list[i]);
    if (cal.n == 0 || n_allele < min_reads) { free(cal.cal); return ; }

    /* Rules: (1) totoal span reads > 10  (2) top1 ratio > 0.7) */
    qsort(cal.cal, cal.n, sizeof(cal1_t), _allele_descend);
    if (is_haplotype) {
        if (cal.cal[0].nf / n_allele > 0.7 ) {
            eva->allele[0] = cal.cal[0].allele;
            haplotype_stutter(a_list, n_allele, eva->allele[0], eva->eva_stu);
            eva->prob = 1.0; eva->status = 1; 
        }
        free(cal.cal); return ;
    }
    float top1_num = cal.cal[0].nf;
    float top2_num = cal.n > 1 ? cal.cal[1].nf : cal.cal[0].nf;
    eva->allele[0] = cal.cal[0].allele;

    /* 1. homozygous allele */
    if (cal.n==1 || (top2_num/top1_num < 0.3 && top1_num/n_allele > 0.7))
        eva->allele[1] = cal.cal[0].allele;
    /* 2. heterozygous allele */
    else if (top2_num/top1_num >= 0.3 && (top1_num+top2_num)/n_allele > 0.7)
        eva->allele[1] = cal.cal[1].allele;
    /* 3. disorganized allele distribution */
    else { free(cal.cal); return ; }

    autosome_stutter(a_list, n_allele, eva->allele, eva->eva_stu);
    eva->prob = 1.0; eva->status = 1; free(cal.cal);
}


static void evaluate_freq(cal_t *in_freq, cal_t *cal_freq, int n_idv)
{
    /* evaluate the parameters of frequency */
    if (in_freq->n != 0) { /* frequency is given by user */
        for (int i=0; i < cal_freq->n; ++i) 
            freq_push(in_freq, cal_freq->cal[i].allele, EPS);
    }
    else { /* frequency of this locus is not given */
        if (n_idv < 20) { 
            for (int i=0; i < cal_freq->n; ++i)
                freq_push(in_freq, cal_freq->cal[i].allele, 1.0);
        }
        else { /* abundant samples available */
            int allele_num=0;
            for (int i=0; i < cal_freq->n; ++i) allele_num += cal_freq->cal[i].nf;
            for (int i=0; i < cal_freq->n; ++i)
                freq_push(in_freq, cal_freq->cal[i].allele, cal_freq->cal[i].nf/allele_num);
        }
    } return ;
}


static void evaluate_stutter(stutter1_t *in_stu, int *cal_stu, int n_idv)
{
    /* evaluate the parameters of stutter */
    int span_reads = 0;
    float value = 0;
    float basic_stu[5] = {0.85, 0.03, 0.02, 0.03, 0.07};

    if (in_stu->num == 0) { /* stutter is not given by user */
        for (int i=0; i < 5; ++i) span_reads += cal_stu[i];
        if (span_reads >= 100) {
            for (int i=0; i < 5; ++i) {
                value = cal_stu[i] / (float)span_reads;
                in_stu->s[i] = value > EPS ? value : EPS;
            }
        } else memcpy(in_stu->s, basic_stu, 5 * sizeof(float));
        in_stu->num = 5;
    }
    /* pre-calculate the value of log10(p_stutter) */
    float v1, v2;
    
    for (int i=0; i < 6; ++i) {
        if (i != 5) /* for log10(stutter) */
            in_stu->slog[i] = log10(in_stu->s[i]);

        v1 = (i==5) ? EPS : in_stu->s[i];
        for (int j=0; j < 6; ++j) { /* for log10(s1+s2) */
            v2 = (j==5) ? EPS : in_stu->s[j];
            in_stu->sslog[i][j] = log10(v1 + v2);
        }
    } return ;
}


void person_destroy(person_t *person)
{
    person1_t *idv;
    int n_idv;

    n_idv = person->n_idv;

    for (int i=0; i < n_idv; ++i) {
        idv = &person->idv[i]; 
        free(idv->sample); free(idv->allele_list);
    }
    free(person->idv); free(person);
}


void params_destroy(params_t *params)
{
    free(params->f_loci->cal); /* allele frequency */
    free(params->s_loci); /* stutter */
    free(params);
}


cal_t *freq_locus_search(freq_t *freq, char *locus)
{
    cal_t *cal;

    for (int idx=0; idx < freq->n; ++idx) {
        if (strcmp(locus, freq->core[idx].name) == 0) 
            return &freq->core[idx];
    }
    return NULL;
}


stutter1_t *stu_locus_search(stutter_t *stutter, char *locus)
{
    stutter1_t *stu;

    for (int idx=0; idx < stutter->n; ++idx) {
        if (strcmp(locus, stutter->core[idx].name) == 0) 
            return &stutter->core[idx];
    }
    return NULL;
}


params_t *locus_model(person1_t *idvs, int n_idv, cal_t *f_loci, stutter1_t *s_loci, int min_reads, int is_haplotype)
 {
     cal_t cal_freq = {0}; /* allele (from each person) distribution */
     int cal_stu[5] = {0}; /* 0(N), 1(add1), 2(add2), 3(del2), 4(del1) */
     evaluate_t eva;

    for (int i=0; i < n_idv; ++i) { /* thread unsafe loop */
        person1_t *idv = &idvs[i];
        if (idv->n == 0) continue; /* no available allele */

        person_process(idv->allele_list, idv->n, &eva, min_reads, is_haplotype); 
        if (eva.status == 0) continue; /* no high confidence allele call */

        if (is_haplotype) 
            allele_count(&cal_freq, eva.allele[0]);
        else { /* autosome chromosome */
            allele_count(&cal_freq, eva.allele[0]); /* first allele */
            allele_count(&cal_freq, eva.allele[1]); /* second allele */
        }
        for (int j=0; j < 5; ++j) cal_stu[j] += eva.eva_stu[j];
    }
    params_t *param;
    err_calloc(param, 1, params_t);
    param->is_haplotype = is_haplotype;

    /* evaluate parameters with both input and called params */
    evaluate_freq(f_loci, &cal_freq, n_idv);
    evaluate_stutter(s_loci, cal_stu, n_idv);
    param->f_loci = f_loci; param->s_loci = s_loci;
    free(cal_freq.cal); /* deallocate cal_freq.cal */

    return param;
 }


static float haplotype_prob(float allele, float *a_list, int n_allele, stutter1_t *stu, float c_prob)
{
    float like_value = 0.0;
    float diff;
    int idx;

    if (n_allele == 0) /* no available allele in the a_list */
        return like_value;

    for (int i=0; i < n_allele; ++i) {
        diff = fabs(a_list[i] - allele);
        if (diff == 0.0||diff==1.0||diff==2.0) {
            idx = a_list[i] < allele ? (int)(-1*diff) + 5 : (int)diff;
            like_value += stu->slog[idx];
        }
        else like_value += -3.0; /* log10(EPS) */ 
    }
    return log10(c_prob) + like_value;
}


static float autosome_prob(float allele1, float allele2, float *a_list, int n_allele, stutter1_t *stu, float c_prob)
{
    float like_value = 0.0;
    float diff;
    int idx1, idx2;

    if (n_allele == 0) /* no available allele in the a_list */
        return like_value;

    for (int i=0; i < n_allele; ++i) { /* calculate the likelihood */
        /* for the 1-st allele */
        diff = fabs(a_list[i] - allele1);
        if (diff == 0.0||diff==1.0||diff==2.0) 
            idx1 = a_list[i] < allele1 ? (int)(-1*diff) + 5 : (int)diff;
        else idx1 = 5;

        /* for the 2-st allele */
        diff = fabs(a_list[i] - allele2);
        if (diff == 0.0||diff==1.0||diff==2.0) 
            idx2 = a_list[i] < allele2 ? (int)(-1*diff) + 5 : (int)diff;
        else idx2 = 5;

        like_value += stu->sslog[idx1][idx2];
    }
    return log10(c_prob) + like_value;
}


static allele_t *gen_candidate(float *a_list, int n_allele, int is_haplotype, int is_amplicon)
{
    allele_t *c_list;
    cal_t cal = {0};
    int n_top, n_comb;

    err_calloc(c_list, 1, allele_t);
    for (int i=0; i < n_allele; ++i) allele_count(&cal, a_list[i]);
    qsort(cal.cal, cal.n, sizeof(cal1_t), _allele_descend);
    n_top = cal.n > 4 ? 4 : cal.n; /* get the top 4 alleles */

    /* 1. haplotype condition (is_haplotype = 1) */
    if (is_haplotype) {
        c_list->n = n_top;
        err_calloc(c_list->core, n_top, allele1_t);
        for (int i=0; i < n_top; ++i) c_list->core[i].allele1 = cal.cal[i].allele;
        free(cal.cal); return c_list;
    }

    /* 2. autosome condition (is_haplotype = 0) */
    allele1_t *cur;
    float ratio;
    
    n_comb = n_top * (n_top + 1) / 2; /* maximum number of combinations */
    err_calloc(c_list->core, n_comb, allele1_t);

    for (int i=0; i < n_top; ++i) {
        for (int j=i; j < n_top; ++j) {
            ratio = cal.cal[i].nf / cal.cal[j].nf;
            cur = &c_list->core[c_list->n];

            if (is_amplicon) { /* like amplicon/target seuening */
                if (ratio > 0.3 && ratio < 3.3) {
                    cur->allele1 = cal.cal[i].allele; cur->allele2 = cal.cal[j].allele;
                }
                else continue; /* skip this candidate allele */
            }
            else { cur->allele1 = cal.cal[i].allele; cur->allele2 = cal.cal[j].allele; }
            c_list->n++;
        }
    } free(cal.cal); return c_list;
}


static void gen_allele(allele_t *c_list)
{
    float diff, value_sum = 0.0;
    float max_value = (float)(-1 << 30);

    for (int i=0; i < c_list->n; ++i) 
        max_value = c_list->core[i].prob > max_value ? c_list->core[i].prob : max_value;

    for (int i=0; i < c_list->n; ++i) {
        diff = c_list->core[i].prob - max_value;
        if (diff > -4.0) {
            c_list->core[i].prob = pow(10.0, diff);
            value_sum += c_list->core[i].prob;
        }
        else c_list->core[i].prob = 0.0;
    }

    for (int i=0; i < c_list->n; ++i) c_list->core[i].prob /= value_sum;
    qsort(c_list->core, c_list->n, sizeof(allele1_t), _prob_descend);
    
    return ;
}


allele1_t get_allele(float *a_list, int n_allele, params_t *params, int min_reads, int is_amplicon)
{
    allele1_t decide_allele = {0}, *a_cur;
    allele_t *c_list; /* c_list: candidate list */

    if (n_allele < min_reads) 
        return decide_allele; /* not enough reads to genotype */

    int idx; /* allele index in the freq */
    float p1, p2, lk_value;

    c_list = gen_candidate(a_list, n_allele, params->is_haplotype, is_amplicon);
    for (int i=0; i < c_list->n; ++i) {
        a_cur = &c_list->core[i]; /* current allele structure */
        idx = allele_idx_search(params->f_loci->cal, params->f_loci->n, a_cur->allele1);
        p1 = (idx==-1) ? 1.0 : params->f_loci->cal[idx].nf;

        if (params->is_haplotype) {
            lk_value = haplotype_prob(a_cur->allele1, a_list, n_allele, params->s_loci, p1);
        }
        else { /* autosome chromosome */
            idx = allele_idx_search(params->f_loci->cal, params->f_loci->n, a_cur->allele2);
            p2 = (idx==-1) ? 1.0 : params->f_loci->cal[idx].nf;
            lk_value = autosome_prob(a_cur->allele1, a_cur->allele2, a_list, n_allele, params->s_loci, p1 * p2);
        }
        a_cur->prob = lk_value;
    }
    gen_allele(c_list); /* get the most likely allele */
    memcpy(&decide_allele, &c_list->core[0], sizeof(allele1_t));
    free(c_list->core);

    return decide_allele;
}