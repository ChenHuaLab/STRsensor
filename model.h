/*************************************************************************
    > File Name: model.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年10月31日 星期六 20时00分00秒
 ************************************************************************/

#ifndef MODEL_H
#define MODEL_H

#include <math.h>

#define EPS 0.001


/*! @typedef person1_t
 @abstract structure of person info at specific locus
 @field sample                 sample name (eg. sample1.bam)[char *]
 @field allele                 determined allele by MAP algorithm (eg. 13.2 or 12/14.3) [char *]
 @field allele_prob            allele probability (eg. 0.9888) [float]
 @field total_reads            total reads within this region [int]
 @field valid_reads            number of reads past the filter rule [int]
 @field span_reads             number of reads entirely span the allele region [int]
 @field n, m                   number of allele and max memory allocated [int]
 @field allele_list            alleles detected by span reads [float *]
*/
typedef struct person1_t {
    char *sample;
    char allele[0x10];
    float allele_prob;
    int total_reads;
    int valid_reads;
    int span_reads;
    int n, m;
    float *allele_list;
} person1_t;

/*! @typedef person_t
 @abstract structure of all person info at specific locus
 @field locus             STR locus name [char *]
 @field n_idv             number of person [int]
 @field idv               structure for each idv [person1_t]
*/
typedef struct person_t {
    char locus[0x40];
    int n_idv;
    person1_t *idv;
} person_t;


/*! @typedef sample_t
 @abstract structure of sample file path list
 @field n             number of samples [int]
 @field name          sample basename of the given file path [char **]
 @field path          file path for each given sample [char **]
*/
typedef struct sample_t {
    int n;
    char **name;
    char **path;
} sample_t;


/*! @typedef cal1_t
 @abstract structure of allele count like -> {13.2:553, 14:508, 15:25}
 @field allele           allele from spanning reads [float]
 @field nf               n: allele number; f: allele frequency [float]
*/
typedef struct cal1_t {
    float allele;
    float nf;
} cal1_t;

/*! @typedef cal_t
 @abstract shared structure of calculated allele count/allele frequency
 @field name             name of STR locus [char *]
 @field n                number of allele count group, like -> 13.2: 553 [int]
 @field m                the maximum number of memory allocated [int]
 @field cal              pointer to each allele count group [cal1_t *]
*/
typedef struct cal_t {
    char name[0x30];
    int n, m;
    cal1_t *cal;
} cal_t;


/*! @typedef freq_t
 @abstract structure for STR candidate allele frequency 
 @field n                 number of STR locus [int]
 @field m                 recorde the memory allocated freq1_t [int]
 @field core              structure for each locus [freq1_t]
*/
typedef struct freq_t {
    int n, m;
    cal_t *core;
} freq_t;


/*! @typedef stutter1_t
 @abstract structure of STR stutter parameters for each locus 
 @field name            name of STR locus [char *]
 @field num             number of stutter parmeter
 @field s               stutter ratio for 0, 1, 2, -2, -1 [float *]
 @field slog            log10 for each stutter for htplotype STR locus [float *]
 @field sslog           log10 (s1 + s2) for autosome STR locus [float **]
*/
typedef struct stutter1_t {
    char name[0x40];
    int num;
    float s[5];
    float slog[5];
    float sslog[6][6];
} stutter1_t;

/*! @typedef stutter_t
 @abstract structure of STR stutter
 @field n                 number of STR locus [int]
 @field m                 recorde the memory allocated stutter1_t [int]
 @field core              structure for each locus [stutter1_t]
*/
typedef struct stutter_t {
    int n, m;
    stutter1_t *core;
} stutter_t;


/*! @typedef params_t
 @abstract structure of params for one given locus
 @field f_loci            frequency for one STR locus [cal_t *]
 @field s_loci            stutter params for one STR locus [stutter1_t *]
 @field core              the locus type, which could be haplotype or not [int]
*/
typedef struct params_t {
    cal_t *f_loci;
    stutter1_t *s_loci;
    int is_haplotype;
} params_t;


/*! @typedef evaluate_t
 @abstract structure of evalutated stutters and allele
 @field allele           determined alleles. [float *]
 @field prob             the probability of the evaluated stutter [floag]
 @field eva_stu          evaluated stutter by given allele that extracted from spanning reads [int *]
 @field status           success or not in evaluating the allele
*/
typedef struct evaluate_t {
    float allele[2];
    float prob;
    int eva_stu[5];
    int status;
} evaluate_t;


/*! @typedef allele1_t
 @abstract structure of allele target, like -> [(11.2, 13.2), 0.937]
 @field allele1          determined allele1 [float]
 @field allele2          determined allele2 [float]
 @field prob             probability/likelihood value [float]
*/
typedef struct allele1_t {
    float allele1;
    float allele2; /* autosome chromosome */
    float prob;
} allele1_t;

/*! @typedef allele_t
 @abstract structure of allele count/allele frequency
 @field n                number of allele target [int]
 @field m                the maximum number of memory allocated [int]
 @field core             pointer to the allele target [allele1_t *]
*/
typedef struct allele_t {
    int n, m;
    allele1_t *core;
} allele_t;


static int _allele_descend(const void *a, const void *b) {
    /* sorted by descending order */
    return ((cal1_t*)b)->nf - ((cal1_t*)a)->nf;
}

static int _candidate_ascend(const void *a, const void *b) {
    /* sorted by descending order */
    return ((cal1_t*)a)->allele - ((cal1_t*)b)->allele;
}

static int _prob_descend(const void *a, const void *b) {
    /* sorte the probability by descending order */
    return ((allele1_t*)b)->prob - ((allele1_t*)a)->prob;
}


/* push an allele to allele list and allocate memory if necessary */
#define allele_push(idv, allele) do {                                                         \
    if ((idv)->n == (idv)->m) {                                                               \
        (idv)->m = (idv)->m ? (idv)->m<<1 : 128;                                              \
        (idv)->allele_list = (float*)realloc((idv)->allele_list, sizeof(float) * (idv)->m);   \
    }                                                                                         \
    (idv)->allele_list[(idv)->n++] = (allele);                                                \
} while (0)


params_t *locus_model(person1_t *idvs, int n_idv, cal_t *f_loci, stutter1_t *s_loci, int min_reads, int is_haplotype);
void person_destroy(person_t *person);
cal_t *freq_locus_search(freq_t *freq, char *locus);
stutter1_t *stu_locus_search(stutter_t *stutter, char *locus);
void allele_count(cal_t *cal, float allele);
allele1_t get_allele(float *a_list, int n_allele, params_t *params, int min_reads, int is_amplicon);


#endif // MODEL_H