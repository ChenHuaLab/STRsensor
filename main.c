/*************************************************************************
    > File Name: main.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2020年10月27日 星期一 22时22分22秒
    > Updated Time: 2021年01月08日 星期五 14时51分18秒
 ************************************************************************/

#include "share.h"
#include "utils.h"
#include "getseq.h"
#include "model.h"
#include "bamio.h"
#include "kmer.h"
#include "cigar.h"

#include <time.h>


static int locus_search(char *locus, char **loci, int n_loci)
{
    int index = -1;

    for (int idx=0; idx < n_loci; ++idx) {
        if (strcmp(locus, loci[idx]) == 0) return idx;
    }
    return index;
}


region_t *read_region(char *region_file, char *fasta_file, uint32_t index_kmer)
{
    region_t *region;
    fasta_t *fa_obj;

    err_calloc(region, 1, region_t);
    fa_obj = fa_init(fasta_file);

    char haps[8], buf[0xff];
    FILE *reg_fp;
    locus_t *l;
    int is_header = 1;

    err_open(reg_fp, region_file, "r");
    while (fgets(buf, 0xff, reg_fp)) {
        if (is_header) { is_header = 0; continue; }

        if (region->n == region->m) {
            region->m = region->m ? region->m + 4 : 4;
            region->locus = (locus_t *)realloc(region->locus, sizeof(locus_t) * region->m);
            region->loci = (char **)realloc(region->loci, sizeof(char*) * region->m);
        }

        l = &(region->locus[region->n]);
        sscanf(buf, "%s%s%d%d%d%d%s", 
                l->name, l->chr, &l->start, &l->end, &l->motif_len, &l->exclude_base, haps);
        l->is_haplotype = (!strcmp(haps, "Yes")) ? 1 : 0;
        l->index_kmer = index_kmer;
        l->fa_obj = fa_obj;
        region->loci[region->n++] = x_strcopy(l->name);
        
        index_t *index;

        for (int i=0; i < 2; ++i) {
            /* i:0 -> FLANK5; i:1 -> FLANK3 */
            index = flanks_index(l->chr, l->start, l->end, fa_obj, l->motif_len, index_kmer, i);
            if (i == FLANK5) l->index5 = index;
            else l->index3 = index;
        }
    } return region;
}


freq_t *read_freq(char *freq_file, char **loci, int n_loci)
{
    freq_t *freq;
    char c_str[0xf], f_str[0xf];

    err_calloc(freq, 1, freq_t);
    { /* pre-allocate memory for loci from STR region */
        freq->n = n_loci;
        err_calloc(freq->core, freq->n, cal_t);
    }
    for (int i=0; i < freq->n; ++i) strcpy(freq->core[i].name, loci[i]);
    if (!freq_file) return freq; /* frequency file is not given */

    gzFile fp;
    cal_t *f;
    kstring_t ks = {0, 0, NULL};
    int *field, field_n, idx;

    err_gzopen(fp, freq_file, "r");
    while (!kgetline(&ks, (kgets_func *)x_gzgets, fp)) {
        if (ks.s[0]=='\r' || ks.s[0]=='\n') continue; /* skip blank line */
        
        field = ksplit(&ks, 0, &field_n);
        idx = locus_search(ks.s+field[0], loci, n_loci);
        if (idx == -1) continue; /* only keep the freq occured in the region locus */

        f = &(freq->core[idx]);
        f->n = f->m = field_n - 1;
        kroundup32(f->m); /* 3->4; 5->8; 9->16; ...*/
        err_calloc(f->cal, f->m, cal1_t);

        for (int i=1; i < field_n; ++i) {
            sscanf(ks.s+field[i], "%[0-9.]:%[0-9.]", c_str, f_str); /* eg. 25.3:0.12345 */
            f->cal[i-1].allele = atof(c_str); f->cal[i-1].nf = atof(f_str);
        }
        free(field); ks.l = 0;
    }
    return freq;
}


stutter_t *read_stutter(char *stutter_file, char **loci, int n_loci)
{
    stutter_t *stutter;

    err_calloc(stutter, 1, stutter_t);
    { /* pre-allocate memory for loci from STR region */
        stutter->n = n_loci;
        err_calloc(stutter->core, stutter->n, stutter1_t);
    }
    for (int i=0; i < stutter->n; ++i) strcpy(stutter->core[i].name, loci[i]);
    if (!stutter_file) return stutter; /* frequency file is not given */

    gzFile fp;
    stutter1_t *s;
    kstring_t ks = {0, 0, NULL};
    int *field, field_n, idx;

    err_gzopen(fp, stutter_file, "r");
    while (!kgetline(&ks, (kgets_func *)x_gzgets, fp)) {
        if (ks.s[0]=='\r' || ks.s[0]=='\n') continue; /* skip blank line */
        
        field = ksplit(&ks, 0, &field_n);
        idx = locus_search(ks.s+field[0], loci, n_loci);
        if (idx == -1) continue; /* only keep the stutter occured in the locus file */

        s = &(stutter->core[idx]);
        s->num = 5; /* indicate the stutter file contains this locus */
        for (int i=1; i < field_n; ++i) sscanf(ks.s+field[i], "%f", &s->s[i-1]);
        free(field); ks.l = 0;
    }
    return stutter;
}


sample_t *read_sample(char *sample_list)
{
    sample_t *sample;
    FILE *fp;
    char buf[0x1000];
    int cur_idx = 0;

    err_calloc(sample, 1, sample_t);
    sample->n = x_linenum(sample_list);

    if (sample->n == 0) {
        fprintf(stderr, "[main/read_sample] Error: NO sample [*.bam] is detected in %s\n", sample_list);
        return NULL;
    }
    err_calloc(sample->name, sample->n, char*);
    err_calloc(sample->path, sample->n, char*);
    err_open(fp, sample_list, "r");

    while (fgets(buf, 0x1000, fp)) {
        if (buf[0]=='\r' || buf[0]=='\n') continue; /* skip blank line */
        sample->path[cur_idx] = x_strstrip(buf);
        sample->name[cur_idx] = x_strcopy(x_basename(sample->path[cur_idx]));
        cur_idx++;
    }
    fclose(fp); return sample;
}


void write_params(stutter_t *stutter, freq_t *freq, char *out_path)
{
    FILE *freq_fp, *stu_fp;
    char fname[0x200];

    /* write the frequency file */
    cal_t *f;
    sprintf(fname, "%s/Frequency.txt", out_path);
    err_open(freq_fp, fname, "w");

    for (int i=0; i < freq->n; ++i) {
        f = &freq->core[i]; 
        if (f->n == 0) continue; /* no frequency available for the locus */

        fprintf(freq_fp, "%s\t", f->name); /* locus name */
        for (int j=0; j < f->n; ++j)
            fprintf(freq_fp, "%.1f:%.6f\t", f->cal[j].allele, f->cal[j].nf);
        fprintf(freq_fp, "\n");
    }
    /* write the parameters */
    stutter1_t *s;
    sprintf(fname, "%s/Stutter.txt", out_path);
    err_open(stu_fp, fname, "w");

    for (int i=0; i < stutter->n; ++i) {
        s = &stutter->core[i]; 
        if (s->num == 0) continue; /* no stutter available for the locus */

        fprintf(stu_fp, "%s\t", s->name); /* locus name */
        for (int j=0; j < s->num; ++j) fprintf(stu_fp, "%.6f\t", s->s[j]);
        fprintf(stu_fp, "\n");
    }
    fclose(freq_fp); fclose(stu_fp);
}


int keep_read_check(bam1_t *aln, xcigar_t *cigar, int allow_dup)
{
    int match_bases = 0;
    int keep_read = 0;

    if (allow_dup == 0) { /* not allow duplicate (WGS)*/
        if (xbam_is_secondary(aln)) return (keep_read);
        if (xbam_is_supplementary(aln)) return (keep_read);
    }
    for (int i=0; i < cigar->n; ++i) {
        if (cigar->cigar_op[i] == 'H') return (keep_read); /* hard clip */
        if (cigar->cigar_op[i] == 'M') match_bases += cigar->cigar_len[i];
    }
    return keep_read + 1;
}


void locus_process_core(xbam_t *bam, person1_t *idv, locus_t *lobj, int frag, int allow_dup, int mis_match)
{
    int ret=0, keep=0;
    float allele;
    kstring_t ks = {0, 0, NULL};
    xcigar_t cigar = {0, 0, NULL, NULL};

    while (ret = xbam_next(bam->bam_hd, bam->iter[frag], bam->aln) >= 0) {
        idv->total_reads++;
        xbam_cigar(bam->aln, &cigar); /* parse the cigar value */
        keep = keep_read_check(bam->aln, &cigar, allow_dup);
        if (keep == 0) continue; /* not keep the read */

        idv->valid_reads++;
        xbam_seq(bam->aln, &ks); /* parse the read sequence */

        allele = cigar_allele(ks.s, &cigar, bam->aln->core.pos+1, lobj, mis_match);
        if (allele != -1.0) { /* extract allele by cigar algorithm */
            idv->span_reads++; allele_push(idv, allele);
            continue;
        }

        allele = kmer_allele(ks.s, lobj, mis_match);
        if (allele != -1.0) { /* extract allele by kmer algorithm */
            idv->span_reads++; allele_push(idv, allele);
        }
    }
    free(ks.s); 
    xbam_cigar_destroy(&cigar);
}


person_t *locus_process(sample_t *sample, char *loci_name, region_t *region, int allow_dup, int mis_match)
{
    person_t *person;

    err_calloc(person, 1, person_t);
    person->n_idv = sample->n;
    strcpy(person->locus, loci_name);
    err_calloc(person->idv, person->n_idv, person1_t);

    char name[0x40]; /* used for 2th-fragment locus */
    locus_t *locus[2]; /* locus[2] is used in 2th-fragment locus */
    int n_locus=0, idx, nl;

    nl = strlen(loci_name);
    idx = locus_search(loci_name, region->loci, region->n);
    if (idx == -1) return NULL; /* no possible to be here */ 
    locus[n_locus++] = &region->locus[idx];

    if (loci_name[nl-1] == 'a') { /* 2-fragment locus: DYS385a etc. */
        strcpy(name, loci_name); name[nl-1] = 'b';
        idx = locus_search(name, region->loci, region->n);
        if (idx == -1) return NULL; /* no possible to be here */ 
        locus[n_locus++] = &region->locus[idx];
    }
    xbam_t *bam;
    person1_t *idv;

    for (int i=0; i < sample->n; ++i) {
        bam = xbam_init(sample->path[i]); idv = &person->idv[i];
        if (bam == NULL) continue; /* something wrong with the bam file */
        idv->sample = x_strcopy(sample->name[i]);

        for (int j=0; j < n_locus; ++j) { /* get the iter for each fragment */
            int ret = xbam_fetch(locus[j]->chr, locus[j]->start, locus[j]->end, bam);
            if (ret == -1) break; /* failed to fetch the location */
            locus_process_core(bam, idv, locus[j], j, allow_dup, mis_match);
        } xbam_destroy(bam);
    }
    return person;
}


void gen_allele_string(float *a_list, int n_allele, char *a_string, cal_t *cal)
{
    char buf[0x40];
    a_string[0] = '\0'; /* makes strcat available */

    /* calculate the allele distribute. eg. 11.0:549 */
    for (int i=0; i < n_allele; ++i) allele_count(cal, a_list[i]);
    qsort(cal->cal, cal->n, sizeof(cal1_t), _allele_descend);

    for (int i=0; i < cal->n; ++i) {
        if (i != cal->n -1) sprintf(buf, "%.1f:%d,", cal->cal[i].allele, (int)cal->cal[i].nf);
        else sprintf(buf, "%.1f:%d", cal->cal[i].allele, (int)cal->cal[i].nf);
        strcat(a_string, buf);
    }
    cal->n = 0; /* maintain the memory */
}


void write_locus_allele(arg_t *args, person_t *idvs, params_t *params, locus_t *locus)
{
    char buf[0x200];
    FILE *file_fp;

    if (idvs->locus[strlen(idvs->locus)-1] == 'a') strcat(idvs->locus, "b");
    sprintf(buf, "%s/%s.txt", args->outpath, idvs->locus);
    err_open(file_fp, buf, "w");

    /* basic information of the locus */
    char haps_mark[0x2][0x8] = {"No", "Yes"};
    fprintf(file_fp, "#Locus\t%s\n", idvs->locus);
    fprintf(file_fp, "#Position\t%s:%d-%d\n", locus->chr, locus->start, locus->end);
    fprintf(file_fp, "#MotifLength\t%d\n", locus->motif_len);
    fprintf(file_fp, "#IsHaplotype\t%s\n", haps_mark[locus->is_haplotype]);
    fprintf(file_fp, "#MinimalReads\t%d\n", args->min_reads);
    fprintf(file_fp, "#MinimalProbability\t%.2f\n", args->min_prob);

    /* evaluated candidate alleles and stutter parameters */
    fprintf(file_fp, "#CandidateAllele\t");
    cal_t *c = params->f_loci;
    qsort(c->cal, c->n, sizeof(cal1_t), _candidate_ascend);
    for (int i=0; i < c->n; ++i) {
        if (i != c->n -1) fprintf(file_fp, "%.1f:%.6f,", c->cal[i].allele, c->cal[i].nf);
        else fprintf(file_fp, "%.1f:%.6f\n", c->cal[i].allele, c->cal[i].nf);
    }
    stutter1_t *s = params->s_loci;
    fprintf(file_fp, "#StutterParams\td_2:%.6f,d_1:%.6f,", s->s[3], s->s[4]);
    fprintf(file_fp, "n:%.6f,u_1:%.6f,u_2:%.6f\n", s->s[0], s->s[1], s->s[2]);

    /* headers for output information and detected alleles */
    cal_t cal = {0};
    person1_t *idv;
    fprintf(file_fp, "#Sample\tAllele\tProbability\tTotalReads\tValidReads\tSpanReads\tAlleleList\n");
    for (int i=0; i < idvs->n_idv; ++i) {
        idv = &idvs->idv[i];
        gen_allele_string(idv->allele_list, idv->n, buf, &cal);
        fprintf(file_fp, "%s\t%s\t%.6f\t", idv->sample, idv->allele, idv->allele_prob);
        fprintf(file_fp, "%d\t%d\t%d\t%s\n", idv->total_reads, idv->valid_reads, idv->span_reads, buf);
    }
    free(cal.cal); fclose(file_fp); 
}


void locus_allele_summary(arg_t *args, person_t *idvs, params_t *params, locus_t *locus)
{
    person1_t *idv;
    allele1_t al;

    for (int i=0; i < idvs->n_idv; ++i) {
        idv = &idvs->idv[i];
        al = get_allele(idv->allele_list, idv->n, params, args->min_reads, args->allow_dup);
        idv->allele_prob = al.prob;

        if (params->is_haplotype) { /* haplotype (X, Y)*/
            if (al.prob > args->min_prob) sprintf(idv->allele, "%.1f", al.allele1);
            else idv->allele[0] = '.';
        }
        else { /* autosome chromosome */
            float x_tmp;
            if (al.allele1 > al.allele2) {
                x_tmp = al.allele2; al.allele2 = al.allele1; al.allele1 = x_tmp;
            }
            if (al.prob > args->min_prob) sprintf(idv->allele, "%.1f/%.1f", al.allele1, al.allele2);
            else strcpy(idv->allele, "./.");
        }
    } write_locus_allele(args, idvs, params, locus);
}


int debug_main(int argc, char **argv)
{
    char *reg_file = "../Materials/Locus.txt";
    // char *freq_file = "../Materials/Frequency.txt";
    char *freq_file = NULL;
    // char *stu_file = "../Materials/Stutter.txt";
    char *stu_file = NULL;
    char *sample_list = "../Materials/bamlist.txt";
    char *fa_file = "/Users/xlzh/AppData/UCSC/Fasta/hg19.fa";

    arg_t args;
    args.allow_dup = 1;
    args.kmer = 12;
    args.mis_match = 2; args.min_reads = 10; args.min_prob = 0.0;
    args.outpath = x_strcopy("/Users/xlzh/Downloads/CodeTest/Test");

    region_t *region = read_region(reg_file, fa_file, args.kmer);
    freq_t *freq = read_freq(freq_file, region->loci, region->n);
    stutter_t *stutter = read_stutter(stu_file, region->loci, region->n);
    sample_t *sample = read_sample(sample_list);
    {
        char *name = "D2S1338";
        locus_t *locus;
        int idx = locus_search(name, region->loci, region->n);
        locus = &region->locus[idx];

        person_t *idvs = locus_process(sample, name, region, args.allow_dup, args.mis_match);
        cal_t *f_loci = freq_locus_search(freq, name);
        stutter1_t *s_loci = stu_locus_search(stutter, name);
        params_t *params = locus_model(idvs->idv, idvs->n_idv, f_loci, s_loci, args.min_reads, locus->is_haplotype);
        
        locus_allele_summary(&args, idvs, params, locus);
        person_destroy(idvs);
    }
    write_params(stutter, freq, args.outpath);
}


int main(int argc, char **argv)
{
    time_t start, end;
    arg_t *args;
    region_t *region; /* locus information and index of flanks kmers */
    freq_t *freq; /* frequency for each STR locus */
    stutter_t *stutter; /* stutter probability for each STR locus */
    sample_t *sample; /* sample basename and path */

    time(&start);
    args = args_parse(argc, argv);
    omp_set_num_threads(args->threads);

    region = read_region(args->region, args->fasta, args->kmer);
    freq = read_freq(args->frequency, region->loci, region->n);
    stutter = read_stutter(args->stutter, region->loci, region->n);
    sample = read_sample(args->infile);

    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i < region->n; ++i) {
        locus_t *locus; /* structure for the given locus */
        person_t *idvs; /* allele list detected for each individual */
        cal_t *f_loci; /* frequency pointer to the str locus */
        stutter1_t *s_loci; /* stutter pointer to the STR locus */
        params_t *params; /* parameters (freq & stutter) of the STR locus */
        char *l_name; /* locus name of the STR. eg. DYS390 */
        
        l_name = region->loci[i];
        if (l_name[strlen(l_name)-1] == 'b') continue; /* eg. DYS385b */
        int idx = locus_search(l_name, region->loci, region->n);
        if (idx == -1) {
            fprintf(stderr, "[main] Error: unexpected condition occured at %s!\n", l_name);
            continue;
        }
        locus = &region->locus[i];
        fprintf(stderr, "[%s] Get the allele of each read for each sample ...\n", l_name);
        idvs = locus_process(sample, l_name, region, args->allow_dup, args->mis_match);

        f_loci = freq_locus_search(freq, l_name);
        s_loci = stu_locus_search(stutter, l_name);
        fprintf(stderr, "[%s] Evaluate the candidate allele and stutter parameters ...\n", l_name);
        params = locus_model(idvs->idv, idvs->n_idv, f_loci, s_loci, args->min_reads, locus->is_haplotype);

        fprintf(stderr, "[%s] Calculate the most likely alleles for each sample ...\n", l_name);
        locus_allele_summary(args, idvs, params, locus); person_destroy(idvs);
        fprintf(stderr, "[%s] Done\n\n", l_name);
    }
    if (args->params_out)
        write_params(stutter, freq, args->outpath);
    
    time(&end);
    fprintf(stderr, "Total time consume: %.1f(s)\n\n", difftime(end, start));
}

