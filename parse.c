/*************************************************************************
    > File Name: parse.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2020年10月27日 星期一 16时18分30秒
    > Updated Time: 2021年01月08日 星期五 14时52分18秒
 ************************************************************************/

#include "share.h"
#include "utils.h"


static void show_usage(void)
{
    char *usage =
        "Usage: STRsensor -i <bamlist.txt> -r <region.txt> -f <genome.fa> -o <output_dir> [options]\n"
        "Options:\n"
        "       -h|--help             print help infomation\n\n"
        "[Required]\n"
        "       -i|--infile     FILE  bam file list, one sample per line [.txt]\n" 
        "       -r|--region     FILE  region file of the STR locus [.txt]\n"
        "       -f|--fasta      FILE  reference genome (must be same as the mapping process used) [.fa]\n"
        "       -o|--outpath    PATH  the path for output file\n\n"
        "[Optional input parameters]\n"
        "       -s|--stutter    FILE  the stutter parameters file for each STR locus [.txt]\n"
        "       -q|--frequency  FILE  the frequency file for each STR locus on each potential allele [.txt]\n\n"
        "[Optional output parameters]\n"
        "       -p|--params_out       output parameters learned by given samples (including\n"
        "                             stutter and allele frequency)\n\n"
        "[Optional filter rules]\n"
        "       -c|--min_prob   FLOAT the minimal probability allowed to call an allele [0.00]\n"
        "       -n|--min_reads  INT   the minimum reads needed to genotype a STR locus for an individual [10]\n"
        "       -m|--mis_match  INT   the maximum mismatch(mismatch/insert/delete) bases allowed in both\n"
        "                             3' and 5' flanking sequence [2]\n\n"
        "[Optional]\n"
        "       -a|--allow_dup        duplication is allowed in target and amplicon sequencing, but not in WGS\n"
        "       -t|--threads    INT   the number of threads used [1]\n\n";

    fprintf(stderr, "Program: STRsensor (v%s)\n", __STRSENSOR_VERSION__);
    fprintf(stderr, "CreateDate: %s\n", __STRSENSOR_CREATE_DATE__);
    fprintf(stderr, "UpdateDate: %s\n", __STRSENSOR_UPDATE_DATE__);
    fprintf(stderr, "Author: XiaolongZhang (xiaolongzhang2015@163.com)\n\n");
    fprintf(stderr, "%s", usage);
    exit(-1);
}


static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "infile", required_argument, NULL, 'i' },
    { "region", required_argument, NULL, 'r' },
    { "fasta", required_argument, NULL, 'f' },
    { "outpath", required_argument, NULL, 'o' },
    { "stutter", optional_argument, NULL, 's' },
    { "frequency", optional_argument, NULL, 'q' },
    { "allow_dup", no_argument, NULL, 'a' },
    { "params_out", no_argument, NULL, 'p' },
    { "threads", optional_argument, NULL, 't' },
    { "min_prob", optional_argument, NULL, 'c' },
    { "mis_match", optional_argument, NULL, 'm' },
    { "min_reads", optional_argument, NULL, 'n' },
    { NULL, 0, NULL, 0 }
};


arg_t *args_parse(int argc, char **argv)
{
    int opt = 0, opterr = 0;
    arg_t *args;

    err_calloc(args, 1, arg_t);

    /* set the default parameters*/
    args->threads = 1; args->kmer = 8;
    args->mis_match = 2; args->min_reads = 10; args->min_prob = 0.0;

    if (argc < 5) show_usage();
    while ( (opt = getopt_long(argc, argv, "i:r:f:o:s:q:t:c:m:n:hap", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': args->help = 1; break;
            case 'i': args->infile = x_strcopy(optarg); break;
            case 'r': args->region = x_strcopy(optarg); break;
            case 'f': args->fasta = x_strcopy(optarg); break;
            case 'o': args->outpath = x_strcopy(optarg); break;
            case 's': args->stutter = x_strcopy(optarg); break;
            case 'q': args->frequency = x_strcopy(optarg); break;
            case 'a': args->allow_dup = 1; args->kmer = 12; break;
            case 'p': args->params_out = 1; break;
            case 't': args->threads = atoi(optarg); break;
            case 'c': args->min_prob = atof(optarg); break;
            case 'm': args->mis_match = atoi(optarg); break;
            case 'n': args->min_reads = atoi(optarg); break;
            case '?': fprintf(stderr, "[Err:parse_opt] Option error occour!\n\n"); args->help = 1;
        }
    }

    if (args->threads > 1) {
        int thr = omp_get_num_procs();
        args->threads = args->threads > thr ? thr : args->threads;
    }
    /* required parameters */
    if (args->help||!args->infile||!args->region||!args->fasta||!args->outpath)
        show_usage();

    return args;
}


