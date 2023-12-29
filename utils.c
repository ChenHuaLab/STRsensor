/*************************************************************************
    > File Name: utils.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2020年11月04日 星期三 10时34分26秒
 ************************************************************************/

#include "utils.h"

/* open gzip file */
char *x_gzgets(char *str, int n, gzFile fp)
{
    return gzgets(fp, str, n);
}


/* copy a string and allocate necessary memory */
char *x_strcopy(const char *str)
{
    char *new_str;

    size_t str_len = strlen(str) + 1;
    new_str = (char *)malloc(str_len * sizeof(char));
    if (!new_str) {
        fprintf(stderr, "[err:kstrdup] failed to malloc memory!\n");
        return NULL;
    }
    memcpy(new_str, str, str_len);

    return new_str;
}


/* copy n elements from src and allocate necessary memory */
char *x_strncopy(char *src, size_t n)
{
    char *des;

    des = (char *)malloc((n+1) * sizeof(char));
    if (!des) {
        fprintf(stderr, "\n[err:%s] failed to malloc memory!\n", __func__);
        exit(-1);
    }
    strncpy(des, src, n); des[n] = '\0';

    return des;
}

/* strip the '\r\n' and '\n' at the end of the string */
char *x_strstrip(char *src_str)
{
    char *des_str;
    int l_str, i;

    l_str = strlen(src_str) + 1;
    des_str = (char *)malloc(l_str * sizeof(char));

    for (i=0; i < l_str; ++i) {
        if (src_str[i] == '\r' || src_str[i] == '\n') break;
        des_str[i] = src_str[i];
    }
    des_str[i] = '\0'; 
    return des_str;
}


/* get the basename of the given path */
char *x_basename(char *file_path)
{
    #ifdef _WIN32
      #define _PDELIM_ 92 // '\'
    #else
      #define _PDELIM_ 47 // '/'
    #endif

    char *start;

    start = strrchr(file_path, _PDELIM_);
    if (start != NULL) return start + 1;

    return file_path;
}


/* get number of lines from given text file */
int x_linenum(char *file_name)
{
    FILE *fp;
    int ln = 0;
    char buf[0x1000];

    err_open(fp, file_name, "r");
    while (fgets(buf, 0x1000, fp)) {
        /* skip the line start with '#' and empty line */
        if (buf[0]!='#' && buf[0]!='\r' && buf[0]!='\n') ln++;
    }
    fclose(fp); return ln;
}