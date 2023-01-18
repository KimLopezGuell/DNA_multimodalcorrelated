#include "schna_ar_p.h"

void nrerror(char error_text[])
     /* SCHNAaP standard error handler */
{
    fprintf(stderr, "SCHNAaP run-time error: <%s>\n", error_text);
    exit(1);
}

long *lvector(long nl, long nh)
     /* allocate an long vector with subscript range v[nl..nh] */
{
    long *v;

    v = (long *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

char *cvector(long nl, long nh)
     /* allocate an char vector with subscript range v[nl..nh] */
{
    char *v;

    v = (char *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(char)));
    if (!v)
        nrerror("allocation failure in cvector()");
    return v - nl + NR_END;
}

double *dvector(long nl, long nh)
     /* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in dvector()");
    return v - nl + NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
     /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

long **lmatrix(long nrl, long nrh, long ncl, long nch)
     /* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    long **m;

    /* allocate pointers to rows */
    m = (long **) malloc((size_t) ((nrow + NR_END) * sizeof(long *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (long *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(long)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
     /* allocate an char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    char **m;

    /* allocate pointers to rows */
    m = (char **) malloc((size_t) ((nrow + NR_END) * sizeof(char *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(char)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

void free_lvector(long *v, long nl, long nh)
     /* free an long vector allocated with lvector() */
{
    free((FREE_ARG) (v + nl - NR_END));
}

void free_cvector(char *v, long nl, long nh)
     /* free an char vector allocated with cvector() */
{
    free((FREE_ARG) (v + nl - NR_END));
}

void free_dvector(double *v, long nl, long nh)
     /* free a double vector allocated with dvector() */
{
    free((FREE_ARG) (v + nl - NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
     /* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}

void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch)
     /* free an long matrix allocated by lmatrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
     /* free an char matrix allocated by cmatrix() */
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}
