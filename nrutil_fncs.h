
/* nrutil.c */
void nrerror(char error_text[]);
long *lvector(long nl, long nh);
char *cvector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
long **lmatrix(long nrl, long nrh, long ncl, long nch);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_lvector(long *v, long nl, long nh);
void free_cvector(char *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
