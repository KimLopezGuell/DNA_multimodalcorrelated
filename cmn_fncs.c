/*<<<
 * ------------------------------------------------------------------------
 * This file is part of the SCHNAarP software package for the analysis and
 * rebuilding of double helical nucleic acid structures. Detailed description
 * of the underlying algorithms can be found in the following three papers:
 * 
 *    (1) Xiang-Jun Lu, M. A. El Hassan and C. A. Hunter (1997), "Structure
 *        and Conformation of Helical Nucleic Acids: Analysis Program
 *        (SCHNAaP)", J. Mol. Biol. 273, 668--680.
 * 
 *    (2) Xiang-Jun Lu, M. A. El Hassan and C. A. Hunter (1997), "Structure
 *        and Conformation of Helical Nucleic Acids: Rebuilding Program
 *        (SCHNArP)", J. Mol. Biol. 273, 681--691.
 * 
 *    (3) M. A. El Hassan and C. R. Calladine (1995), "The Assessment of the
 *        Geometry of Dinucleotide Steps in Double-Helical DNA; a New Local
 *        Calculation Scheme", J. Mol. Biol. 251, 648--664.
 * 
 * Copyright (C) 1996-2002 Xiang-Jun Lu. Permission to use, copy, and modify
 * this software is hereby granted to all academic and not-for-profit
 * institutions without fee, provided that the above papers (or part of) are
 * properly cited in related publications.  However, the right to use this
 * software in conjunction with for profit activities, and the right to
 * distribute the software or modified or extended versions thereof for profit
 * are *NOT* granted except by prior arrangement and written consent of the
 * copyright holders.
 * 
 * The software is provided "AS-IS" and without warranty of any kind, express,
 * implied or otherwise, including without limitation, any warranty of
 * merchantability or fitness for a particular purpose.
 * 
 * The software was written by Xiang-Jun Lu (3dna.lu@gmail.com). It started
 * off his PhD thesis work at Dr. Chris Hunter Laboratory at University of
 * Sheffield, U.K. The source code should be valuable for understanding the
 * elegant CEHS algorithm, and many other aspects of nucleic acid structures.
 * For real world applications, however, you are strongly recommended to
 * switch to 3DNA (http://3dna.rutgers.edu).
 * 
 * Revision history -- please see file 'Change.log'
 * 
 * $LastChangedDate: 2009-06-01 18:34:25 -0400 (Mon, 01 Jun 2009) $
 * $LastChangedRevision: 1556 $
 * $LastChangedBy: xiangjun $
 * ------------------------------------------------------------------------
>>>*/

#include "schna_ar_p.h"

/* >>>adapt the following 4 functions from  NR in C (2nd Edition) <<< */

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void jacobi(double **a, long n, double *d, double **v, long *nrot)
{
    long j, iq, ip, i;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

/* BEGIN INSERT jacobi.c HERE */

/* the first two lines would be:
    b = dvector(1, n);
    z = dvector(1, n);
*/

/* (float) ----> (double)
    if (i > 4 && (double) (fabs(d[ip]) + g) == (double) fabs(d[ip])
                  ^^^^^^
*/

/* end with
    nrerror("Too many iterations.");
*/

	b=dvector(1,n);
	z=dvector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z,1,n);
			free_dvector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");

/* END INSERT jacobi.c HERE */
}

#undef ROTATE

void eigsrt(double *d, double **v, long n)
{
    long k, j, i;
    double p;

/* BEGIN INSERT eigsrt.c HERE */
for (i=1;i<n;i++) { 
		p=d[k=i]; 
		for (j=i+1;j<=n;j++) 
			if (d[j] >= p) p=d[k=j]; 
		if (k != i) { 
			d[k]=d[i]; 
			d[i]=p; 
			for (j=1;j<=n;j++) { 
				p=v[j][i]; 
				v[j][i]=v[j][k]; 
				v[j][k]=p; 
			} 
		} 
	} 
/* END INSERT eigsrt.c HERE */
}

#define TINY 1.0e-20

void dludcmp(double **a, long n, long *indx, double *d)
{
    long i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;

/* BEGIN INSERT ludcmp.c HERE */

/* start with
    vv = dvector(1, n);
    *d = 1.0;
*/

/* ends with
    free_dvector(vv, 1, n);
*/



        vv=dvector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_dvector(vv,1,n);

/* END INSERT ludcmp.c HERE */

}

#undef TINY

void dlubksb(double **a, long n, long *indx, double *b)
{
    long i, ii = 0, ip, j;
    double sum;

/* BEGIN INSERT lubksb.c HERE */
        for (i=1;i<=n;i++) {
                ip=indx[i];
                sum=b[ip];
                b[ip]=b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii=i;
                b[i]=sum;
        }
        for (i=n;i>=1;i--) {
                sum=b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i]=sum/a[i][i];
        }

/* END INSERT lubksb.c HERE */
}

/* >>> END of 4 functions from NR in C (2nd Edition) <<< */

double dist_ab(double a, double b)
{
    return (sqrt(a * a + b * b));
}

void zero_dmatrix(double **a, long nr, long nc)
{
    long i, j;

    for (i = 1; i <= nr; i++) {
        for (j = 1; j <= nc; j++)
            a[i][j] = 0.0;
    }
}

void one_dmatrix(double **a, long nr, long nc)
{
    long i, j;

    for (i = 1; i <= nr; i++) {
        for (j = 1; j <= nc; j++)
            a[i][j] = 1.0;
    }
}

void zero_dvector(double *a, long n)
{
    long i;

    for (i = 1; i <= n; i++)
        a[i] = 0.0;
}

void zero_lvector(long *a, long n)
{
    long i;

    for (i = 1; i <= n; i++)
        a[i] = 0;
}

void identity_dmatrix(double **a, long n)
{
    long i, j;

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++)
            a[i][j] = 0.0;
        a[i][i] = 1.0;
    }
}

double deg2rad(double theta)
{
    return (theta * PI / 180.0);
}

double rad2deg(double theta)
{
    return (theta * 180.0 / PI);
}

void rotx(double **rotmat, double theta)
     /* "theta" MUST be in degree */
{
    double c, s, ang;

    ang = deg2rad(theta);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = 1.0;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = c;
    rotmat[2][3] = -s;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = s;
    rotmat[3][3] = c;
}

void roty(double **rotmat, double theta)
     /* "theta" MUST be in degree */
{
    double c, s, ang;

    ang = deg2rad(theta);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = s;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = 1.0;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = -s;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = c;
}

void rotz(double **rotmat, double theta)
     /* "theta" MUST be in degree */
{
    double c, s, ang;

    ang = deg2rad(theta);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = -s;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = s;
    rotmat[2][2] = c;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = 1.0;
}

void mul_dmatrix(double **o, double **a, double **b,
                 long nr_a, long nc_a, long nr_b, long nc_b)
{
    double tmp;
    long i, j, k;

    if (nc_a != nr_b)
        nrerror("Matrices a and b do NOT conform!");

    for (i = 1; i <= nr_a; i++) {
        for (j = 1; j <= nc_b; j++) {
            tmp = 0.0;
            for (k = 1; k <= nc_a; k++)
                tmp += a[i][k] * b[k][j];
            o[i][j] = tmp;
        }
    }
}

void tr_dmatrix(double **o, double **a, long nr, long nc)
     /* Matrix transpose */
{
    long i, j;

    for (i = 1; i <= nc; i++) {
        for (j = 1; j <= nr; j++)
            o[i][j] = a[j][i];
    }
}

void cov_dmatrix(double **cov_mtx, double **xyz, long nr, long nc)
     /* Get the covariance matrix of XYZ */
{
    long i, j;
    double **xyzT, **o, **oT, **tmp1, **tmp2, **tmp3, **tmp4;

    xyzT = dmatrix(1, nc, 1, nr);
    o = dmatrix(1, nr, 1, 1);
    oT = dmatrix(1, 1, 1, nr);
    tmp1 = dmatrix(1, nc, 1, nc);
    tmp2 = dmatrix(1, nc, 1, 1);
    tmp3 = dmatrix(1, 1, 1, nc);
    tmp4 = dmatrix(1, nc, 1, nc);

    one_dmatrix(o, nr, 1);
    tr_dmatrix(oT, o, nr, 1);
    tr_dmatrix(xyzT, xyz, nr, nc);

    mul_dmatrix(tmp1, xyzT, xyz, nc, nr, nr, nc);
    mul_dmatrix(tmp2, xyzT, o, nc, nr, nr, 1);
    mul_dmatrix(tmp3, oT, xyz, 1, nr, nr, nc);
    mul_dmatrix(tmp4, tmp2, tmp3, nc, 1, 1, nc);

    for (i = 1; i <= nc; i++) {
        for (j = 1; j <= nc; j++)
            tmp4[i][j] = tmp4[i][j] / nr;
    }

    for (i = 1; i <= nc; i++) {
        for (j = 1; j <= nc; j++)
            cov_mtx[i][j] = (tmp1[i][j] - tmp4[i][j]) / (nr - 1);
    }

    free_dmatrix(o, 1, nr, 1, 1);
    free_dmatrix(oT, 1, 1, 1, nr);
    free_dmatrix(tmp2, 1, nc, 1, 1);
    free_dmatrix(tmp3, 1, 1, 1, nc);
    free_dmatrix(xyzT, 1, nc, 1, nr);
    free_dmatrix(tmp1, 1, nc, 1, nc);
    free_dmatrix(tmp4, 1, nc, 1, nc);
}

void mean_dmatrix(double *o, double **a, long nr, long nc)
     /* Columnwise mean of matrix a */
{
    long i, j;
    double dsum;

    for (i = 1; i <= nc; i++) {
        dsum = 0;
        for (j = 1; j <= nr; j++)
            dsum += a[j][i];
        o[i] = dsum / nr;
    }
}

void std_dmatrix(double *o, double **a, long nr, long nc)
     /* Columnwise std of matrix a */
{
    long i, j;
    double tmp, dsum, *ave_a;

    if (nr < 1)
        nrerror("Number of samples < 1");

    ave_a = dvector(1, nc);

    mean_dmatrix(ave_a, a, nr, nc);

    for (i = 1; i <= nc; i++) {
        dsum = 0;
        for (j = 1; j <= nr; j++) {
            tmp = a[j][i] - ave_a[i];
            dsum += DSQR(tmp);
        }
        o[i] = sqrt(dsum / (nr - 1));
    }
    free_dvector(ave_a, 1, nc);
}

double mean_dvector(double *a, long n)
{
    long i;
    double dsum = 0;

    for (i = 1; i <= n; i++)
        dsum += a[i];
    return (dsum / n);
}

double std_dvector(double *a, long n)
{
    long i;
    double mean_a, dsum = 0, tmp;

    if (n < 1)
        nrerror("Number of samples < 1");

    mean_a = mean_dvector(a, n);

    for (i = 1; i <= n; i++) {
        tmp = a[i] - mean_a;
        dsum += DSQR(tmp);
    }

    return (sqrt(dsum / (n - 1)));
}

double dot_dvector(double *a, double *b, long n)
{
    long i;
    double dsum = 0;

    for (i = 1; i <= n; i++)
        dsum += a[i] * b[i];

    return dsum;
}

void cross_dvector(double *o, double *a, double *b)
{
    o[1] = a[2] * b[3] - a[3] * b[2];
    o[2] = a[3] * b[1] - a[1] * b[3];
    o[3] = a[1] * b[2] - a[2] * b[1];
}

double mag_dvector(double *a, long n)
     /* Get the magnitude of a vector */
{
    return (sqrt(dot_dvector(a, a, n)));
}

void norm_dvector(double *o, double *a, long n)
     /* Normalize a vector */
{
    long i;
    double mag;

    mag = mag_dvector(a, n);
    if (mag <= EPS)
        zero_dvector(o, n);  /* a must be a zero vector */
    else {
        for (i = 1; i <= n; i++)
            o[i] = a[i] / mag;
    }
}

void lsplane(double *vec_normal, double *org_dist, double *pnts_dist,
             double **xyz, long nr)
     /* Get the LS plane */
{
    long i, j, nrot, nc = 3;
    double *d, **v, **cov_mtx;

    d = dvector(1, nc);
    v = dmatrix(1, nc, 1, nc);
    cov_mtx = dmatrix(1, nc, 1, nc);

    cov_dmatrix(cov_mtx, xyz, nr, nc);
    jacobi(cov_mtx, nc, d, v, &nrot);
    eigsrt(d, v, nc);

    for (i = 1; i <= nc; i++)
        vec_normal[i] = v[i][nc];
    if ((d[1] - d[nc]) <= EPS) {
        vec_normal[1] = 0;
        vec_normal[2] = 0;
        vec_normal[3] = 1;
    }

    if (vec_normal[nc] < 0) {
        for (i = 1; i <= nc; i++)
            vec_normal[i] = -vec_normal[i];
    }

    mean_dmatrix(d, xyz, nr, nc);
    *org_dist = dot_dvector(d, vec_normal, nc);

    for (i = 1; i <= nr; i++) {
        for (j = 1; j <= nc; j++)
            pnts_dist[i] += xyz[i][j] * vec_normal[j];
        pnts_dist[i] -= *org_dist;
    }

    free_dvector(d, 1, nc);
    free_dmatrix(v, 1, nc, 1, nc);
    free_dmatrix(cov_mtx, 1, nc, 1, nc);
}

void lsline(double *vec_line, double *xyzpoint, double *pnts_dist, double **xyz, long nr)
     /* Get the LS line */
{
    long i, j, nrot, nc = 3;
    double *d, **v, **cov_mtx, *o;

    d = dvector(1, nc);
    v = dmatrix(1, nc, 1, nc);
    cov_mtx = dmatrix(1, nc, 1, nc);
    o = dvector(1, nc);

    cov_dmatrix(cov_mtx, xyz, nr, nc);
    jacobi(cov_mtx, nc, d, v, &nrot);
    eigsrt(d, v, nc);

    for (i = 1; i <= nc; i++)
        vec_line[i] = v[i][1];
    if ((d[1] - d[nc]) <= EPS) {
        vec_line[1] = 0;
        vec_line[2] = 0;
        vec_line[3] = 1;
    }

    if (vec_line[nc] < 0) {
        for (i = 1; i <= nc; i++)
            vec_line[i] = -vec_line[i];
    }

    mean_dmatrix(xyzpoint, xyz, nr, nc);

    for (i = 1; i <= nr; i++) {
        for (j = 1; j <= nc; j++)
            d[j] = xyz[i][j] - xyzpoint[j];
        cross_dvector(o, d, vec_line);
        pnts_dist[i] = mag_dvector(o, nc);
    }

    free_dvector(d, 1, nc);
    free_dmatrix(v, 1, nc, 1, nc);
    free_dmatrix(cov_mtx, 1, nc, 1, nc);
    free_dvector(o, 1, nc);
}

void alignz(double **rotmat, double *rot_axis)
{
    double mag, v12, c, s, *dircos, **Rz, **Rx;

    dircos = dvector(1, 3);
    Rz = dmatrix(1, 3, 1, 3);
    Rx = dmatrix(1, 3, 1, 3);

    mag = mag_dvector(rot_axis, 3);
    if (mag <= EPS) {
        identity_dmatrix(rotmat, 3);
        return;
    } else
        norm_dvector(dircos, rot_axis, 3);

    /* Rotate about z-axis to make <dircos> in yz-plane */
    v12 = dist_ab(dircos[1], dircos[2]);
    if (v12 <= EPS)
        identity_dmatrix(Rz, 3);
    else {
        c = dircos[2] / v12;
        s = dircos[1] / v12;
        Rz[1][1] = c;
        Rz[1][2] = -s;
        Rz[1][3] = 0.0;
        Rz[2][1] = s;
        Rz[2][2] = c;
        Rz[2][3] = 0.0;
        Rz[3][1] = 0.0;
        Rz[3][2] = 0.0;
        Rz[3][3] = 1.0;
    }

    /* Rotate about x-axis to make <dircos> align with +z-axis */
    Rx[1][1] = 1.0;
    Rx[1][2] = 0.0;
    Rx[1][3] = 0.0;
    Rx[2][1] = 0.0;
    Rx[2][2] = dircos[3];
    Rx[2][3] = -v12;
    Rx[3][1] = 0.0;
    Rx[3][2] = v12;
    Rx[3][3] = dircos[3];

    /* Final rotation matrix */
    mul_dmatrix(rotmat, Rx, Rz, 3, 3, 3, 3);

    free_dvector(dircos, 1, 3);
    free_dmatrix(Rz, 1, 3, 1, 3);
    free_dmatrix(Rx, 1, 3, 1, 3);
}

void arbrot(double **rotmat, double *rot_axis, double theta)
{
    double mag, v12, c, s, ang;
    double *dircos, **Rz, **Rx, **R, **tmp, **tmpT;

    ang = deg2rad(theta);

    dircos = dvector(1, 3);
    Rz = dmatrix(1, 3, 1, 3);
    Rx = dmatrix(1, 3, 1, 3);
    R = dmatrix(1, 3, 1, 3);
    tmp = dmatrix(1, 3, 1, 3);
    tmpT = dmatrix(1, 3, 1, 3);

    mag = mag_dvector(rot_axis, 3);
    if (mag <= EPS) {
        identity_dmatrix(rotmat, 3);
        return;
    } else
        norm_dvector(dircos, rot_axis, 3);

    /* Rotate about z-axis to make <dircos> in yz-plane */
    v12 = dist_ab(dircos[1], dircos[2]);
    if (v12 <= EPS)
        identity_dmatrix(Rz, 3);
    else {
        c = dircos[2] / v12;
        s = dircos[1] / v12;
        Rz[1][1] = c;
        Rz[1][2] = -s;
        Rz[1][3] = 0.0;
        Rz[2][1] = s;
        Rz[2][2] = c;
        Rz[2][3] = 0.0;
        Rz[3][1] = 0.0;
        Rz[3][2] = 0.0;
        Rz[3][3] = 1.0;
    }

    /* Rotate about x-axis to make <dircos> align with +z-axis */
    Rx[1][1] = 1.0;
    Rx[1][2] = 0.0;
    Rx[1][3] = 0.0;
    Rx[2][1] = 0.0;
    Rx[2][2] = dircos[3];
    Rx[2][3] = -v12;
    Rx[3][1] = 0.0;
    Rx[3][2] = v12;
    Rx[3][3] = dircos[3];

    /* Real rotation is about z-axis by <ang> */
    c = cos(ang);
    s = sin(ang);
    R[1][1] = c;
    R[1][2] = -s;
    R[1][3] = 0.0;
    R[2][1] = s;
    R[2][2] = c;
    R[2][3] = 0.0;
    R[3][1] = 0.0;
    R[3][2] = 0.0;
    R[3][3] = 1.0;

    /* Final rotation matrix */
    mul_dmatrix(tmp, Rx, Rz, 3, 3, 3, 3);
    mul_dmatrix(rotmat, R, tmp, 3, 3, 3, 3);
    tr_dmatrix(tmpT, Rx, 3, 3);
    mul_dmatrix(tmp, tmpT, rotmat, 3, 3, 3, 3);
    tr_dmatrix(tmpT, Rz, 3, 3);
    mul_dmatrix(rotmat, tmpT, tmp, 3, 3, 3, 3);

    free_dvector(dircos, 1, 3);
    free_dmatrix(Rz, 1, 3, 1, 3);
    free_dmatrix(Rx, 1, 3, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
    free_dmatrix(tmp, 1, 3, 1, 3);
    free_dmatrix(tmpT, 1, 3, 1, 3);
}

double magang(double *veca, double *vecb)
     /* Get the magnitude of angle in degrees between two vectors */
{
    double *na, *nb, dval, ang;

    na = dvector(1, 3);
    nb = dvector(1, 3);

    /* Get the normalized vectors */
    norm_dvector(na, veca, 3);
    norm_dvector(nb, vecb, 3);

    /* Get the angle in degrees */
    dval = dot_dvector(na, nb, 3);

    /* Special case: due to round off error, abs(dval) could > 1.0!!! */
    if (dval > 1.0)
        dval = 1.0;
    if (dval < -1.0)
        dval = -1.0;

    ang = rad2deg(acos(dval));

    free_dvector(na, 1, 3);
    free_dvector(nb, 1, 3);

    return ang;
}

double vec_ang(double *veca, double *vecb, double *vec_ref)
     /* <veca> & <vecb> will be normalized in "magang"
        <vec_ref> is used for sign control and does NOT need to be normalized
        ang_deg> in the range of [-180, +180] */
{
    double ang, *aXb;

    aXb = dvector(1, 3);

    ang = magang(veca, vecb);
    cross_dvector(aXb, veca, vecb);
    if (dot_dvector(aXb, vec_ref, 3) < 0)
        ang = -ang;

    free_dvector(aXb, 1, 3);

    return ang;
}

void bpfunc(double **orien, double *pos, double *param)
     /* The order of the 6 base-pair parameters should be:
        %          x       y        z       x        y         z
        %        SHEAR  STRETCH  STAGGER  BUCKLE PROPELLER  OPENING
      */
{
    long i, j;
    double propeller, opening, buckle, stretch, stagger, shear;
    double buckleopening, phi, dsum;
    double *a, *b, *r, **tmp1, **tmp2, **mbp_orien;

    a = dvector(1, 3);
    b = dvector(1, 3);
    r = dvector(1, 3);
    tmp1 = dmatrix(1, 3, 1, 3);
    tmp2 = dmatrix(1, 3, 1, 3);
    mbp_orien = dmatrix(1, 3, 1, 3);

    shear = param[1];
    stretch = param[2];
    stagger = param[3];
    buckle = param[4];
    propeller = param[5];
    opening = param[6];

    buckleopening = dist_ab(opening, buckle);

    a[1] = buckle;
    a[2] = 0;
    a[3] = opening;
    b[1] = 1.0;
    b[2] = 0.0;
    b[3] = 0.0;
    r[1] = 0.0;
    r[2] = 1.0;
    r[3] = 0.0;

    phi = vec_ang(a, b, r);

    roty(orien, 0.5 * propeller - phi);
    rotx(tmp1, buckleopening);
    mul_dmatrix(tmp2, orien, tmp1, 3, 3, 3, 3);
    roty(tmp1, 0.5 * propeller + phi);
    mul_dmatrix(orien, tmp2, tmp1, 3, 3, 3, 3);

    roty(mbp_orien, 0.5 * propeller - phi);
    rotx(tmp1, 0.5 * buckleopening);
    mul_dmatrix(tmp2, mbp_orien, tmp1, 3, 3, 3, 3);
    roty(tmp1, phi);
    mul_dmatrix(mbp_orien, tmp2, tmp1, 3, 3, 3, 3);

    a[1] = shear;
    a[2] = stretch;
    a[3] = stagger;
    tr_dmatrix(tmp1, mbp_orien, 3, 3);

    for (i = 1; i <= 3; i++) {
        dsum = 0;
        for (j = 1; j <= 3; j++)
            dsum += a[j] * tmp1[j][i];
        pos[i] = dsum;
    }

    free_dvector(a, 1, 3);
    free_dvector(b, 1, 3);
    free_dvector(r, 1, 3);
    free_dmatrix(tmp1, 1, 3, 1, 3);
    free_dmatrix(tmp2, 1, 3, 1, 3);
    free_dmatrix(mbp_orien, 1, 3, 1, 3);
}

void stepfunc(double **orien, double **mst, double *pos, double *param)
     /* The order of the 6 step parameters should be:
        %          x     y    z       x    y    z
        %        SHIFT SLIDE RISE   TILT ROLL TWIST
      */
{
    long i, j;
    double twist, roll, tilt, rise, slide, shift;
    double rolltilt, phi, dsum;
    double *a, *b, *r, **tmp1, **tmp2;

    a = dvector(1, 3);
    b = dvector(1, 3);
    r = dvector(1, 3);
    tmp1 = dmatrix(1, 3, 1, 3);
    tmp2 = dmatrix(1, 3, 1, 3);

    shift = param[1];
    slide = param[2];
    rise = param[3];
    tilt = param[4];
    roll = param[5];
    twist = param[6];

    rolltilt = dist_ab(tilt, roll);

    a[1] = tilt;
    a[2] = roll;
    a[3] = 0.0;
    b[1] = 0.0;
    b[2] = 1.0;
    b[3] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    r[3] = 1.0;

    phi = vec_ang(a, b, r);

    rotz(orien, 0.5 * twist - phi);
    roty(tmp1, rolltilt);
    mul_dmatrix(tmp2, orien, tmp1, 3, 3, 3, 3);
    rotz(tmp1, 0.5 * twist + phi);
    mul_dmatrix(orien, tmp2, tmp1, 3, 3, 3, 3);

    rotz(mst, 0.5 * twist - phi);
    roty(tmp1, 0.5 * rolltilt);
    mul_dmatrix(tmp2, mst, tmp1, 3, 3, 3, 3);
    rotz(tmp1, phi);
    mul_dmatrix(mst, tmp2, tmp1, 3, 3, 3, 3);

    a[1] = shift;
    a[2] = slide;
    a[3] = rise;
    tr_dmatrix(tmp1, mst, 3, 3);

    for (i = 1; i <= 3; i++) {
        dsum = 0;
        for (j = 1; j <= 3; j++)
            dsum += a[j] * tmp1[j][i];
        pos[i] = dsum;
    }

    free_dvector(a, 1, 3);
    free_dvector(b, 1, 3);
    free_dvector(r, 1, 3);
    free_dmatrix(tmp1, 1, 3, 1, 3);
    free_dmatrix(tmp2, 1, 3, 1, 3);
}

void cirshift(long *o, long *a, long n, long m)
{
    long i, j, ipm;
    for (i = 1; i <= n; i++) {
        ipm = i + m;
        if (ipm <= 0)
            ipm = ipm + n;
        if (ipm <= n)
            o[ipm] = a[i];
        else {
            j = ipm % n;
            o[j] = a[i];
        }
    }
}

void vec_orth(double *o, double *veca, double *vecb)
{
    long n = 3;
    double *nb, *tmp1, *tmp2;

    nb = dvector(1, n);
    tmp1 = dvector(1, n);
    tmp2 = dvector(1, n);

    norm_dvector(nb, vecb, n);
    cross_dvector(tmp1, veca, nb);
    cross_dvector(tmp2, nb, tmp1);
    norm_dvector(o, tmp2, n);

    free_dvector(nb, 1, n);
    free_dvector(tmp1, 1, n);
    free_dvector(tmp2, 1, n);
}

double val_ang(double **abcxyz)
{
    long i;
    double *ba, *bc, ang;

    ba = dvector(1, 3);
    bc = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        ba[i] = abcxyz[1][i] - abcxyz[2][i];
        bc[i] = abcxyz[3][i] - abcxyz[2][i];
    }

    ang = magang(ba, bc);

    free_dvector(ba, 1, 3);
    free_dvector(bc, 1, 3);

    return ang;
}

double torsion(double **abcdxyz)
{
    long i;
    double *ab, *bc, *cd, *p_abc, *p_bcd, ang;

    ab = dvector(1, 3);
    bc = dvector(1, 3);
    cd = dvector(1, 3);
    p_abc = dvector(1, 3);
    p_bcd = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        ab[i] = abcdxyz[2][i] - abcdxyz[1][i];
        bc[i] = abcdxyz[3][i] - abcdxyz[2][i];
        cd[i] = abcdxyz[4][i] - abcdxyz[3][i];
    }

    cross_dvector(p_abc, ab, bc);
    cross_dvector(p_bcd, bc, cd);

    ang = vec_ang(p_abc, p_bcd, bc);

    free_dvector(ab, 1, 3);
    free_dvector(bc, 1, 3);
    free_dvector(cd, 1, 3);
    free_dvector(p_abc, 1, 3);
    free_dvector(p_bcd, 1, 3);

    return ang;
}

void i_diff(long *o, long *a, long n)
{
    long i;
    for (i = 1; i < n; i++)
        o[i] = a[i + 1] - a[i];
}

void logical_not(long *o, long *a, long n)
{
    long i;
    for (i = 1; i <= n; i++)
        if (a[i] == 0)
            o[i] = 1;
        else
            o[i] = 0;
}

void isort(long n, long *ra)
{
    long l, j, ir, i;
    long rra;

    l = (n >> 1) + 1;
    ir = n;
    for (;;) {
        if (l > 1)
            rra = ra[--l];
        else {
            rra = ra[ir];
            ra[ir] = ra[1];
            if (--ir == 1) {
                ra[1] = rra;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j + 1])
                ++j;
            if (rra < ra[j]) {
                ra[i] = ra[j];
                j += (i = j);
            } else
                j = ir + 1;
        }
        ra[i] = rra;
    }
}

void dsort(long n, double *ra)
{
    long l, j, ir, i;
    double rra;

    l = (n >> 1) + 1;
    ir = n;
    for (;;) {
        if (l > 1)
            rra = ra[--l];
        else {
            rra = ra[ir];
            ra[ir] = ra[1];
            if (--ir == 1) {
                ra[1] = rra;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j + 1])
                ++j;
            if (rra < ra[j]) {
                ra[i] = ra[j];
                j += (i = j);
            } else
                j = ir + 1;
        }
        ra[i] = rra;
    }
}

/*  Find inverse of matrix 'a' and return result as 'y' */
void dinverse(double **a, long n, double **y)
{
    double d, *col;
    long i, j, *indx;

    col = dvector(1, n);
    indx = lvector(1, n);

    dludcmp(a, n, indx, &d);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            col[i] = 0.0;
        col[j] = 1.0;
        dlubksb(a, n, indx, col);
        for (i = 1; i <= n; i++)
            y[i][j] = col[i];
    }

    free_lvector(indx, 1, n);
    free_dvector(col, 1, n);
}

/* o=axb: a--[1..m], b--[1..m, 1..n], o--[1..n] */
void dvec_x_dmtx(double *o, double *a, double **b, long na, long nr, long nc)
{
    long i, j;
    double dsum;

    if (na != nr)
        nrerror("Elements # in <a> does not match # of row in b");
    for (i = 1; i <= nc; i++) {
        dsum = 0.0;
        for (j = 1; j <= na; j++)
            dsum += a[j] * b[j][i];
        o[i] = dsum;
    }
}

void duplicate_dmtx(double **o, double **a, long nr, long nc)
     /* Make a copy of matrix b to a */
{
    long i, j;

    for (i = 1; i <= nr; i++)
        for (j = 1; j <= nc; j++)
            o[i][j] = a[i][j];
}

double max_dvector(double *a, long n)
     /* Find the maximum value of vector <a> */
{
    long i;
    double max_d;

    max_d = 1.0e-10;

    for (i = 1; i <= n; i++)
        max_d = DMAX(max_d, a[i]);

    return max_d;
}

double min_dvector(double *a, long n)
     /* Find the minimum value of vector <a> */
{
    long i;
    double min_d;

    min_d = 1.0e+10;

    for (i = 1; i <= n; i++)
        min_d = DMIN(min_d, a[i]);

    return min_d;
}

/* ====================== String Functions ====================== */

void upper(char *a, long n)
{
    long i;

    for (i = 0; i < n; i++)
        a[i] = toupper(a[i]);
}

void read_stdin_str(char *strval)
{
    char str[1024] = "";

    strcpy(strval, "");  /* initialize it */
    fgets(str, BUF512, stdin);
    sscanf(str, "%s", strval);
}

void comp_base(char *o, char *a, long n)
{
    long i, j;
    char *base1 = "ATCG", *compb = "TAGC";

    upper(a, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < 4; j++)
            if (a[i] == base1[j])
                break;
        if (j < 4)
            o[i] = compb[j];
        else
            nrerror("Illegal base exists!");
    }

    o[n] = '\0';
    a[n] = '\0';  /* Just to make sure */
}

FILE *open_file(char *filnam, char *filmod)
{
    FILE *f1;
    errno = 0;

    if ((f1 = fopen(filnam, filmod)) == NULL) {
        printf("Open file <%s> error: %s\n", filnam, strerror(errno));
        nrerror("Failed to open file ...");
    }

    return f1;
}

void get_sequence(char **bpseq, long *nbp)
{
    FILE *f1;
    long iseq, i, n_rep;
    char str[BUF512], tmp[BUF512];

    printf("\nInput the base sequence:\n");
    printf("1. From a data file (complete sequence)\n");
    printf("2. From the keyboard (enter only the repeating sequence)\n");
    printf("Your choice (1 or 2, Dft 2): ");
    read_stdin_str(str);
    iseq = atoi(str);

    if (iseq == 1) {
        printf("Name of the base sequence file: ");
        read_stdin_str(str);
        f1 = open_file(str, "r");

        *nbp = 0;
        while (fscanf(f1, "%s", str) != EOF) {
            *nbp += strlen(str);
            strcat(bpseq[1], str);
        }
        fclose(f1);
    } else {
        printf("\nRepeat unit (Dft polyA): ");
        read_stdin_str(str);
        if (strlen(str) == 0)
            strcpy(str, "A\0");

        printf("Number of repeats (Dft 10): ");
        read_stdin_str(tmp);
        n_rep = atoi(tmp);
        if (n_rep == 0)
            n_rep = 10;

        *nbp = n_rep * strlen(str);
        for (i = 0; i < n_rep; i++)
            strcat(bpseq[1], str);
    }

    comp_base(bpseq[2], bpseq[1], *nbp);
}

void nrptx(char *o, char x, long n)
{
    long i;

    for (i = 0; i < n; i++)
        o[i] = x;
    o[n] = '\0';
}

void prt_sep(FILE * f1, char x, long n)
{
    char str[BUF512];

    nrptx(str, x, n);
    fprintf(f1, "%s\n", str);
}

void rdalc(long *num, long *nbond, char **asym, long **lkg, double **xyz, char *filnam)
{
    FILE *f1;
    long i;
    char str[BUF512];

    f1 = open_file(filnam, "r");

    fgets(str, sizeof(str), f1);
    sscanf(str, "%ld %*s %ld", num, nbond);

    for (i = 1; i <= *num; i++) {
        fgets(str, sizeof(str), f1);
        sscanf(str, "%*d %s %lf %lf %lf", asym[i], &xyz[i][1], &xyz[i][2], &xyz[i][3]);
    }

    for (i = 1; i <= *nbond; i++) {
        fgets(str, sizeof(str), f1);
        sscanf(str, "%*d %ld %ld", &lkg[i][1], &lkg[i][2]);
    }

    fclose(f1);
}

void wrtalc(long num, long nbond, char **asym, long **lkg, double **xyz, char *filnam)
{
    FILE *f1;
    long i;

    f1 = open_file(filnam, "w");

    fprintf(f1, "%5ld ATOMS, %5ld BONDS\n", num, nbond);

    for (i = 1; i <= num; i++)
        fprintf(f1, "%5ld %2s   %9.4f%9.4f%9.4f\n", i, asym[i],
                xyz[i][1], xyz[i][2], xyz[i][3]);

    for (i = 1; i <= nbond; i++)
        fprintf(f1, "%5ld%6ld%6ld\n", i, lkg[i][1], lkg[i][2]);

    fclose(f1);
}

void rdpdb(long *num, char **asym, char *btype, char *strand,
           long *bnum, double **xyz, char *filnam)
{
    FILE *f1;
    long n = 0;
    char str[BUF512], tmp[BUF512];

    f1 = open_file(filnam, "r");

    while (fgets(str, sizeof(str), f1) != NULL) {
        upper(str, strlen(str));
        if (strncmp(str, "END", 3) == 0)
            break;
        if (strncmp(str, "ATOM", 4) == 0) {
            n++;
            strncpy(asym[n], str + 13, 3);
            asym[n][3] = '\0';
            btype[n] = str[19];
            strand[n] = str[21];
            tmp[0] = '\0';
            bnum[n] = atoi(strncat(tmp, str + 23, 3));
            tmp[0] = '\0';
            xyz[n][1] = atof(strncat(tmp, str + 30, 8));
            tmp[0] = '\0';
            xyz[n][2] = atof(strncat(tmp, str + 38, 8));
            tmp[0] = '\0';
            xyz[n][3] = atof(strncat(tmp, str + 46, 8));
        }
    }

    fclose(f1);

    *num = n;
}

void wrtpdb(long num, char **asym, char *btype, char *strand,
            long *bnum, double **xyz, char *filnam)
{
    FILE *f1;
    long i;

    f1 = open_file(filnam, "w");

    fprintf(f1,
            "HEADER  Nucleic Acids Structure Generated by Juergen Walther\n");
    fprintf(f1, "AUTHOR  Program written by Juergen Walther\n");

    for (i = 1; i <= num; i++)
        fprintf(f1, "ATOM  %5ld  %3s  D%c %c%4ld    %8.3f%8.3f%8.3f\n", i, asym[i],
                btype[i], strand[i], bnum[i], xyz[i][1], xyz[i][2], xyz[i][3]);

    fprintf(f1, "END\n");
    fclose(f1);
}

void wrt_my_pdb(long num, char **asym, char *btype, char *strand,
            long *bnum, double **xyz, char *filnam, double **xyz_float)
{
    FILE *f1;
    long i,j;

    f1 = open_file(filnam, "w");

    fprintf(f1,
            "HEADER  Nucleic Acids Structure Generated by Juergen Walther\n");
    fprintf(f1, "AUTHOR  Program written by Juergen Walther\n");

    for (i = 1; i <= num; i++)
        fprintf(f1, "ATOM  %5ld  %3s  D%c %c%4ld    %8.3f%8.3f%8.3f\n", i, asym[i],
                btype[i], strand[i], bnum[i], xyz[i][1], xyz[i][2], xyz[i][3]);

	j=0;
   for (i=num+1;i<=num+30;i++)
	{
        fprintf(f1, "ATOM  %5ld  N1   DC A   1    %8.3f%8.3f%8.3f\n", i, xyz_float[j][0], xyz_float[j][1], xyz_float[j][2]);
	j=j+1;
	}

    fprintf(f1, "END\n");
    fclose(f1);
}


void wrt_my_pdb_from_pdbtot(int num, char **pdb_tot, char *filnam, double **xyz_float, int num_fl)
{
    FILE *f1;
    long i,j;

    f1 = open_file(filnam, "w");

    fprintf(f1,
            "HEADER  Nucleic Acids Structure Generated by Juergen Walther\n");
    fprintf(f1, "AUTHOR  Program written by Juergen Walther\n");

    for (i = 1; i <= num; i++)
        fprintf(f1, "%s", pdb_tot[i]);

	j=0;
   for (i=num+1;i<=num+num_fl;i++)
	{
        fprintf(f1, "ATOM  %5ld  N1   DC A   1    %8.3f%8.3f%8.3f\n", i, xyz_float[j][0], xyz_float[j][1], xyz_float[j][2]);
	j=j+1;
	}

    fprintf(f1, "END\n");
    fclose(f1);
}


void wrt_my_pdb_from_pdbtot_phosphate(int num, char **pdb_tot, char *filnam, double **pos_ph, int num_bp, char *baset_ph, char *strand_ph, int *base_num_ph)
{
    FILE *f1;
    long i,j,k;
    char atom[4],base[4];	//*atom, *base;
    char strand;
    int base_num, base_num_buf;
    double x,y,z;

    f1 = open_file(filnam, "w");

    fprintf(f1,
            "HEADER  Nucleic Acids Structure Generated by Juergen Walther\n");
    fprintf(f1, "AUTHOR  Program written by Juergen Walther\n");

	j=1;
    for (i = 1; i <= num; i++){

	base_num_buf = base_num;	// store previous value
	sscanf(pdb_tot[i], "%s %s %c %d %lf %lf %lf\n", atom, base, &strand, &base_num, &x,&y,&z);

	if(base_num == base_num_buf + 1)
	{
	for(k=0;k<2*(num_bp-1);k++){
		if(base_num_ph[k] == base_num)
		{
		fprintf(f1, "ATOM  %5ld  P    D%c %c%4d    %8.3f%8.3f%8.3f\n", j, baset_ph[k], strand_ph[k], base_num_ph[k], pos_ph[k][0], pos_ph[k][1], pos_ph[k][2]);
		j++;
		}
					}												
	}			

		        fprintf(f1, "ATOM  %5ld  %s", j, pdb_tot[i]);
			j++;
				    
						}

    fprintf(f1, "END\n");
    fclose(f1);
}


void str_idx(long *n_match, long *idx, char **matx, char *y, long nr)
{
    long i, num = 0;

    for (i = 1; i <= nr; i++) {
        if (strcmp(matx[i], y) == 0) {
            num++;
            idx[num] = i;
        }
    }

    *n_match = num;
}

void set_bp(long *num1, long *num2, char **asym1, char **asym2,
            double **xyz1, double **xyz2, char *bpname,
            double *param, long xdir, char *dir_name)
{
    const long B_NUM = 100;
    long i, j, num, n_match, *idx, posYnum, negYnum;
    double **orien, *pos, **bpxyz, **xyz, **tmp1;
    double *x_axis, *y_axis, *z_axis, *bporg, *bp_normal;
    double dvalue, *pnts_dist;

    char *btype, *strand, str[BUF512];
    long n1, n2, *bnum;

    btype = cvector(1, B_NUM);
    strand = cvector(1, B_NUM);
    bnum = lvector(1, B_NUM);

    strcpy(str, dir_name);
    i = strlen(str);
    str[i] = bpname[0];
    str[i + 1] = '\0';
    strcat(str, ".pdb");
    rdpdb(&n1, asym1, btype, strand, bnum, xyz1, str);
    str[i] = bpname[1];
    rdpdb(&n2, asym2, btype, strand, bnum, xyz2, str);

    *num1 = n1;
    *num2 = n2;
    num = n1 + n2;

    if (xdir == 1) {
        for (i = 1; i <= n1; i++) {
            xyz1[i][1] = -xyz1[i][1];
            xyz1[i][3] = -xyz1[i][3];
        }
        for (i = 1; i <= n2; i++) {
            xyz2[i][1] = -xyz2[i][1];
            xyz2[i][2] = -xyz2[i][2];
        }
    } else {
        for (i = 1; i <= n2; i++) {
            xyz2[i][2] = -xyz2[i][2];
            xyz2[i][3] = -xyz2[i][3];
        }
    }

    orien = dmatrix(1, 3, 1, 3);
    pos = dvector(1, 3);
    tmp1 = dmatrix(1, 3, 1, 3);
    bpxyz = dmatrix(1, num, 1, 3);

    bpfunc(orien, pos, param);
    tr_dmatrix(tmp1, orien, 3, 3);
    mul_dmatrix(bpxyz, xyz1, tmp1, n1, 3, 3, 3);

    for (i = 1; i <= n1; i++) {
        for (j = 1; j <= 3; j++)
            bpxyz[i][j] += pos[j];
    }
    for (i = 1; i <= n2; i++) {
        for (j = 1; j <= 3; j++)
            bpxyz[i + n1][j] = xyz2[i][j];
    }

    idx = lvector(1, num);

    str_idx(&n_match, idx, asym1, "N9 ", n1);
    if (n_match == 1)
        str_idx(&n_match, idx, asym1, "C8 ", n1);
    else
        str_idx(&n_match, idx, asym1, "C6 ", n1);
    posYnum = idx[1];

    str_idx(&n_match, idx, asym2, "N9 ", n2);
    if (n_match == 1)
        str_idx(&n_match, idx, asym2, "C8 ", n2);
    else
        str_idx(&n_match, idx, asym2, "C6 ", n2);
    negYnum = n1 + idx[1];

    y_axis = dvector(1, 3);
    bporg = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        y_axis[i] = bpxyz[posYnum][i] - bpxyz[negYnum][i];
        bporg[i] = 0.5 * (bpxyz[posYnum][i] + bpxyz[negYnum][i]);
    }
    norm_dvector(y_axis, y_axis, 3);

    n_match = 0;
    for (i = 1; i <= n1; i++) {
        if ((asym1[i][0] != 'H') && (asym1[i][2] != '\'') &&
            (asym1[i][0] != 'P') && (asym1[i][2] != 'P')) {
            n_match++;
            idx[n_match] = i;
        }
    }

    for (i = 1; i <= n2; i++) {
        if ((asym2[i][0] != 'H') && (asym2[i][2] != '\'') &&
            (asym2[i][0] != 'P') && (asym2[i][2] != 'P')) {
            n_match++;
            idx[n_match] = i + n1;
        }
    }

    xyz = dmatrix(1, num, 1, 3);

    for (i = 1; i <= n_match; i++) {
        for (j = 1; j <= 3; j++)
            xyz[i][j] = bpxyz[idx[i]][j];
    }

    bp_normal = dvector(1, 3);
    z_axis = dvector(1, 3);
    x_axis = dvector(1, 3);
    pnts_dist = dvector(1, n_match);

    lsplane(bp_normal, &dvalue, pnts_dist, xyz, n_match);

    vec_orth(z_axis, bp_normal, y_axis);
    cross_dvector(x_axis, y_axis, z_axis);

    for (i = 1; i <= num; i++) {
        for (j = 1; j <= 3; j++)
            bpxyz[i][j] -= bporg[j];
        dvalue = 0;
        for (j = 1; j <= 3; j++)
            dvalue += bpxyz[i][j] * x_axis[j];
        xyz[i][1] = dvalue;
        dvalue = 0;
        for (j = 1; j <= 3; j++)
            dvalue += bpxyz[i][j] * y_axis[j];
        xyz[i][2] = dvalue;
        dvalue = 0;
        for (j = 1; j <= 3; j++)
            dvalue += bpxyz[i][j] * z_axis[j];
        xyz[i][3] = dvalue;
    }

    for (i = 1; i <= n1; i++) {
        for (j = 1; j <= 3; j++)
            xyz1[i][j] = xyz[i][j];
    }

    for (i = 1; i <= n2; i++) {
        for (j = 1; j <= 3; j++)
            xyz2[i][j] = xyz[i + n1][j];
    }

    free_lvector(idx, 1, num);

    free_dmatrix(orien, 1, 3, 1, 3);
    free_dvector(pos, 1, 3);
    free_dmatrix(bpxyz, 1, num, 1, 3);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_dmatrix(tmp1, 1, 3, 1, 3);

    free_dvector(x_axis, 1, 3);
    free_dvector(y_axis, 1, 3);
    free_dvector(z_axis, 1, 3);
    free_dvector(bporg, 1, 3);
    free_dvector(bp_normal, 1, 3);

    free_dvector(pnts_dist, 1, n_match);

    free_cvector(btype, 1, B_NUM);
    free_cvector(strand, 1, B_NUM);
    free_lvector(bnum, 1, B_NUM);
}
