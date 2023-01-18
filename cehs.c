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

int main(long argc, char *argv[])
{
    FILE *fo;
    char inpfil[BUF512], filstr[BUF512], outfil[BUF512], *ptr;
    long ipos;

    /* Atom characterization */
    long *RY_idx, *base_num, num_atom, num_bp, n_base;
    char **atom_id, **bp_seq, *base_type;
    double **xyz;

    /* Get atomic listing */
    long *C68, *base;

    /* Reference triads */
    double **bx, **by, **bz, **borg;

    if (argc == 1) {
        printf("Input file name: ");
        read_stdin_str(inpfil);
    } else
        strcpy(inpfil, argv[1]);

    ptr = strchr(inpfil, '.');
    if (ptr != NULL) {
        ipos = ptr - inpfil + 1;
        strncpy(filstr, inpfil, (size_t) ipos);
        filstr[ipos] = '\0';
    } else {
        strcpy(filstr, inpfil);
        strcat(filstr, ".\0");
    }

    strcpy(outfil, filstr);
    fo = open_file(strcat(outfil, "out"), "w");

    RY_idx = lvector(1, MAX_ATOM);
    atom_id = cmatrix(1, MAX_ATOM, 1, 4);
    bp_seq = cmatrix(1, MAX_ATOM, 1, 3);
    xyz = dmatrix(1, MAX_ATOM, 1, 3);
    base_type = cvector(1, MAX_ATOM);
    base_num = lvector(1, MAX_ATOM);

    /* Read in the data file and get the information about each atom */
    get_data_lite(RY_idx, atom_id, bp_seq, xyz, base_type, base_num,
                  &num_atom, &num_bp, inpfil, fo);

    /* Get atomic list */
    C68 = lvector(1, 2 * num_bp);
    base = lvector(1, num_atom);  /* Nothing special */

    /* Get the sorted atom listing */
    get_list_lite(C68, base, &n_base, RY_idx, atom_id, num_bp, num_atom);

    /* Get reference triad of each base and base-pair */
    bx = dmatrix(1, num_bp, 1, 9);
    by = dmatrix(1, num_bp, 1, 9);
    bz = dmatrix(1, num_bp, 1, 9);
    borg = dmatrix(1, num_bp, 1, 9);
    get_frame_lite(bx, by, bz, borg, RY_idx, C68, num_bp, base, n_base,
                   base_num, atom_id, xyz, bp_seq, fo);

    /* Get local parameters based on CEHS */
    CEHS_par_lite(bx, by, bz, borg, num_bp, bp_seq, filstr, fo);

    fclose(fo);

    free_lvector(RY_idx, 1, MAX_ATOM);
    free_cvector(base_type, 1, MAX_ATOM);
    free_lvector(base_num, 1, MAX_ATOM);
    free_lvector(C68, 1, 2 * num_bp);
    free_lvector(base, 1, num_atom);
    free_cmatrix(atom_id, 1, MAX_ATOM, 1, 4);
    free_cmatrix(bp_seq, 1, MAX_ATOM, 1, 3);
    free_dmatrix(xyz, 1, MAX_ATOM, 1, 3);
    free_dmatrix(bx, 1, num_bp, 1, 9);
    free_dmatrix(by, 1, num_bp, 1, 9);
    free_dmatrix(bz, 1, num_bp, 1, 9);
    free_dmatrix(borg, 1, num_bp, 1, 9);

    return 0;
}

void get_data_lite(long *RY_idx, char **atom_id, char **bp_seq, double **xyz,
                   char *base_type, long *base_num, long *num_atom, long *num_bp,
                   char *inpfil, FILE * fo)
{
    FILE *fi;
    char str[BUF512], asym[4], btype;
    long bnum, n = 0, ns = 0, ne = 0, N_purine, i, j, k;
    double x, y, z;

    long *bidx, *bsidx, *beidx, *N9_idx, *tmp1, *tmp2;
    long ncircle, inum, icol, ist, max_base = 21, tnum = 100;
    char *tmp_str;

    /* Print out the date & time of running this program */
    time_t run_time;

    fi = open_file(inpfil, "r");

    while (fgets(str, sizeof(str), fi) != NULL) {
        upper(str, strlen(str));
        sscanf(str, "%s %c %ld %lf %lf %lf", asym, &btype, &bnum, &x, &y, &z);

        if (asym[0] != 'H') {
            n++;
            i = strlen(asym);
            if (i < 1 || i > 3)
                nrerror("Wrong atomic idenfication!");
            if (i == 1)
                strcat(asym, "  \0");
            if (i == 2)
                strcat(asym, " \0");

            if (asym[2] == '*')
                asym[2] = '\'';
            if (strcmp(asym, "O1'") == 0)
                strcpy(asym, "O4'");
            if (strcmp(asym, "OL ") == 0)
                strcpy(asym, "O1P");
            if (strcmp(asym, "OR ") == 0)
                strcpy(asym, "O2P");

            strcpy(atom_id[n], asym);
            base_type[n] = btype;
            base_num[n] = bnum;
            xyz[n][1] = x;
            xyz[n][2] = y;
            xyz[n][3] = z;
        }
    }
    fclose(fi);

    if (!n)
        nrerror("Empty input file!");

    bidx = lvector(1, n);
    bsidx = lvector(1, n);
    beidx = lvector(1, n);
    N9_idx = lvector(1, n);
    tmp1 = lvector(1, n - 1);
    tmp2 = lvector(1, n + 1);

    for (i = 1; i <= n; i++)
        bidx[i] = base_num[i] + 100 * base_type[i];  /* char as integer & use bidx temp */

    i_diff(tmp1, bidx, n);
    logical_not(tmp1, tmp1, n - 1);

    tmp2[1] = 0;
    tmp2[n + 1] = 0;
    for (i = 1; i <= n - 1; i++)
        tmp2[i + 1] = tmp1[i];

    i_diff(bidx, tmp2, n + 1);

    for (i = 1; i <= n; i++) {
        if (bidx[i] == 1) {
            ns++;
            bsidx[ns] = i;
        }
        if (bidx[i] == -1) {
            ne++;
            beidx[ne] = i;
        }
    }

    if (ns % 2)
        nrerror("The two strands are not of equal length!");

    str_idx(&N_purine, N9_idx, atom_id, "N9 ", n);

    zero_lvector(RY_idx, n);
    for (i = 1; i <= N_purine; i++) {
        for (j = 1; j <= ns; j++) {
            if ((N9_idx[i] >= bsidx[j]) && (N9_idx[i] <= beidx[j]))
                break;
        }
        for (k = bsidx[j]; k <= beidx[j]; k++)
            RY_idx[k] = 1;
    }

    *num_bp = ns / 2;
    *num_atom = n;
    k = *num_bp;  /* K used as num_bp for short */

    if (n > tnum)
        tnum = n;
    tmp_str = cvector(0, tnum);

    for (i = 1; i <= k; i++) {
        bp_seq[i][1] = base_type[bsidx[i]];
        bp_seq[i][2] = base_type[bsidx[2 * k - i + 1]];
    }

    /* Write a short note about the program and the structure */
    prt_sep(fo, '*', 76);
    nrptx(tmp_str, ' ', 14);
    fprintf(fo, "%sMain output file from SCHNAaP (Version 1.1, Aug. 1998)\n\n", tmp_str);
    fprintf(fo, "          An analysis program for double-helical");
    fprintf(fo, " nucleic acid structures\n");
    nrptx(tmp_str, ' ', 37);
    fprintf(fo, "\n%sby\n\n", tmp_str);
    nrptx(tmp_str, ' ', 14);
    fprintf(fo, "%sXiang-Jun Lu & C. A. Hunter", tmp_str);
    fprintf(fo, " (University of Sheffield)\n\n");
    nrptx(tmp_str, ' ', 20);
    fprintf(fo, "%sM. A. El Hassan (University of Cambridge)\n\n", tmp_str);
    prt_sep(fo, '*', 76);

    fprintf(fo, "1. The list of the parameters given below are along the 5' to 3' \n");
    fprintf(fo, "   direction of strand I and 3' to 5' direction of strand II\n\n");
    fprintf(fo, "2. All torsion angles are given in the range of [-180, +180]\n\n");

    prt_sep(fo, '*', 76);

    run_time = time(NULL);

    fprintf(fo, "Name of input data file: %s\n", inpfil);
    fprintf(fo, "Date and time of this analysis: ");
    fprintf(fo, "%s", ctime(&run_time));

    fprintf(fo, "\nNo. of heavy atoms: %4ld\n", n);
    fprintf(fo, "No. of base-pairs:  %4ld\n", k);
    fprintf(fo, "\nBase pair sequence: \n");

    ncircle = ceil(k / (double) max_base);  /* max_base (21) bases per line maximum */

    for (inum = 1; inum <= ncircle; inum++) {
        if (inum == ncircle) {
            icol = k % max_base;  /* Remainder */
            if (!icol)
                icol = max_base;
        } else
            icol = max_base;

        ist = (inum - 1) * max_base;  /* Starting column */

        fprintf(fo, "           ");
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%3ld", i);
        fprintf(fo, "\n");

        fprintf(fo, "Strand  I: ");
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%3c", bp_seq[i][1]);
        fprintf(fo, "\n");

        fprintf(fo, "           ");
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%3c", '|');
        fprintf(fo, "\n");

        fprintf(fo, "Strand II: ");
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%3c", bp_seq[i][2]);
        fprintf(fo, "\n\n");
    }

    free_lvector(bidx, 1, n);
    free_lvector(bsidx, 1, n);
    free_lvector(beidx, 1, n);
    free_lvector(N9_idx, 1, n);
    free_lvector(tmp1, 1, n - 1);
    free_lvector(tmp2, 1, n + 1);
    free_cvector(tmp_str, 0, tnum);
}

void get_list_lite(long *C68, long *base, long *n_base, long *RY_idx, char **atom_id,
                   long num_bp, long num_atom)
{
    long i, j, lbn, nR = 0, nY = 0;
    long *Ridx, *Yidx;
    char **Ratom, **Yatom;

    long nRC8, *RC8_idx, *RC8;  /* For RC8 atoms numbering */
    long nYC6, *YC6_idx, *YC6;  /* For YC6 atoms numbering */
    long nC68;

    lbn = 2 * num_bp;  /* The sequential no# of last base */

    *n_base = 0;

    /* Get R and Y atom symbols */
    Ridx = lvector(1, num_atom);
    Yidx = lvector(1, num_atom);
    Ratom = cmatrix(1, num_atom, 1, 4);
    Yatom = cmatrix(1, num_atom, 1, 4);

    for (i = 1; i <= num_atom; i++) {
        if (RY_idx[i] == 1) {
            nR++;
            Ridx[nR] = i;
            strcpy(Ratom[nR], atom_id[i]);
        } else {
            nY++;
            Yidx[nY] = i;
            strcpy(Yatom[nY], atom_id[i]);
        }
    }

    /* Get the atom seq# list for and RC8/YC6 */
    RC8_idx = lvector(1, nR);
    str_idx(&nRC8, RC8_idx, Ratom, "C8 ", nR);
    RC8 = lvector(1, nRC8);
    for (i = 1; i <= nRC8; i++)
        RC8[i] = Ridx[RC8_idx[i]];

    YC6_idx = lvector(1, nY);
    str_idx(&nYC6, YC6_idx, Yatom, "C6 ", nY);
    YC6 = lvector(1, nYC6);
    for (i = 1; i <= nYC6; i++)
        YC6[i] = Yidx[YC6_idx[i]];

    nC68 = nRC8 + nYC6;
    if (nC68 != lbn)
        nrerror("RC8/YC6 atom #s does not match #s of base-pairs!");
    for (i = 1; i <= nRC8; i++)
        C68[i] = RC8[i];
    j = nRC8;
    for (i = 1; i <= nYC6; i++)
        C68[j + i] = YC6[i];
    isort(nC68, C68);

    /* =============== base atoms =============== */
    for (i = 1; i <= num_atom; i++) {
        if ((atom_id[i][0] != 'P') &&
            (atom_id[i][2] != 'P') && (atom_id[i][2] != '\'') && (atom_id[i][0] != 'H')) {
            (*n_base)++;
            base[*n_base] = i;
        }
    }

    /* Free all the temporary integer vectors */
    free_lvector(Ridx, 1, num_atom);
    free_lvector(Yidx, 1, num_atom);
    free_lvector(RC8_idx, 1, nR);
    free_lvector(RC8, 1, nRC8);
    free_lvector(YC6_idx, 1, nY);
    free_lvector(YC6, 1, nYC6);
    free_cmatrix(Ratom, 1, num_atom, 1, 4);
    free_cmatrix(Yatom, 1, num_atom, 1, 4);
}

void get_frame_lite(double **bx, double **by, double **bz, double **borg,
                    long *RY_idx, long *C68, long num_bp, long *base, long n_base,
                    long *base_num, char **atom_id, double **xyz, char **bp_seq,
                    FILE * fo)
/* GET_FRAME_LITE get the reference frame for each base & base-pair
        bx,by,bz--(x,y,z)-axes of the reference triads
        borg--xyz coordinates of the reference triads */
{
    long i, j, k, lbn;
    long *tmp1, *tmp2, *tmp3, *idx;
    long *base_bn, *base_en, ns = 0, ne = 0;
    long *base_bn_1, *base_en_1, *base_bn_2, *base_en_2;

    double **bnorm_1, **bnorm_2, **bpnorm, *tmp_vec;  /* normal vector */
    double **org_1, **org_2, **bporg;  /* orgin */
    double **bx_1, **bx_2, **bpx;  /* x-axis */
    double **by_1, **by_2, **bpy;  /* y-axis */
    double **bz_1, **bz_2, **bpz;  /* z-axis */

    long *bnum_1, *bnum_2, n1, n2, n;
    double **b1xyz, **b2xyz, **bpxyz, *tmp_dist, mrise;

    long *hexag_1, *hexag_2;

    char *non_WC;
    double *bp1_dir, *bp2_dir;
    long sum_sign = 0;

    /* --------------------------------------------------------------------
       This initialization methods works for ANSI-C. But for unknown
       reasons, it does NOT work on HP computers ......
       char *hexag_atom[7]={"XXX", "N1 ", "C2 ", "N3 ", "C4 ", "C5 ", "C6 "};
       -------------------------------------------------------------------- */

    /* The following stupid method should always work */
    char *hexag_atom[7];
    hexag_atom[0] = "XXX";
    hexag_atom[1] = "N1 ";
    hexag_atom[2] = "C2 ";
    hexag_atom[3] = "N3 ";
    hexag_atom[4] = "C4 ";
    hexag_atom[5] = "C5 ";
    hexag_atom[6] = "C6 ";
    /* Stop here ====================================== */

    lbn = 2 * num_bp;  /* Last base number */

    /* Get the indexing of the first base-pair & its normal vector */
    tmp1 = lvector(1, n_base - 1);
    tmp2 = lvector(1, n_base + 1);
    tmp3 = lvector(1, n_base);
    idx = lvector(1, n_base);

    /* Good for handling bases with continuous sequential # */
    for (i = 1; i <= n_base; i++)
        tmp3[i] = base[i] + 100 * base_num[base[i]];

    i_diff(tmp1, tmp3, n_base);
    for (i = 1; i <= n_base - 1; i++)
        tmp1[i] = tmp1[i] - 1;
    logical_not(tmp1, tmp1, n_base - 1);
    tmp2[1] = 0;
    tmp2[n_base + 1] = 0;
    for (i = 1; i <= n_base - 1; i++)
        tmp2[i + 1] = tmp1[i];
    i_diff(idx, tmp2, n_base + 1);

    base_bn = lvector(1, lbn);
    base_en = lvector(1, lbn);

    for (i = 1; i <= n_base; i++) {
        if (idx[i] == 1) {
            ns++;
            if (ns > lbn)
                nrerror("Something wrong with naming in you PDB file!");
            base_bn[ns] = base[i];
        }
        if (idx[i] == -1) {
            ne++;
            if (ne > lbn)
                nrerror("Something wrong with naming in you PDB file!");
            base_en[ne] = base[i];
        }
    }

    /* Note: ns should always equal to ne */

    if (ns != lbn)
        nrerror("Something wrong with naming in you PDB file!");

    /* Find the atomic indexing for each base on strands I & II */
    base_bn_1 = lvector(1, num_bp);
    base_en_1 = lvector(1, num_bp);
    base_bn_2 = lvector(1, num_bp);
    base_en_2 = lvector(1, num_bp);

    for (i = 1; i <= num_bp; i++) {
        k = lbn - i + 1;  /* Matching pair: i <===> 2*num_bp-i+1 */
        base_bn_1[i] = base_bn[i];
        base_en_1[i] = base_en[i];
        base_bn_2[i] = base_bn[k];
        base_en_2[i] = base_en[k];
    }

    /* Initialize the arrays for base/base-pair normals, y-axes and origins */
    bnorm_1 = dmatrix(1, num_bp, 1, 3);
    bnorm_2 = dmatrix(1, num_bp, 1, 3);
    bpnorm = dmatrix(1, num_bp, 1, 3);

    org_1 = dmatrix(1, num_bp, 1, 3);
    org_2 = dmatrix(1, num_bp, 1, 3);
    bporg = dmatrix(1, num_bp, 1, 3);

    bx_1 = dmatrix(1, num_bp, 1, 3);
    bx_2 = dmatrix(1, num_bp, 1, 3);
    bpx = dmatrix(1, num_bp, 1, 3);

    by_1 = dmatrix(1, num_bp, 1, 3);
    by_2 = dmatrix(1, num_bp, 1, 3);
    bpy = dmatrix(1, num_bp, 1, 3);

    bz_1 = dmatrix(1, num_bp, 1, 3);
    bz_2 = dmatrix(1, num_bp, 1, 3);
    bpz = dmatrix(1, num_bp, 1, 3);

    hexag_1 = lvector(1, 6);
    hexag_2 = lvector(1, 6);
    tmp_vec = dvector(1, 3);

    non_WC = cvector(1, num_bp);
    bp1_dir = dvector(1, 3);
    bp2_dir = dvector(1, 3);

    /* Note: Set the origin of base-triad to be the mid-point
       of N1--C4 of R and N3--C6 of Y. */
    for (i = 1; i <= num_bp; i++) {
        /* Strand I base */
        n1 = base_en_1[i] - base_bn_1[i] + 1;  /* Number of atoms in strand I base */
        bnum_1 = lvector(1, n1);
        for (j = base_bn_1[i]; j <= base_en_1[i]; j++) {
            k = j - base_bn_1[i] + 1;  /* Indexing from 1...n1 */
            bnum_1[k] = j;
        }
        b1xyz = dmatrix(1, n1, 1, 3);  /* Get its xyz coordinates */
        for (j = 1; j <= n1; j++)
            for (k = 1; k <= 3; k++)
                b1xyz[j][k] = xyz[bnum_1[j]][k];

        /* Strand II base */
        n2 = base_en_2[i] - base_bn_2[i] + 1;  /* Number of atoms in strand II base */
        bnum_2 = lvector(1, n2);
        for (j = base_bn_2[i]; j <= base_en_2[i]; j++) {
            k = j - base_bn_2[i] + 1;  /* Indexing from 1...n2 */
            bnum_2[k] = j;
        }
        b2xyz = dmatrix(1, n2, 1, 3);  /* Get its xyz coordinates */
        for (j = 1; j <= n2; j++)
            for (k = 1; k <= 3; k++)
                b2xyz[j][k] = xyz[bnum_2[j]][k];

        /* Get xyz coordinates for the base-pair */
        n = n1 + n2;
        bpxyz = dmatrix(1, n, 1, 3);
        for (j = 1; j <= n1; j++)  /* base 1 */
            for (k = 1; k <= 3; k++)
                bpxyz[j][k] = b1xyz[j][k];
        for (j = 1; j <= n2; j++)  /* base 2 */
            for (k = 1; k <= 3; k++)
                bpxyz[n1 + j][k] = b2xyz[j][k];

        /* Get base/base-pair normals */
        tmp_dist = dvector(1, n);
        lsplane(bnorm_1[i], &mrise, tmp_dist, b1xyz, n1);
        lsplane(bnorm_2[i], &mrise, tmp_dist, b2xyz, n2);
        lsplane(bpnorm[i], &mrise, tmp_dist, bpxyz, n);

        /* Get the listing of hexagonal atoms in bases */
        for (j = 1; j <= 6; j++) {
            for (k = 1; k <= n1; k++)
                if (!strcmp(atom_id[bnum_1[k]], hexag_atom[j]))
                    break;
            if (k > n1)
                nrerror("There are missing hexagonal atoms with strand I base.");
            hexag_1[j] = bnum_1[k];

            for (k = 1; k <= n2; k++)
                if (!strcmp(atom_id[bnum_2[k]], hexag_atom[j]))
                    break;
            if (k > n2)
                nrerror("There are missing hexagonal atoms with strand II base.");
            hexag_2[j] = bnum_2[k];
        }

        /* Strand I base */
        if (RY_idx[base_bn_1[i]] == 1) {  /* purine base: N1---> C4 */
            for (j = 1; j <= 3; j++)  /* y-axis */
                tmp_vec[j] = xyz[hexag_1[4]][j] - xyz[hexag_1[1]][j];
            norm_dvector(by_1[i], tmp_vec, 3);
            for (j = 1; j <= 3; j++)  /* origin */
                org_1[i][j] = 0.5 * (xyz[hexag_1[4]][j] + xyz[hexag_1[1]][j]);
            for (j = 1; j <= 3; j++)  /* C2--->N1 vector */
                tmp_vec[j] = xyz[hexag_1[1]][j] - xyz[hexag_1[2]][j];
        } else {  /* pyrimidine base: N3--->C6 */
            for (j = 1; j <= 3; j++)  /* y-axis */
                tmp_vec[j] = xyz[hexag_1[6]][j] - xyz[hexag_1[3]][j];
            norm_dvector(by_1[i], tmp_vec, 3);
            for (j = 1; j <= 3; j++)  /* origin */
                org_1[i][j] = 0.5 * (xyz[hexag_1[6]][j] + xyz[hexag_1[3]][j]);
            for (j = 1; j <= 3; j++)  /* C2--->N3 vector */
                tmp_vec[j] = xyz[hexag_1[3]][j] - xyz[hexag_1[2]][j];
        }
        /* Decide the direction of base I */
        cross_dvector(bp1_dir, tmp_vec, by_1[i]);
        if (dot_dvector(bp1_dir, bnorm_1[i], 3) < 0)
            for (j = 1; j <= 3; j++)
                bnorm_1[i][j] = -bnorm_1[i][j];
        if (dot_dvector(bp1_dir, bpnorm[i], 3) < 0)  /* BP follows base I */
            for (j = 1; j <= 3; j++)
                bpnorm[i][j] = -bpnorm[i][j];

        /* Strand II base */
        if (RY_idx[base_bn_2[i]] == 1) {  /* purine base: C4--->N1 */
            for (j = 1; j <= 3; j++)  /* y-axis */
                tmp_vec[j] = xyz[hexag_2[1]][j] - xyz[hexag_2[4]][j];
            norm_dvector(by_2[i], tmp_vec, 3);
            for (j = 1; j <= 3; j++)  /* origin */
                org_2[i][j] = 0.5 * (xyz[hexag_2[1]][j] + xyz[hexag_2[4]][j]);
            for (j = 1; j <= 3; j++)  /* N1--->C2 vector */
                tmp_vec[j] = xyz[hexag_2[2]][j] - xyz[hexag_2[1]][j];
        } else {  /* pyrimidine base: C6--->N3 */
            for (j = 1; j <= 3; j++)  /* y-axis */
                tmp_vec[j] = xyz[hexag_2[3]][j] - xyz[hexag_2[6]][j];
            norm_dvector(by_2[i], tmp_vec, 3);
            for (j = 1; j <= 3; j++)  /* origin */
                org_2[i][j] = 0.5 * (xyz[hexag_2[3]][j] + xyz[hexag_2[6]][j]);
            for (j = 1; j <= 3; j++)  /* N3--->C2 vector */
                tmp_vec[j] = xyz[hexag_2[2]][j] - xyz[hexag_2[3]][j];
        }
        /* Decide the direction of base II */
        cross_dvector(bp2_dir, by_2[i], tmp_vec);
        if (dot_dvector(bp2_dir, bnorm_2[i], 3) < 0)
            for (j = 1; j <= 3; j++)
                bnorm_2[i][j] = -bnorm_2[i][j];

        /* Check the faces of the two bases */
        if (dot_dvector(bp1_dir, bp2_dir, 3) < 0)
            non_WC[i] = '*';
        else
            non_WC[i] = ' ';

        /* Free temporary array & matrix */
        free_lvector(bnum_1, 1, n1);
        free_lvector(bnum_2, 1, n2);
        free_dmatrix(b1xyz, 1, n1, 1, 3);
        free_dmatrix(b2xyz, 1, n2, 1, 3);
        free_dmatrix(bpxyz, 1, n, 1, 3);
        free_dvector(tmp_dist, 1, n);
    }

    /* Get the base-pair's y-axis and origin */
    for (i = 1; i <= num_bp; i++) {
        k = lbn - i + 1;
        for (j = 1; j <= 3; j++) {
            tmp_vec[j] = xyz[C68[i]][j] - xyz[C68[k]][j];
            bporg[i][j] = 0.5 * (xyz[C68[i]][j] + xyz[C68[k]][j]);
        }
        norm_dvector(bpy[i], tmp_vec, 3);
    }

    /* Check to see if BASE-PAIR NORMAL is along 5'--->3' of strand I
       Make sure that the program works for Z-DNA as before */
    for (i = 1; i < num_bp; i++) {
        for (j = 1; j <= 3; j++)
            tmp_vec[j] = bporg[i + 1][j] - bporg[i][j];
        if (dot_dvector(tmp_vec, bpnorm[i], 3) < 0)
            sum_sign += -1;
        else
            sum_sign += 1;
    }

    if (sum_sign < 0) {  /* No clear idea yet: Change ALL signs (Z-DNA) */
        for (i = 1; i <= num_bp; i++)
            for (j = 1; j <= 3; j++) {
                bnorm_1[i][j] = -bnorm_1[i][j];
                bnorm_2[i][j] = -bnorm_2[i][j];
                bpnorm[i][j] = -bpnorm[i][j];
            }
    }

    /* Print out the base-pair origin and normal vector */
    prt_sep(fo, '*', 76);
    fprintf(fo,
            "The origin (Ox, Oy, Oz) and normal vector (Nx, Ny, Nz) of each base-pair\n");
    fprintf(fo, "       Note: A * denotes a non-WC base-pair (%5ld)\n\n", sum_sign);
    fprintf(fo, "                Ox        Oy        Oz        Nx        Ny        Nz\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fo, " %2ld %c-%c%c ", i, bp_seq[i][1], bp_seq[i][2], non_WC[i]);
        for (j = 1; j <= 3; j++)
            fprintf(fo, "%10.2f", bporg[i][j]);
        for (j = 1; j <= 3; j++)
            fprintf(fo, "%10.2f", bpnorm[i][j]);
        fprintf(fo, "\n");
    }
    fprintf(fo, "\n");

    /* Make sure y-axis and the normal are orthogonal */
    for (i = 1; i <= num_bp; i++) {
        vec_orth(bz_1[i], bnorm_1[i], by_1[i]);
        cross_dvector(bx_1[i], by_1[i], bz_1[i]);
        vec_orth(bz_2[i], bnorm_2[i], by_2[i]);
        cross_dvector(bx_2[i], by_2[i], bz_2[i]);
        vec_orth(bpz[i], bpnorm[i], bpy[i]);
        cross_dvector(bpx[i], bpy[i], bpz[i]);
    }

    /* Combine these vectors together:
       Bases on strand I at columns 1..3
       Bases on strand II at columns 4..6
       Base-pairs at columns: 7..9 */
    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++)
            bx[i][j] = bx_1[i][j];
        for (j = 1; j <= 3; j++)
            bx[i][3 + j] = bx_2[i][j];
        for (j = 1; j <= 3; j++)
            bx[i][6 + j] = bpx[i][j];

        for (j = 1; j <= 3; j++)
            by[i][j] = by_1[i][j];
        for (j = 1; j <= 3; j++)
            by[i][3 + j] = by_2[i][j];
        for (j = 1; j <= 3; j++)
            by[i][6 + j] = bpy[i][j];

        for (j = 1; j <= 3; j++)
            bz[i][j] = bz_1[i][j];
        for (j = 1; j <= 3; j++)
            bz[i][3 + j] = bz_2[i][j];
        for (j = 1; j <= 3; j++)
            bz[i][6 + j] = bpz[i][j];

        for (j = 1; j <= 3; j++)
            borg[i][j] = org_1[i][j];
        for (j = 1; j <= 3; j++)
            borg[i][3 + j] = org_2[i][j];
        for (j = 1; j <= 3; j++)
            borg[i][6 + j] = bporg[i][j];
    }

    free_lvector(tmp1, 1, n_base - 1);
    free_lvector(tmp2, 1, n_base + 1);
    free_lvector(tmp3, 1, n_base);
    free_lvector(idx, 1, n_base);
    free_lvector(base_bn, 1, lbn);
    free_lvector(base_en, 1, lbn);
    free_lvector(base_bn_1, 1, num_bp);
    free_lvector(base_en_1, 1, num_bp);
    free_lvector(base_bn_2, 1, num_bp);
    free_lvector(base_en_2, 1, num_bp);
    free_lvector(hexag_1, 1, 6);
    free_lvector(hexag_2, 1, 6);
    free_dvector(tmp_vec, 1, 3);
    free_dmatrix(bnorm_1, 1, num_bp, 1, 3);
    free_dmatrix(bnorm_2, 1, num_bp, 1, 3);
    free_dmatrix(bpnorm, 1, num_bp, 1, 3);
    free_dmatrix(org_1, 1, num_bp, 1, 3);
    free_dmatrix(org_2, 1, num_bp, 1, 3);
    free_dmatrix(bporg, 1, num_bp, 1, 3);
    free_dmatrix(bx_1, 1, num_bp, 1, 3);
    free_dmatrix(bx_2, 1, num_bp, 1, 3);
    free_dmatrix(bpx, 1, num_bp, 1, 3);
    free_dmatrix(by_1, 1, num_bp, 1, 3);
    free_dmatrix(by_2, 1, num_bp, 1, 3);
    free_dmatrix(bpy, 1, num_bp, 1, 3);
    free_dmatrix(bz_1, 1, num_bp, 1, 3);
    free_dmatrix(bz_2, 1, num_bp, 1, 3);
    free_dmatrix(bpz, 1, num_bp, 1, 3);
    free_cvector(non_WC, 1, num_bp);
    free_dvector(bp1_dir, 1, 3);
    free_dvector(bp2_dir, 1, 3);
}

void CEHS_par_lite(double **bx, double **by, double **bz, double **borg,
                   long num_bp, char **bp_seq, char *filstr, FILE * fo)
     /* CEHS_PAR get the local CEHS base-pair and step parameters */
{
    long i, j;

    double **xyz_lb, **xyz_rb, *org_lb, *org_rb;
    double **xyz_ubp, **xyz_lbp, *org_ubp, *org_lbp;
    double **bp_par, **step_par, **mstxyz;

    double *ave_par, *std_par;

    char filnam[BUF512];
    FILE *fc;

    /* Get the base-pair parameters */
    bp_par = dmatrix(1, num_bp, 1, 6);
    xyz_rb = dmatrix(1, 3, 1, 3);
    xyz_lb = dmatrix(1, 3, 1, 3);
    org_rb = dvector(1, 3);
    org_lb = dvector(1, 3);

    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++) {
            xyz_rb[j][1] = bx[i][j + 3];  /* x-axis of base II */
            xyz_rb[j][2] = by[i][j + 3];  /* y-axis of base II */
            xyz_rb[j][3] = bz[i][j + 3];  /* z-axis of base II */
            org_rb[j] = borg[i][j + 3];  /* origin of base II */

            xyz_lb[j][1] = bx[i][j];  /* x-axis of base I */
            xyz_lb[j][2] = by[i][j];  /* y-axis of base I */
            xyz_lb[j][3] = bz[i][j];  /* z-axis of base I */
            org_lb[j] = borg[i][j];  /* origin of base I */
        }
        get_bp_par(bp_par[i], xyz_rb, org_rb, xyz_lb, org_lb);
    }

    /* Get the step parameters */
    step_par = dmatrix(1, num_bp - 1, 1, 6);
    xyz_ubp = dmatrix(1, 3, 1, 3);
    xyz_lbp = dmatrix(1, 3, 1, 3);
    org_ubp = dvector(1, 3);
    org_lbp = dvector(1, 3);

    mstxyz = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_bp - 1; i++) {
        for (j = 1; j <= 3; j++) {
            xyz_lbp[j][1] = bx[i][j + 6];  /* x-axis of lower base-pair */
            xyz_lbp[j][2] = by[i][j + 6];  /* y-axis of lower base-pair */
            xyz_lbp[j][3] = bz[i][j + 6];  /* z-axis of lower base-pair */
            org_lbp[j] = borg[i][j + 6];  /* origin of lower base-pair */

            xyz_ubp[j][1] = bx[i + 1][j + 6];  /* x-axis of upper base-pair */
            xyz_ubp[j][2] = by[i + 1][j + 6];  /* y-axis of upper base-pair */
            xyz_ubp[j][3] = bz[i + 1][j + 6];  /* z-axis of upper base-pair */
            org_ubp[j] = borg[i + 1][j + 6];  /* origin of upper base-pair */
        }

        get_step_par(step_par[i], mstxyz, xyz_lbp, org_lbp, xyz_ubp, org_ubp);
    }

    /* Print out the base-pair and step parameters */
    prt_sep(fo, '*', 76);
    fprintf(fo, "Local CEHS base-pair and step parameters by El Hassan & Calladine\n");

    fprintf(fo, "\nBase-pair parameters\n");
    fprintf(fo, "   bp        Shear    Stretch   Stagger");
    fprintf(fo, "    Buckle  Propeller  Opening\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fo, " %2ld %c-%c ", i, bp_seq[i][1], bp_seq[i][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%10.2f", bp_par[i][j]);
        fprintf(fo, "\n");
    }
    fprintf(fo, "         ");
    prt_sep(fo, '~', 59);

    ave_par = dvector(1, 6);
    std_par = dvector(1, 6);

    mean_dmatrix(ave_par, bp_par, num_bp, 6);
    std_dmatrix(std_par, bp_par, num_bp, 6);

    fprintf(fo, "    ave.");
    for (i = 1; i <= 6; i++)
        fprintf(fo, "%10.2f", ave_par[i]);
    fprintf(fo, "\n    s.d.");
    for (i = 1; i <= 6; i++)
        fprintf(fo, "%10.2f", std_par[i]);
    fprintf(fo, "\n");

    fprintf(fo, "\nStep parameters\n");
    fprintf(fo, "  step       Shift     Slide      Rise");
    fprintf(fo, "      Tilt      Roll     Twist\n");
    for (i = 1; i <= num_bp - 1; i++) {
        fprintf(fo, "%2ld %c%c/%c%c", i, bp_seq[i][1], bp_seq[i + 1][1],
                bp_seq[i + 1][2], bp_seq[i][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%10.2f", step_par[i][j]);
        fprintf(fo, "\n");
    }

    if (num_bp > 2) {
        fprintf(fo, "         ");
        prt_sep(fo, '~', 59);
        mean_dmatrix(ave_par, step_par, num_bp - 1, 6);
        std_dmatrix(std_par, step_par, num_bp - 1, 6);
        fprintf(fo, "    ave.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%10.2f", ave_par[i]);
        fprintf(fo, "\n    s.d.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%10.2f", std_par[i]);
        fprintf(fo, "\n\n");
    }

    /* Write bp & step parameter data file for rebuilding! */
    strcpy(filnam, filstr);
    fc = open_file(strcat(filnam, "ceh"), "w");
    fprintf(fc, "%4ld base-pairs\n", num_bp);
    fprintf(fc, "        Shear    Str    Stag    Buck    Prop    Open");
    fprintf(fc, "   Shift   Slide    Rise    Tilt    Roll   Twist\n");

    i = 1;
    fprintf(fc, "%c-%c ", bp_seq[i][1], bp_seq[i][2]);
    for (j = 1; j <= 6; j++)
        fprintf(fc, "%8.2f", bp_par[i][j]);
    for (j = 1; j <= 6; j++)
        fprintf(fc, "%8.2f", 0.0);
    fprintf(fc, "\n");

    for (i = 2; i <= num_bp; i++) {
        fprintf(fc, "%c-%c ", bp_seq[i][1], bp_seq[i][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fc, "%8.2f", bp_par[i][j]);
        for (j = 1; j <= 6; j++)
            fprintf(fc, "%8.2f", step_par[i - 1][j]);
        fprintf(fc, "\n");
    }
    fclose(fc);

    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dvector(org_lb, 1, 3);
    free_dvector(org_rb, 1, 3);
    free_dvector(org_ubp, 1, 3);
    free_dvector(org_lbp, 1, 3);
    free_dvector(ave_par, 1, 6);
    free_dvector(std_par, 1, 6);
    free_dmatrix(xyz_lb, 1, 3, 1, 3);
    free_dmatrix(xyz_rb, 1, 3, 1, 3);
    free_dmatrix(step_par, 1, num_bp - 1, 1, 6);
    free_dmatrix(xyz_ubp, 1, 3, 1, 3);
    free_dmatrix(xyz_lbp, 1, 3, 1, 3);
    free_dmatrix(mstxyz, 1, 3, 1, 3);
}
