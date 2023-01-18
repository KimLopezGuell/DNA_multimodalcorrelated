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

void get_data(long *RY_idx, char **atom_id, char **bp_seq, double **xyz,
              char *base_type, char *strand, long *base_num,
              long *num_atom, long *num_bp, char *inpfil, FILE * fo)
{
    FILE *fi;
    char str[BUF512], asym[4], btype;
    long bnum, n = 0, ns = 0, ne = 0, N_purine, i, j, k;
    double x, y, z;

    long *bidx, *bsidx, *beidx, *N9_idx, *tmp1, *tmp2;
    long numA, numB, ncircle, inum, icol, ist, max_base = 21, tnum = 100;
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

    numA = bsidx[k + 1] - 1;
    numB = n - numA;

    if (n > tnum)
        tnum = n;
    tmp_str = cvector(0, tnum);

    nrptx(tmp_str, 'A', numA + 1);
    strcpy(strand, tmp_str);
    nrptx(tmp_str, 'B', numB + 1);
    strcat(strand, tmp_str);

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
    fprintf(fo, "%s\n", ctime(&run_time));
    fprintf(fo, "\nNo. of heavy atoms: %4ld\n", n);
    fprintf(fo, "No. of base-pairs:  %4ld\n", k);

    fprintf(fo, "\nBase pair sequence: \n");

    /* max_base (21) bases per line maximum */
    ncircle = (long) ceil(k / (double) max_base);

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

void get_list(long *C1, long *N, long *P, long *O4, long *C68, long *chi,
              long *mc, long *nmc, long *sugar, long *base, long *n_base,
              long *RY_idx, char **atom_id, long num_bp, long num_atom, double **xyz)
{
    long i, j, idx, lbn, nR = 0, nY = 0;
    long *Ridx, *Yidx;
    char **Ratom, **Yatom;

    long nRC1, *RC1_idx, *RC1;  /* For RC1' atoms numbering */
    long nYC1, *YC1_idx, *YC1;  /* For YC1' atoms numbering */
    long nRN9, *RN9_idx, *RN9;  /* For RN9 atoms numbering */
    long nYN1, *YN1_idx, *YN1;  /* For YN1 atoms numbering */
    long nRC8, *RC8_idx, *RC8;  /* For RC8 atoms numbering */
    long nYC6, *YC6_idx, *YC6;  /* For YC6 atoms numbering */
    long nC1, nN, nC68;

    /* Main chain atoms */
    long nO5, *O5, nC5, *C5, nC4, *C4, nC3, *C3, nO3, *O3, nP;

    /* chi torsion angle */
    long nRC4, *RC4_idx, *RC4;  /* For RC4 atoms numbering */
    long nYC2, *YC2_idx, *YC2;  /* For YC2 atoms numbering */
    long *C24;  /* For RC4/YC2 atoms numbering */

    /* sugar torsion angle */
    long nO4, nC2, *C2;

    /* To check for special O3' linkage case */
    double *dd, aveC3_O3, *tmp_vec;
    char yn, str[BUF512];
    long iA, iB;

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

    /* Get the atom seq# list for RC1/YC1, RN9/YN1 and RC8/YC6 */

    /* RC1' atoms */
    RC1_idx = lvector(1, nR);
    str_idx(&nRC1, RC1_idx, Ratom, "C1'", nR);
    RC1 = lvector(1, nRC1);
    for (i = 1; i <= nRC1; i++)
        RC1[i] = Ridx[RC1_idx[i]];

    /* YC1' atoms */
    YC1_idx = lvector(1, nY);
    str_idx(&nYC1, YC1_idx, Yatom, "C1'", nY);
    YC1 = lvector(1, nYC1);
    for (i = 1; i <= nYC1; i++)
        YC1[i] = Yidx[YC1_idx[i]];

    /* Sort all C1' atoms together */
    nC1 = nRC1 + nYC1;
    if (nC1 != lbn)
        nrerror("C1' atom #s does not match #s of base-pairs!");
    for (i = 1; i <= nRC1; i++)
        C1[i] = RC1[i];
    j = nRC1;
    for (i = 1; i <= nYC1; i++)
        C1[j + i] = YC1[i];
    isort(nC1, C1);

    /* ---------------------------------------- */
    /* RN9 atoms */
    RN9_idx = lvector(1, nR);
    str_idx(&nRN9, RN9_idx, Ratom, "N9 ", nR);
    RN9 = lvector(1, nRN9);
    for (i = 1; i <= nRN9; i++)
        RN9[i] = Ridx[RN9_idx[i]];

    /* YN1 atoms */
    YN1_idx = lvector(1, nY);
    str_idx(&nYN1, YN1_idx, Yatom, "N1 ", nY);
    YN1 = lvector(1, nYN1);
    for (i = 1; i <= nYN1; i++)
        YN1[i] = Yidx[YN1_idx[i]];

    /* Sort all RN9-YN1 atoms together */
    nN = nRN9 + nYN1;
    if (nN != lbn)
        nrerror("RN9/YN1 atom #s does not match #s of base-pairs!");
    for (i = 1; i <= nRN9; i++)
        N[i] = RN9[i];
    j = nRN9;
    for (i = 1; i <= nYN1; i++)
        N[j + i] = YN1[i];
    isort(nN, N);

    /* ---------------------------------------- */
    /* RC8 atoms */
    RC8_idx = lvector(1, nR);
    str_idx(&nRC8, RC8_idx, Ratom, "C8 ", nR);
    RC8 = lvector(1, nRC8);
    for (i = 1; i <= nRC8; i++)
        RC8[i] = Ridx[RC8_idx[i]];

    /* YC6 atoms */
    YC6_idx = lvector(1, nY);
    str_idx(&nYC6, YC6_idx, Yatom, "C6 ", nY);
    YC6 = lvector(1, nYC6);
    for (i = 1; i <= nYC6; i++)
        YC6[i] = Yidx[YC6_idx[i]];

    /* Sort all RC8-YC1 atoms together */
    nC68 = nRC8 + nYC6;
    if (nC68 != lbn)
        nrerror("RC8/YC6 atom #s does not match #s of base-pairs!");
    for (i = 1; i <= nRC8; i++)
        C68[i] = RC8[i];
    j = nRC8;
    for (i = 1; i <= nYC6; i++)
        C68[j + i] = YC6[i];
    isort(nC68, C68);

    /* =============== Main chain: O5'-C5'-C4'-C3'-O3'-P =============== */

    O5 = lvector(1, lbn);
    C5 = lvector(1, lbn);
    C4 = lvector(1, lbn);
    C3 = lvector(1, lbn);
    O3 = lvector(1, lbn);

    str_idx(&nO5, O5, atom_id, "O5'", num_atom);
    str_idx(&nC5, C5, atom_id, "C5'", num_atom);
    str_idx(&nC4, C4, atom_id, "C4'", num_atom);
    str_idx(&nC3, C3, atom_id, "C3'", num_atom);
    str_idx(&nO3, O3, atom_id, "O3'", num_atom);
    str_idx(&nP, P, atom_id, "P  ", num_atom);

    if (nO5 != lbn)
        nrerror("O5' atom #s does not match #s of base-pairs!");
    if (nC5 != lbn)
        nrerror("C5' atom #s does not match #s of base-pairs!");
    if (nC4 != lbn)
        nrerror("C4' atom #s does not match #s of base-pairs!");
    if (nC3 != lbn)
        nrerror("C3' atom #s does not match #s of base-pairs!");
    if (nO3 != lbn)
        nrerror("O3' atom #s does not match #s of base-pairs!");

    /*% In some cases, O3(i) is linked to P(i) instead of P(i+1) (Convention)
       %   Here I use the distance between O3'(i)--C3'(i) as a measure to
       %   distinguish between the two, and give a warning message. */
    dd = dvector(1, lbn);
    tmp_vec = dvector(1, 3);
    for (i = 1; i <= lbn; i++) {
        for (j = 1; j <= 3; j++)
            tmp_vec[j] = xyz[O3[i]][j] - xyz[C3[i]][j];
        dd[i] = sqrt(dot_dvector(tmp_vec, tmp_vec, 3));
    }
    aveC3_O3 = mean_dvector(dd, lbn);

    if (aveC3_O3 > 2.8) {
        printf("\nThe the mean C3'(i)--O3'(i) length is (%4.1f) > 2.8 A\n", aveC3_O3);
        printf("It seems this structure has the linkage O3(i)--P(i)\n");
        printf("instead of the conventional (IUPAC-IUB) O3(i)--P(i+1)\n\n");
        printf("Is this the case [y(Dft)/n]: ");
        read_stdin_str(str);
        yn = toupper(str[0]);
        if (yn != 'N') {
            iA = O3[1];
            iB = O3[num_bp + 1];  /* Move the 1st to the last */
            for (i = 1; i < num_bp; i++)
                O3[i] = O3[i + 1];
            O3[num_bp] = iA;
            for (i = num_bp + 2; i < lbn; i++)
                O3[i] = O3[i + 1];
            O3[lbn] = iB;
            printf("\nI have re-assigned O3'' atoms according to IUPAC-IUB\n");
            printf("NB: The backbone torsion angles--alpha & delta--which are\n");
            printf("    associated with O3'' atom, for the 3'' terminal residue\n");
            printf("    along each strand are MEANINGLESS in the output file.\n");
        }
    }
    /* End of special handling of O3' issue */

    *nmc = nO5 + nC5 + nC4 + nC3 + nO3 + nP;  /* 5*lbn+nP */

    if (nP == lbn) {  /* With terminal P */
        for (i = 1; i <= lbn; i++) {
            j = (i - 1) * 6;
            mc[j + 1] = P[i];
            mc[j + 2] = O5[i];
            mc[j + 3] = C5[i];
            mc[j + 4] = C4[i];
            mc[j + 5] = C3[i];
            mc[j + 6] = O3[i];
        }
    } else if (nP == lbn - 2) {
        if (P[1] < O5[1])  /* Missing P in the middle */
            nrerror("There are missing P atoms in the middle of the strucutre!");

        /* For strand I */
        for (i = 1; i <= num_bp; i++) {
            j = (i - 1) * 6;
            mc[j + 1] = O5[i];
            mc[j + 2] = C5[i];
            mc[j + 3] = C4[i];
            mc[j + 4] = C3[i];
            mc[j + 5] = O3[i];
            if (i != num_bp)
                mc[j + 6] = P[i];
        }
        idx = j + 5;

        /* For strand II */
        for (i = num_bp + 1; i <= lbn; i++) {
            j = idx + (i - num_bp - 1) * 6;
            mc[j + 1] = O5[i];
            mc[j + 2] = C5[i];
            mc[j + 3] = C4[i];
            mc[j + 4] = C3[i];
            mc[j + 5] = O3[i];
            if (i != lbn)
                mc[j + 6] = P[i - 1];
        }
    } else
        nrerror("P atom #s does not match #s of base-pairs!");

    /* =============== chi angle atoms =============== */

    /* chi of R: O4'-C1'-N9-C4, Y: O4'-C1'-N1-C2 */
    RC4_idx = lvector(1, nR);
    str_idx(&nRC4, RC4_idx, Ratom, "C4 ", nR);
    RC4 = lvector(1, nRC4);
    for (i = 1; i <= nRC4; i++)
        RC4[i] = Ridx[RC4_idx[i]];

    YC2_idx = lvector(1, nY);
    str_idx(&nYC2, YC2_idx, Yatom, "C2 ", nY);
    YC2 = lvector(1, nYC2);
    for (i = 1; i <= nYC2; i++)
        YC2[i] = Yidx[YC2_idx[i]];

    /* Combine and sort RC4 and YC2 */
    if (nRC4 + nYC2 != lbn)
        nrerror("RC4'/YC2' atom #s does not match #s of base-pairs!");
    C24 = lvector(1, lbn);
    for (i = 1; i <= nRC4; i++)
        C24[i] = RC4[i];
    for (i = 1; i <= nYC2; i++)
        C24[nRC4 + i] = YC2[i];
    isort(lbn, C24);

    str_idx(&nO4, O4, atom_id, "O4'", num_atom);

    if (nO4 != lbn)
        nrerror("O4' atom #s does not match #s of base-pairs!");

    for (i = 1; i <= lbn; i++) {
        j = (i - 1) * 4;
        chi[j + 1] = O4[i];
        chi[j + 2] = C1[i];
        chi[j + 3] = N[i];
        chi[j + 4] = C24[i];
    }

    /* =============== sugar ring atoms =============== */
    C2 = lvector(1, lbn);
    str_idx(&nC2, C2, atom_id, "C2'", num_atom);

    if (nC2 != lbn)
        nrerror("C2' atom #s does not match #s of base-pairs!");

    for (i = 1; i <= lbn; i++) {
        j = (i - 1) * 5;
        sugar[j + 1] = C4[i];
        sugar[j + 2] = O4[i];
        sugar[j + 3] = C1[i];
        sugar[j + 4] = C2[i];
        sugar[j + 5] = C3[i];
    }

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
    free_lvector(RC1_idx, 1, nR);
    free_lvector(RC1, 1, nRC1);
    free_lvector(YC1_idx, 1, nY);
    free_lvector(YC1, 1, nYC1);
    free_lvector(RN9_idx, 1, nR);
    free_lvector(RN9, 1, nRN9);
    free_lvector(YN1_idx, 1, nY);
    free_lvector(YN1, 1, nYN1);
    free_lvector(RC8_idx, 1, nR);
    free_lvector(RC8, 1, nRC8);
    free_lvector(YC6_idx, 1, nY);
    free_lvector(YC6, 1, nYC6);
    free_lvector(O5, 1, lbn);
    free_lvector(C5, 1, lbn);
    free_lvector(C4, 1, lbn);
    free_lvector(C3, 1, lbn);
    free_lvector(O3, 1, lbn);
    free_lvector(RC4_idx, 1, nR);
    free_lvector(RC4, 1, nRC4);
    free_lvector(YC2_idx, 1, nY);
    free_lvector(YC2, 1, nYC2);
    free_lvector(C2, 1, lbn);
    free_cmatrix(Ratom, 1, num_atom, 1, 4);
    free_cmatrix(Yatom, 1, num_atom, 1, 4);
    free_lvector(C24, 1, lbn);
    free_dvector(dd, 1, lbn);
    free_dvector(tmp_vec, 1, 3);
}

void set_str(double **xyz, long *C1, long *N, long *C68, long *base,
             long *base_num, long n_base, long num_atom, long num_bp, FILE * fo)
     /* SET_STR set the structure with regard to the global reference frame
        if it is NOT strongly curved (see also: get_frame)

        Use C1' & RN9/YN1 equivalent atom-pair vectors to define the
        "best-fit" helical axis

        Note: xyz--both as input & output */
{
    long i, j, k, n_vec, nbp1, lbn;
    long *C1_I, *C1_II, *N_I, *N_II, *b_vecNum, *e_vecNum;
    long *tmp1, *tmp2, *idx;
    long *base_bn, *base_en, ns = 0, ne = 0;
    long *bp1_idx, num_bp1;

    double **vec_xyz;
    double *z_axis, mrise, *pnts_dist;
    double **rotmat, **rotmat_T, **xyz2;
    double **bp1_xyz, *bp1_normal, *bp1x, *bp1y, *bp1z, **bp1xyz;
    double *hinge, ang_deg, **glbxyz;
    double **b_xy, **e_xy, **dxy, *g, exy2, bxy2;
    double **t2x2, **dxy_T, **inv_t2x2, *org_xyz;

    nbp1 = num_bp - 1;
    lbn = 2 * num_bp;  /* Last base number */
    if (num_bp < 2)
        nrerror("Not enough vectors to define helical axis!");

    C1_I = lvector(1, num_bp);
    C1_II = lvector(1, num_bp);
    N_I = lvector(1, num_bp);
    N_II = lvector(1, num_bp);

    /* Get vector indexing for each strand in 5'-->3' of strand I */
    for (i = 1; i <= num_bp; i++) {
        j = lbn - i + 1;
        C1_I[i] = C1[i];
        C1_II[i] = C1[j];
        N_I[i] = N[i];
        N_II[i] = N[j];
    }

    n_vec = 4 * nbp1;
    b_vecNum = lvector(1, n_vec);
    e_vecNum = lvector(1, n_vec);
    vec_xyz = dmatrix(1, n_vec, 1, 3);

    /* Vector beginning numbering */
    j = 0;
    for (i = 1; i <= nbp1; i++)
        b_vecNum[j + i] = C1_I[i];
    j += nbp1;
    for (i = 1; i <= nbp1; i++)
        b_vecNum[j + i] = C1_II[i];
    j += nbp1;
    for (i = 1; i <= nbp1; i++)
        b_vecNum[j + i] = N_I[i];
    j += nbp1;
    for (i = 1; i <= nbp1; i++)
        b_vecNum[j + i] = N_II[i];

    /* Vector ending numbering */
    j = 0;
    for (i = 2; i <= num_bp; i++)
        e_vecNum[j + i - 1] = C1_I[i];
    j += nbp1;
    for (i = 2; i <= num_bp; i++)
        e_vecNum[j + i - 1] = C1_II[i];
    j += nbp1;
    for (i = 2; i <= num_bp; i++)
        e_vecNum[j + i - 1] = N_I[i];
    j += nbp1;
    for (i = 2; i <= num_bp; i++)
        e_vecNum[j + i - 1] = N_II[i];

    /* Get the actual vectors for defining the helical axis */
    for (i = 1; i <= n_vec; i++) {
        for (j = 1; j <= 3; j++)
            vec_xyz[i][j] = xyz[e_vecNum[i]][j] - xyz[b_vecNum[i]][j];
    }

    z_axis = dvector(1, 3);
    pnts_dist = dvector(1, n_vec);

    /* Get the "best-fit" helical axis, i.e. z_axis */
    lsplane(z_axis, &mrise, pnts_dist, vec_xyz, n_vec);

    /* Align the structure along +z_axis (5'-->3' of strand I) */
    rotmat = dmatrix(1, 3, 1, 3);
    rotmat_T = dmatrix(1, 3, 1, 3);
    xyz2 = dmatrix(1, num_atom, 1, 3);

    if (mrise < 0) {
        for (i = 1; i <= 3; i++)
            z_axis[i] = -z_axis[i];
    }

    fprintf(fo, "Helical axis: %9.3f%9.3f%9.3f\n\n", z_axis[1], z_axis[2], z_axis[3]);

    alignz(rotmat, z_axis);
    tr_dmatrix(rotmat_T, rotmat, 3, 3);
    mul_dmatrix(xyz2, xyz, rotmat_T, num_atom, 3, 3, 3);

    /* z_axis is now [0 0 1] */
    z_axis[1] = 0;
    z_axis[2] = 0;
    z_axis[3] = 1;

    /* Get the indexing of the first base-pair & its normal vector */
    tmp1 = lvector(1, n_base - 1);
    tmp2 = lvector(1, n_base + 1);
    idx = lvector(1, n_base);

    /* Good for handling bases with continuous sequential # */
    for (i = 1; i <= n_base; i++)
        idx[i] = base[i] + 100 * base_num[base[i]];

    i_diff(tmp1, idx, n_base);
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
    /* Note: ns & ne should be 2*num_bp */

    if (ns != lbn)
        nrerror("Something wrong with naming in you PDB file!");

    num_bp1 = (base_en[1] - base_bn[1] + 1) + (base_en[lbn] - base_bn[lbn] + 1);
    bp1_idx = lvector(1, num_bp1);
    for (i = base_bn[1]; i <= base_en[1]; i++) {
        j = i - base_bn[1] + 1;
        bp1_idx[j] = i;
    }
    k = j;
    for (i = base_bn[lbn]; i <= base_en[lbn]; i++) {
        j = i - base_bn[lbn] + 1;
        bp1_idx[k + j] = i;
    }

    bp1_xyz = dmatrix(1, num_bp1, 1, 3);
    bp1_normal = dvector(1, 3);
    bp1x = dvector(1, 3);
    bp1y = dvector(1, 3);
    bp1z = dvector(1, 3);
    bp1xyz = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_bp1; i++) {
        k = bp1_idx[i];
        for (j = 1; j <= 3; j++)
            bp1_xyz[i][j] = xyz2[k][j];
    }

    /* Get the first base-pair's xyz-axes */
    free_dvector(pnts_dist, 1, n_vec);  /* Free original */
    pnts_dist = dvector(1, num_bp1);  /* Create a new one */
    lsplane(bp1_normal, &mrise, pnts_dist, bp1_xyz, num_bp1);

    for (i = 1; i <= 3; i++)
        bp1y[i] = xyz2[C68[1]][i] - xyz2[C68[lbn]][i];
    norm_dvector(bp1y, bp1y, 3);
    vec_orth(bp1z, bp1_normal, bp1y);
    cross_dvector(bp1x, bp1y, bp1z);
    for (i = 1; i <= 3; i++) {
        bp1xyz[i][1] = bp1x[i];
        bp1xyz[i][2] = bp1y[i];
        bp1xyz[i][3] = bp1z[i];
    }

    /* Align the first base-pair's z-axis with global z-axis */
    hinge = dvector(1, 3);
    glbxyz = dmatrix(1, 3, 1, 3);

    cross_dvector(hinge, z_axis, bp1z);
    ang_deg = magang(z_axis, bp1z);
    arbrot(rotmat, hinge, -ang_deg);
    mul_dmatrix(glbxyz, rotmat, bp1xyz, 3, 3, 3, 3);

    /* Use a least-square method to define global origin's xy coordinates */
    b_xy = dmatrix(1, n_vec, 1, 2);
    e_xy = dmatrix(1, n_vec, 1, 2);
    dxy = dmatrix(1, n_vec, 1, 2);
    g = dvector(1, n_vec);

    for (i = 1; i <= n_vec; i++) {
        exy2 = 0;
        bxy2 = 0;
        for (j = 1; j <= 2; j++) {
            b_xy[i][j] = xyz2[b_vecNum[i]][j];
            e_xy[i][j] = xyz2[e_vecNum[i]][j];
            dxy[i][j] = 2 * (e_xy[i][j] - b_xy[i][j]);
            exy2 = exy2 + e_xy[i][j] * e_xy[i][j];
            bxy2 = bxy2 + b_xy[i][j] * b_xy[i][j];
        }
        g[i] = exy2 - bxy2;
    }

    t2x2 = dmatrix(1, 2, 1, 2);
    dxy_T = dmatrix(1, 2, 1, n_vec);
    inv_t2x2 = dmatrix(1, 2, 1, 2);
    org_xyz = dvector(1, 3);

    tr_dmatrix(dxy_T, dxy, n_vec, 2);
    mul_dmatrix(t2x2, dxy_T, dxy, 2, n_vec, n_vec, 2);
    dinverse(t2x2, 2, inv_t2x2);
    /* Using b_xy as a temporary array */
    mul_dmatrix(b_xy, dxy, inv_t2x2, n_vec, 2, 2, 2);

    zero_dvector(org_xyz, 3);

    for (i = 1; i <= 2; i++)
        for (j = 1; j <= n_vec; j++)
            org_xyz[i] += g[j] * b_xy[j][i];

    /* Set global origin's z coordinate to be that of the first base-pair */
    org_xyz[3] = 0.5 * (xyz2[C68[1]][3] + xyz2[C68[lbn]][3]);

    /* Reset the whole structure in the global reference frame */
    for (i = 1; i <= num_atom; i++)
        for (j = 1; j <= 3; j++)
            xyz2[i][j] = xyz2[i][j] - org_xyz[j];

    mul_dmatrix(xyz, xyz2, glbxyz, num_atom, 3, 3, 3);

    free_lvector(C1_I, 1, num_bp);
    free_lvector(C1_II, 1, num_bp);
    free_lvector(N_I, 1, num_bp);
    free_lvector(N_II, 1, num_bp);
    free_lvector(b_vecNum, 1, n_vec);
    free_lvector(e_vecNum, 1, n_vec);
    free_dvector(z_axis, 1, 3);
    free_lvector(tmp1, 1, n_base - 1);
    free_lvector(tmp2, 1, n_base + 1);
    free_lvector(idx, 1, n_base);
    free_lvector(base_bn, 1, lbn);
    free_lvector(base_en, 1, lbn);
    free_lvector(bp1_idx, 1, num_bp1);
    free_dvector(bp1_normal, 1, 3);
    free_dvector(bp1x, 1, 3);
    free_dvector(bp1y, 1, 3);
    free_dvector(bp1z, 1, 3);
    free_dvector(pnts_dist, 1, num_bp1);
    free_dvector(hinge, 1, 3);
    free_dvector(g, 1, n_vec);
    free_dvector(org_xyz, 1, 3);
    free_dmatrix(vec_xyz, 1, n_vec, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
    free_dmatrix(rotmat_T, 1, 3, 1, 3);
    free_dmatrix(xyz2, 1, num_atom, 1, 3);
    free_dmatrix(bp1_xyz, 1, num_bp1, 1, 3);
    free_dmatrix(bp1xyz, 1, 3, 1, 3);
    free_dmatrix(glbxyz, 1, 3, 1, 3);
    free_dmatrix(b_xy, 1, n_vec, 1, 2);
    free_dmatrix(e_xy, 1, n_vec, 1, 2);
    free_dmatrix(dxy, 1, n_vec, 1, 2);
    free_dmatrix(t2x2, 1, 2, 1, 2);
    free_dmatrix(dxy_T, 1, 2, 1, n_vec);
    free_dmatrix(inv_t2x2, 1, 2, 1, 2);
}

void wrt_hel_dat(long num_atom, char **atom_id, char *base_type,
                 char *strand, long *base_num, double **xyz,
                 long *base, long n_base, char *filstr)
     /* WRT_HEL_DAT write the helical structural data in PDB format */
{
    char **b_asym, *b_btype, *b_strand, filnam[BUF512];
    long i, j, *b_bnum;
    double **b_xyz;

    /* First: Write out the whole structure */
    strcpy(filnam, filstr);
    wrtpdb(num_atom, atom_id, base_type, strand, base_num, xyz, strcat(filnam, "all"));

    /* Second: Write out the structure with ONLY base */
    b_asym = cmatrix(1, n_base, 1, 4);
    b_btype = cvector(1, n_base);
    b_strand = cvector(1, n_base);
    b_bnum = lvector(1, n_base);
    b_xyz = dmatrix(1, n_base, 1, 3);

    for (i = 1; i <= n_base; i++) {
        strcpy(b_asym[i], atom_id[base[i]]);
        b_btype[i] = base_type[base[i]];
        b_strand[i] = strand[base[i]];
        b_bnum[i] = base_num[base[i]];
        for (j = 1; j <= 3; j++)
            b_xyz[i][j] = xyz[base[i]][j];
    }

    strcpy(filnam, filstr);
    wrtpdb(n_base, b_asym, b_btype, b_strand, b_bnum, b_xyz, strcat(filnam, "bp"));

    free_cmatrix(b_asym, 1, n_base, 1, 4);
    free_cvector(b_btype, 1, n_base);
    free_cvector(b_strand, 1, n_base);
    free_lvector(b_bnum, 1, n_base);
    free_dmatrix(b_xyz, 1, n_base, 1, 3);
}

void main_chain(long *mc, long nmc, long *chi, long num_bp, double **xyz,
                char **bp_seq, FILE * fo)

     /* MAIN_CHAIN get the torsion angles of the sugar-phosphate backbone
        and glycosyl angle chi */
{
    long ip1;  /* Indicator for P atom */
    long nmcta;  /* Number of main-chain torsion angles */
    double *mcta1, *mcta2, *mcta_total, **mcta_matrix, *ep_ze;
    long num2;  /* Beginning number of strand II */
    long i, j, n, k, nt, *tor_idx;
    double **tmc_xyz;
    double *chi1, *chi2;
    char bstr[BUF512];

    /* Get the main chain torsion angles, all in 5'-->3' direction of I */
    ip1 = (nmc / 2) % 2;  /* 0 if first atom is P, 1 for O5' */
    nt = 6 * num_bp;  /* Number of total torsion angles */
    nmcta = nt - 3 - ip1;

    mcta1 = dvector(1, nmcta);
    mcta2 = dvector(1, nmcta);
    mcta_total = dvector(1, nt);
    mcta_matrix = dmatrix(1, num_bp, 1, 6);

    tor_idx = lvector(1, 4);
    tmc_xyz = dmatrix(1, 4, 1, 3);  /* Four sets of xyz coordinates */

    for (i = 0; i <= 3; i++)
        tor_idx[i + 1] = i;
    num2 = nmc / 2;

    for (i = 1; i <= nmcta; i++) {
        for (j = 1; j <= 4; j++)  /* For strand I */
            for (k = 1; k <= 3; k++)
                tmc_xyz[j][k] = xyz[mc[i + tor_idx[j]]][k];
        mcta1[i] = torsion(tmc_xyz);

        for (j = 1; j <= 4; j++)  /* For strand II */
            for (k = 1; k <= 3; k++)
                tmc_xyz[j][k] = xyz[mc[num2 + i + tor_idx[j]]][k];
        mcta2[i] = torsion(tmc_xyz);
    }

    /* Get the chi torsion angles */
    chi1 = dvector(1, num_bp);
    chi2 = dvector(1, num_bp);
    num2 = 4 * num_bp;

    for (i = 1; i <= num_bp; i++) {
        n = (i - 1) * 4 + 1;
        for (j = 1; j <= 4; j++)
            for (k = 1; k <= 3; k++)
                tmc_xyz[j][k] = xyz[chi[n + tor_idx[j]]][k];
        chi1[i] = torsion(tmc_xyz);

        for (j = 1; j <= 4; j++)
            for (k = 1; k <= 3; k++)
                tmc_xyz[j][k] = xyz[chi[num2 + n + tor_idx[j]]][k];
        chi2[i] = torsion(tmc_xyz);
    }

    if (ip1 == 1)
        mcta_total[2] = 0.0;  /* First beta */
    mcta_total[1] = 0.0;  /* First alpha */
    mcta_total[nt - 1] = 0.0;  /* Last epsilon */
    mcta_total[nt] = 0.0;  /* Last zeta */

    ep_ze = dvector(1, num_bp);

    /* For strand I */
    for (i = 1; i <= nmcta; i++)
        mcta_total[i + ip1 + 1] = mcta1[i];
    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 6; j++)
            mcta_matrix[i][j] = mcta_total[(i - 1) * 6 + j];
        ep_ze[i] = mcta_matrix[i][5] - mcta_matrix[i][6];
    }

    /* Print out the main chain torsion angle and chi */
    prt_sep(fo, '*', 76);
    fprintf(fo, "Main chain torsion angles: \n\n");
    fprintf(fo, "Note: alpha:   O3'-P-O5'-C5'\n");
    fprintf(fo, "      beta:    P-O5'-C5'-C4'\n");
    fprintf(fo, "      gamma:   O5'-C5'-C4'-C3'\n");
    fprintf(fo, "      delta:   C5'-C4'-C3'-O3'\n");
    fprintf(fo, "      epsilon: C4'-C3'-O3'-P\n");
    fprintf(fo, "      zeta:    C3'-O3'-P-O5'\n\n");
    fprintf(fo, "      chi for pyrimidines: O4'-C1'-N1-C2\n");
    fprintf(fo, "          chi for purines: O4'-C1'-N9-C4\n\n");

    fprintf(fo, "Strand I\n");
    fprintf(fo, " base    alpha    beta   gamma   delta");
    fprintf(fo, " epsilon    zeta   ep-ze    chi\n");

    strcpy(bstr, "     -- ");

    i = 1;  /* For the first set of torsion angles: No alpha OR beta */
    fprintf(fo, "%2ld%3c ", i, bp_seq[i][1]);
    if (ip1 == 1)
        fprintf(fo, "%s%s", bstr, bstr);
    else
        fprintf(fo, "%s", bstr);
    for (j = ip1 + 2; j <= 6; j++)
        fprintf(fo, "%8.2f", mcta_matrix[i][j]);
    fprintf(fo, "%8.2f%8.2f\n", ep_ze[i], chi1[i]);

    for (i = 2; i <= num_bp - 1; i++) {  /* Middle sets */
        fprintf(fo, "%2ld%3c ", i, bp_seq[i][1]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%8.2f", mcta_matrix[i][j]);
        fprintf(fo, "%8.2f%8.2f\n", ep_ze[i], chi1[i]);
    }

    i = num_bp;  /* The last set: No epsilon & zeta */
    fprintf(fo, "%2ld%3c ", i, bp_seq[i][1]);
    for (j = 1; j <= 4; j++)
        fprintf(fo, "%8.2f", mcta_matrix[i][j]);
    fprintf(fo, "%s%s%s%8.2f\n", bstr, bstr, bstr, chi1[i]);
    fprintf(fo, "       ");
    prt_sep(fo, '~', 63);

    /* For strand II */
    for (i = 1; i <= nmcta; i++)
        mcta_total[i + ip1 + 1] = mcta2[i];
    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 6; j++)
            mcta_matrix[i][j] = mcta_total[(i - 1) * 6 + j];
        ep_ze[i] = mcta_matrix[i][5] - mcta_matrix[i][6];
    }

    fprintf(fo, "\nStrand II\n");
    fprintf(fo, " base    alpha    beta   gamma   delta");
    fprintf(fo, " epsilon    zeta   ep-ze    chi\n");

    /* Note the reversed order! */
    i = num_bp;  /* 1st set of torsion angles: No epsilon & zeta (3'-->5'!!!) */
    fprintf(fo, "%2d%3c ", 1, bp_seq[1][2]);
    for (j = 1; j <= 4; j++)
        fprintf(fo, "%8.2f", mcta_matrix[i][j]);
    fprintf(fo, "%s%s%s%8.2f\n", bstr, bstr, bstr, chi2[i]);

    for (i = num_bp - 1; i >= 2; i--) {  /* Middle sets */
        j = num_bp + 1 - i;
        fprintf(fo, "%2ld%3c ", j, bp_seq[j][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%8.2f", mcta_matrix[i][j]);
        fprintf(fo, "%8.2f%8.2f\n", ep_ze[i], chi2[i]);
    }

    i = 1;  /* For the LAST set of torsion angles: No alpha OR beta */
    fprintf(fo, "%2ld%3c ", num_bp, bp_seq[num_bp][2]);
    if (ip1 == 1)
        fprintf(fo, "%s%s", bstr, bstr);
    else
        fprintf(fo, "%s", bstr);
    for (j = ip1 + 2; j <= 6; j++)
        fprintf(fo, "%8.2f", mcta_matrix[i][j]);
    fprintf(fo, "%8.2f%8.2f\n", ep_ze[i], chi2[i]);
    fprintf(fo, "       ");
    prt_sep(fo, '~', 63);

    fprintf(fo, "\n");

    free_dvector(mcta1, 1, nmcta);
    free_dvector(mcta2, 1, nmcta);
    free_dvector(mcta_total, 1, nt);
    free_lvector(tor_idx, 1, 4);
    free_dvector(chi1, 1, num_bp);
    free_dvector(chi2, 1, num_bp);
    free_dvector(ep_ze, 1, num_bp);
    free_dmatrix(mcta_matrix, 1, num_bp, 1, 6);
    free_dmatrix(tmc_xyz, 1, 4, 1, 3);
}

/* Get the conformational parameters of the sugar ring. */
void sugar_ana(long *sugar, long num_bp, double **xyz, char **bp_seq, FILE * fo)
{
    double **tor1, **tor2, **ang1, **ang2;
    double **abcd_xyz, **abc_xyz;
    long i, j, k, m, num2, *idx, *tor_j;

    double Pconst, dtmp, *P_1, *tm_1, *P_2, *tm_2;

    double *ave_tor, *std_tor, ave_P, std_P, ave_tm, std_tm;

    /* Added June-9-1998  */
    long *sp1_idx, *sp2_idx;  /* Index of sugar packering modes */
    char *sugar_pucker[10];
    sugar_pucker[0] = "C3'-endo";
    sugar_pucker[1] = "C4'-exo ";
    sugar_pucker[2] = "O4'-endo";
    sugar_pucker[3] = "C1'-exo ";
    sugar_pucker[4] = "C2'-endo";
    sugar_pucker[5] = "C3'-exo ";
    sugar_pucker[6] = "C4'-endo";
    sugar_pucker[7] = "O4'-exo ";
    sugar_pucker[8] = "C1'-endo";
    sugar_pucker[9] = "C2'-exo ";
    /* End of addition */

    /* Sugar torsion angles and valence angles */
    tor1 = dmatrix(1, num_bp, 1, 5);
    tor2 = dmatrix(1, num_bp, 1, 5);
    ang1 = dmatrix(1, num_bp, 1, 5);
    ang2 = dmatrix(1, num_bp, 1, 5);
    abcd_xyz = dmatrix(1, 4, 1, 3);
    abc_xyz = dmatrix(1, 3, 1, 3);

    idx = lvector(1, 5);
    tor_j = lvector(1, 5);

    sp1_idx = lvector(1, num_bp);
    sp2_idx = lvector(1, num_bp);

    num2 = 5 * num_bp;
    for (i = 1; i <= 5; i++)
        idx[i] = i;

    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 5; j++) {
            cirshift(tor_j, idx, 5, 1 - j);
            for (m = 1; m <= 4; m++)
                for (k = 1; k <= 3; k++)
                    abcd_xyz[m][k] = xyz[sugar[tor_j[m]]][k];
            tor1[i][j] = torsion(abcd_xyz);

            for (m = 1; m <= 4; m++)
                for (k = 1; k <= 3; k++)
                    abcd_xyz[m][k] = xyz[sugar[num2 + tor_j[m]]][k];
            tor2[i][j] = torsion(abcd_xyz);

            for (m = 1; m <= 3; m++)
                for (k = 1; k <= 3; k++)
                    abc_xyz[m][k] = xyz[sugar[tor_j[m]]][k];
            ang1[i][j] = val_ang(abc_xyz);

            for (m = 1; m <= 3; m++)
                for (k = 1; k <= 3; k++)
                    abc_xyz[m][k] = xyz[sugar[num2 + tor_j[m]]][k];
            ang2[i][j] = val_ang(abc_xyz);
        }
        for (j = 1; j <= 5; j++)
            idx[j] = 5 + idx[j];
    }

    /* Get pseudorotation phase angle P and the amplitude of puckering
       1. tan(P)=((v4+v1)-(v3+v0))/(2*v2*(sin(36)+sin(72)))
       2. tm=v2/cos(P). Wrong in Saenger's book, pp.20
       3. Here index is 1-offset instead of 0-offset as in the eqs. */
    Pconst = sin(PI / 5.0) + sin(PI / 2.5);
    P_1 = dvector(1, num_bp);
    tm_1 = dvector(1, num_bp);
    P_2 = dvector(1, num_bp);
    tm_2 = dvector(1, num_bp);

    for (i = 1; i <= num_bp; i++) {
        dtmp = (tor1[i][5] + tor1[i][2]) - (tor1[i][4] + tor1[i][1]);
        P_1[i] = atan2(dtmp, 2.0 * tor1[i][3] * Pconst);
        tm_1[i] = tor1[i][3] / cos(P_1[i]);
        P_1[i] = rad2deg(P_1[i]);
        if (P_1[i] < 0.0)
            P_1[i] += 360.0;
        sp1_idx[i] = (long) (P_1[i] / 36.0);  /* New adding */

        dtmp = (tor2[i][5] + tor2[i][2]) - (tor2[i][4] + tor2[i][1]);
        P_2[i] = atan2(dtmp, 2.0 * tor2[i][3] * Pconst);
        tm_2[i] = tor2[i][3] / cos(P_2[i]);
        P_2[i] = rad2deg(P_2[i]);
        if (P_2[i] < 0.0)
            P_2[i] += 360.0;
        sp2_idx[i] = (long) (P_2[i] / 36.0);  /* New adding */
    }

    ave_tor = dvector(1, 5);
    std_tor = dvector(1, 5);

    /* Print out sugar torsion angles, P, tm and puckering mode */
    prt_sep(fo, '*', 76);
    fprintf(fo, "Sugar conformation parameters: \n\n");
    fprintf(fo, "Note: v0: C4'-O4'-C1'-C2'\n");
    fprintf(fo, "      v1: O4'-C1'-C2'-C3'\n");
    fprintf(fo, "      v2: C1'-C2'-C3'-C4'\n");
    fprintf(fo, "      v3: C2'-C3'-C4'-O4'\n");
    fprintf(fo, "      v4: C3'-C4'-O4'-C1'\n\n");
    fprintf(fo, "      P:  The phase angle of pseudorotation of the sugar ring\n");
    fprintf(fo, "      tm: The amplitude of pseudorotation of the sugar ring\n\n");

    fprintf(fo, "Strand I:\n");
    fprintf(fo,
            " base      v0      v1      v2      v3      v4      P       tm   Puckering\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fo, "%2ld%3c ", i, bp_seq[i][1]);
        for (j = 1; j <= 5; j++)
            fprintf(fo, "%8.2f", tor1[i][j]);
        fprintf(fo, "%8.2f%8.2f%11s\n", P_1[i], tm_1[i], sugar_pucker[sp1_idx[i]]);
    }
    fprintf(fo, "       ");
    prt_sep(fo, '~', 66);

    mean_dmatrix(ave_tor, tor1, num_bp, 5);
    std_dmatrix(std_tor, tor1, num_bp, 5);
    ave_P = mean_dvector(P_1, num_bp);
    std_P = std_dvector(P_1, num_bp);
    ave_tm = mean_dvector(tm_1, num_bp);
    std_tm = std_dvector(tm_1, num_bp);

    fprintf(fo, "  ave.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%8.2f", ave_tor[i]);
    fprintf(fo, "%8.2f", ave_P);
    fprintf(fo, "%8.2f\n", ave_tm);

    fprintf(fo, "  s.d.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%8.2f", std_tor[i]);
    fprintf(fo, "%8.2f", std_P);
    fprintf(fo, "%8.2f\n\n", std_tm);

    fprintf(fo, "Strand II:\n");
    fprintf(fo,
            " base      v0      v1      v2      v3      v4      P       tm   Puckering\n");
    for (i = num_bp; i >= 1; i--) {
        j = num_bp + 1 - i;
        fprintf(fo, "%2ld%3c ", j, bp_seq[j][2]);
        for (j = 1; j <= 5; j++)
            fprintf(fo, "%8.2f", tor2[i][j]);
        fprintf(fo, "%8.2f%8.2f%11s\n", P_2[i], tm_2[i], sugar_pucker[sp2_idx[i]]);
    }
    fprintf(fo, "       ");
    prt_sep(fo, '~', 66);

    mean_dmatrix(ave_tor, tor2, num_bp, 5);
    std_dmatrix(std_tor, tor2, num_bp, 5);
    ave_P = mean_dvector(P_2, num_bp);
    std_P = std_dvector(P_2, num_bp);
    ave_tm = mean_dvector(tm_2, num_bp);
    std_tm = std_dvector(tm_2, num_bp);

    fprintf(fo, "  ave.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%8.2f", ave_tor[i]);
    fprintf(fo, "%8.2f", ave_P);
    fprintf(fo, "%8.2f\n", ave_tm);

    fprintf(fo, "  s.d.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%8.2f", std_tor[i]);
    fprintf(fo, "%8.2f", std_P);
    fprintf(fo, "%8.2f\n\n", std_tm);

    /* Print out sugar valence angles */
    prt_sep(fo, '-', 76);
    fprintf(fo, "Sugar bond angles for strand I:\n");
    fprintf(fo, " base    C4-O4-C1  O4-C1-C2  C1-C2-C3  C2-C3-C4  C3-C4-O4\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fo, "%2ld%3c ", i, bp_seq[i][1]);
        for (j = 1; j <= 5; j++)
            fprintf(fo, "%10.2f", ang1[i][j]);
        fprintf(fo, "\n");
    }
    fprintf(fo, "       ");
    prt_sep(fo, '~', 49);

    mean_dmatrix(ave_tor, ang1, num_bp, 5);
    std_dmatrix(std_tor, ang1, num_bp, 5);

    fprintf(fo, "  ave.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%10.2f", ave_tor[i]);
    fprintf(fo, "\n");
    fprintf(fo, "  s.d.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%10.2f", std_tor[i]);
    fprintf(fo, "\n");

    fprintf(fo, "\nSugar bond angles for strand II:\n");
    fprintf(fo, " base    C4-O4-C1  O4-C1-C2  C1-C2-C3  C2-C3-C4  C3-C4-O4\n");
    for (i = num_bp; i >= 1; i--) {
        j = num_bp + 1 - i;
        fprintf(fo, "%2ld%3c ", j, bp_seq[j][2]);
        for (j = 1; j <= 5; j++)
            fprintf(fo, "%10.2f", ang2[i][j]);
        fprintf(fo, "\n");
    }
    fprintf(fo, "       ");
    prt_sep(fo, '~', 49);

    mean_dmatrix(ave_tor, ang2, num_bp, 5);
    std_dmatrix(std_tor, ang2, num_bp, 5);

    fprintf(fo, "  ave.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%10.2f", ave_tor[i]);
    fprintf(fo, "\n");
    fprintf(fo, "  s.d.");
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%10.2f", std_tor[i]);
    fprintf(fo, "\n\n");

    free_lvector(idx, 1, 5);
    free_lvector(tor_j, 1, 5);
    free_dvector(P_1, 1, num_bp);
    free_dvector(tm_1, 1, num_bp);
    free_dvector(P_2, 1, num_bp);
    free_dvector(tm_2, 1, num_bp);
    free_dvector(ave_tor, 1, 5);
    free_dvector(std_tor, 1, 5);
    free_dmatrix(tor1, 1, num_bp, 1, 5);
    free_dmatrix(tor2, 1, num_bp, 1, 5);
    free_dmatrix(ang1, 1, num_bp, 1, 5);
    free_dmatrix(ang2, 1, num_bp, 1, 5);
    free_dmatrix(abcd_xyz, 1, 4, 1, 3);
    free_dmatrix(abc_xyz, 1, 3, 1, 3);
    free_lvector(sp1_idx, 1, num_bp);
    free_lvector(sp2_idx, 1, num_bp);
}

void polar_coord(long *alist, long midnum, char *asym, double **xyz,
                 char **bp_seq, FILE * fo)
     /* POLAR_COORD get polar coordinates and relevant parameters */
{
    long i, j;
    double **xyz1, **xyz2, *r1, *phi1, *dphi1, *r2, *phi2, *dphi2;
    double *dist_xyz1, *dist_xy1, *pitch1, *dist_xyz2, *dist_xy2, *pitch2;
    double *dxyz, *dz1, *dz2, *ave6, *std6;
    char bstr[BUF512];

    /* Get xyz coordinates for each strand */
    xyz1 = dmatrix(1, midnum, 1, 3);
    xyz2 = dmatrix(1, midnum, 1, 3);
    for (i = 1; i <= midnum; i++)
        for (j = 1; j <= 3; j++) {
            xyz1[i][j] = xyz[alist[i]][j];
            xyz2[i][j] = xyz[alist[midnum + i]][j];
        }

    /* (1). Polar coordinates (r, phi and z), phi in [-180,+180]
       Note: These parameters are helical axis dependent */
    r1 = dvector(1, midnum);
    r2 = dvector(1, midnum);
    phi1 = dvector(1, midnum);
    phi2 = dvector(1, midnum);
    dphi1 = dvector(1, midnum - 1);
    dphi2 = dvector(1, midnum - 1);

    for (i = 1; i <= midnum; i++) {
        r1[i] = dist_ab(xyz1[i][1], xyz1[i][2]);
        phi1[i] = rad2deg(atan2(xyz1[i][2], xyz1[i][1]));
        r2[i] = dist_ab(xyz2[i][1], xyz2[i][2]);
        phi2[i] = rad2deg(atan2(xyz2[i][2], xyz2[i][1]));
    }
    for (i = 1; i <= midnum - 1; i++) {
        dphi1[i] = phi1[i + 1] - phi1[i];
        if (dphi1[i] < -180)
            dphi1[i] += 360.0;
        if (dphi1[i] > +180)
            dphi1[i] -= 360.0;
        dphi2[i] = phi2[i] - phi2[i + 1];  /* Note the change of order */
        if (dphi2[i] < -180)
            dphi2[i] += 360.0;
        if (dphi2[i] > +180)
            dphi2[i] -= 360.0;
    }

    /* (2). Get distance between atoms, its projections etc */
    dist_xyz1 = dvector(1, midnum - 1);
    dist_xy1 = dvector(1, midnum - 1);
    pitch1 = dvector(1, midnum - 1);
    dist_xyz2 = dvector(1, midnum - 1);
    dist_xy2 = dvector(1, midnum - 1);
    pitch2 = dvector(1, midnum - 1);
    dz1 = dvector(1, midnum - 1);
    dz2 = dvector(1, midnum - 1);
    dxyz = dvector(1, 3);

    for (i = 1; i <= midnum - 1; i++) {
        for (j = 1; j <= 3; j++)
            dxyz[j] = xyz1[i + 1][j] - xyz1[i][j];
        dz1[i] = dxyz[3];
        dist_xyz1[i] = mag_dvector(dxyz, 3);
        dist_xy1[i] = dist_ab(dxyz[1], dxyz[2]);
        pitch1[i] = rad2deg(asin(dxyz[3] / dist_xyz1[i]));

        for (j = 1; j <= 3; j++)
            dxyz[j] = xyz2[i][j] - xyz2[i + 1][j];
        /* Note the change of order */
        dz2[i] = dxyz[3];
        dist_xyz2[i] = mag_dvector(dxyz, 3);
        dist_xy2[i] = dist_ab(dxyz[1], dxyz[2]);
        pitch2[i] = rad2deg(asin(dxyz[3] / dist_xyz2[i]));
    }

    strcpy(bstr, "     -- ");

    /* (3). Print out polar coordinates and distances */
    /* Strand I */
    fprintf(fo, "%s atoms on strand I:", asym);
    fprintf(fo, "\n  base        r     phi       z    d_phi");
    fprintf(fo, "     dz      dxy    dxyz   pitch\n");
    for (i = 1; i <= midnum - 1; i++) {
        fprintf(fo, "%2ld %3c  ", i, bp_seq[i][1]);
        fprintf(fo, "%8.2f%8.2f%8.2f%8.2f", r1[i], phi1[i], xyz1[i][3], dphi1[i]);
        fprintf(fo, "%8.2f%8.2f%8.2f%8.2f\n", dz1[i], dist_xy1[i], dist_xyz1[i],
                pitch1[i]);
    }
    i = midnum;
    fprintf(fo, "%2ld %3c  ", i, bp_seq[i][1]);
    fprintf(fo, "%8.2f%8.2f%8.2f", r1[i], phi1[i], xyz1[i][3]);
    for (i = 1; i <= 5; i++)
        fprintf(fo, "%s", bstr);
    fprintf(fo, "\n");

    ave6 = dvector(1, 6);
    std6 = dvector(1, 6);

    ave6[1] = mean_dvector(r1, midnum);
    std6[1] = std_dvector(r1, midnum);

    ave6[2] = mean_dvector(dphi1, midnum - 1);
    std6[2] = std_dvector(dphi1, midnum - 1);
    ave6[3] = mean_dvector(dz1, midnum - 1);
    std6[3] = std_dvector(dz1, midnum - 1);
    ave6[4] = mean_dvector(dist_xy1, midnum - 1);
    std6[4] = std_dvector(dist_xy1, midnum - 1);
    ave6[5] = mean_dvector(dist_xyz1, midnum - 1);
    std6[5] = std_dvector(dist_xyz1, midnum - 1);
    ave6[6] = mean_dvector(pitch1, midnum - 1);
    std6[6] = std_dvector(pitch1, midnum - 1);

    fprintf(fo, "        ");
    prt_sep(fo, '~', 64);
    fprintf(fo, "    ave.%8.2f    --      --  ", ave6[1]);
    for (i = 2; i <= 6; i++)
        fprintf(fo, "%8.2f", ave6[i]);
    fprintf(fo, "\n");
    fprintf(fo, "    s.d.%8.2f    --      --  ", std6[1]);
    for (i = 2; i <= 6; i++)
        fprintf(fo, "%8.2f", std6[i]);
    fprintf(fo, "\n");

    /* Strand II */
    fprintf(fo, "\n%s atoms on strand II:", asym);
    fprintf(fo, "\n  base        r     phi       z    d_phi");
    fprintf(fo, "     dz      dxy    dxyz   pitch\n");
    for (i = midnum - 1; i >= 1; i--) {
        j = midnum - i;
        fprintf(fo, "%2ld %3c  ", j, bp_seq[j][2]);
        fprintf(fo, "%8.2f%8.2f%8.2f%8.2f", r2[i + 1], phi2[i + 1], xyz2[i + 1][3],
                dphi2[i]);
        fprintf(fo, "%8.2f%8.2f%8.2f%8.2f\n", dz2[i], dist_xy2[i], dist_xyz2[i],
                pitch2[i]);
    }
    i = 1;
    fprintf(fo, "%2ld %3c  ", midnum, bp_seq[midnum][2]);
    fprintf(fo, "%8.2f%8.2f%8.2f", r2[i], phi2[i], xyz2[i][3]);

    for (i = 1; i <= 5; i++)
        fprintf(fo, "%s", bstr);
    fprintf(fo, "\n");

    ave6[1] = mean_dvector(r2, midnum);
    std6[1] = std_dvector(r2, midnum);

    ave6[2] = mean_dvector(dphi2, midnum - 1);
    std6[2] = std_dvector(dphi2, midnum - 1);
    ave6[3] = mean_dvector(dz2, midnum - 1);
    std6[3] = std_dvector(dz2, midnum - 1);
    ave6[4] = mean_dvector(dist_xy2, midnum - 1);
    std6[4] = std_dvector(dist_xy2, midnum - 1);
    ave6[5] = mean_dvector(dist_xyz2, midnum - 1);
    std6[5] = std_dvector(dist_xyz2, midnum - 1);
    ave6[6] = mean_dvector(pitch2, midnum - 1);
    std6[6] = std_dvector(pitch2, midnum - 1);

    fprintf(fo, "        ");
    prt_sep(fo, '~', 64);
    fprintf(fo, "    ave.%8.2f    --      --  ", ave6[1]);
    for (i = 2; i <= 6; i++)
        fprintf(fo, "%8.2f", ave6[i]);
    fprintf(fo, "\n");
    fprintf(fo, "    s.d.%8.2f    --      --  ", std6[1]);
    for (i = 2; i <= 6; i++)
        fprintf(fo, "%8.2f", std6[i]);
    fprintf(fo, "\n\n");

    free_dmatrix(xyz1, 1, midnum, 1, 3);
    free_dmatrix(xyz2, 1, midnum, 1, 3);
    free_dvector(r1, 1, midnum);
    free_dvector(r2, 1, midnum);
    free_dvector(phi1, 1, midnum);
    free_dvector(phi2, 1, midnum);
    free_dvector(dphi1, 1, midnum - 1);
    free_dvector(dphi2, 1, midnum - 1);
    free_dvector(dist_xyz1, 1, midnum - 1);
    free_dvector(dist_xy1, 1, midnum - 1);
    free_dvector(pitch1, 1, midnum - 1);
    free_dvector(dist_xyz2, 1, midnum - 1);
    free_dvector(dist_xy2, 1, midnum - 1);
    free_dvector(pitch2, 1, midnum - 1);
    free_dvector(dz1, 1, midnum - 1);
    free_dvector(dz2, 1, midnum - 1);
    free_dvector(dxyz, 1, 3);
    free_dvector(ave6, 1, 6);
    free_dvector(std6, 1, 6);
}

void get_lambda(long *C1, long *N, long *C68, long num_bp, double **xyz,
                char **bp_seq, FILE * fo)
     /* GET_LAMBDA get the angle between C1'- C1' and the glycosidic bond
        and the distances between C1'-C1', C6-C8 for each bp */
{
    long i, j, k, n;
    double **C1xyz, **Nxyz, **C68xyz;
    double **C1_C1, **diff_C68, **N_C1_1, **N_C1_2;
    double *lambda1, *lambda2, *tmp_vec;
    double *C1_dist, *C68_dist;

    n = 2 * num_bp;
    C1xyz = dmatrix(1, n, 1, 3);
    Nxyz = dmatrix(1, n, 1, 3);
    C68xyz = dmatrix(1, n, 1, 3);

    /* Define the vectors */
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= 3; j++) {
            C1xyz[i][j] = xyz[C1[i]][j];
            Nxyz[i][j] = xyz[N[i]][j];
            C68xyz[i][j] = xyz[C68[i]][j];
        }
    }

    C1_C1 = dmatrix(1, num_bp, 1, 3);
    diff_C68 = dmatrix(1, num_bp, 1, 3);
    N_C1_1 = dmatrix(1, num_bp, 1, 3);
    N_C1_2 = dmatrix(1, num_bp, 1, 3);

    for (i = 1; i <= num_bp; i++) {
        k = n - i + 1;  /* Matching pair: i <===> 2*num_bp-i+1 */
        for (j = 1; j <= 3; j++) {
            C1_C1[i][j] = C1xyz[i][j] - C1xyz[k][j];
            diff_C68[i][j] = C68xyz[i][j] - C68xyz[k][j];
            N_C1_1[i][j] = Nxyz[i][j] - C1xyz[i][j];
            N_C1_2[i][j] = Nxyz[k][j] - C1xyz[k][j];
        }
    }

    lambda1 = dvector(1, num_bp);
    lambda2 = dvector(1, num_bp);
    tmp_vec = dvector(1, 3);

    /* Get the lambda angles */
    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++)
            tmp_vec[j] = -C1_C1[i][j];
        lambda1[i] = magang(tmp_vec, N_C1_1[i]);
        lambda2[i] = magang(N_C1_2[i], C1_C1[i]);
    }

    /* Get the C1'-C1' and C6-C8 distances */
    C1_dist = dvector(1, num_bp);
    C68_dist = dvector(1, num_bp);

    for (i = 1; i <= num_bp; i++) {
        C1_dist[i] = mag_dvector(C1_C1[i], 3);
        C68_dist[i] = mag_dvector(diff_C68[i], 3);
    }

    /* Print out the lambda angles */
    prt_sep(fo, '*', 76);
    fprintf(fo, "Lambda:    the angles between C1'-N1 or C1'-N9 glycosidic bonds\n");
    fprintf(fo, "           and the base-pair C1'-C1' line;\n");
    fprintf(fo, "Distances: between C1'-C1' and RC8-YC6 for each base-pair");
    fprintf(fo, "\n\n   bp    Lambda(I) Lambda(II) C1'-C1'   RC8-YC6\n");
    for (i = 1; i <= num_bp; i++) {
        fprintf(fo, "%2ld %c-%c", i, bp_seq[i][1], bp_seq[i][2]);
        fprintf(fo, "%10.2f%10.2f", lambda1[i], lambda2[i]);
        fprintf(fo, "%10.2f%10.2f\n", C1_dist[i], C68_dist[i]);
    }

    fprintf(fo, "        ");
    prt_sep(fo, '~', 38);
    fprintf(fo, "  ave.");
    fprintf(fo, "%10.2f", mean_dvector(lambda1, num_bp));
    fprintf(fo, "%10.2f", mean_dvector(lambda2, num_bp));
    fprintf(fo, "%10.2f", mean_dvector(C1_dist, num_bp));
    fprintf(fo, "%10.2f\n", mean_dvector(C68_dist, num_bp));

    fprintf(fo, "  s.d.");
    fprintf(fo, "%10.2f", std_dvector(lambda1, num_bp));
    fprintf(fo, "%10.2f", std_dvector(lambda2, num_bp));
    fprintf(fo, "%10.2f", std_dvector(C1_dist, num_bp));
    fprintf(fo, "%10.2f\n\n", std_dvector(C68_dist, num_bp));

    free_dvector(lambda1, 1, num_bp);
    free_dvector(lambda2, 1, num_bp);
    free_dvector(tmp_vec, 1, 3);
    free_dvector(C1_dist, 1, num_bp);
    free_dvector(C68_dist, 1, num_bp);
    free_dmatrix(C1xyz, 1, n, 1, 3);
    free_dmatrix(Nxyz, 1, n, 1, 3);
    free_dmatrix(C68xyz, 1, n, 1, 3);
    free_dmatrix(C1_C1, 1, num_bp, 1, 3);
    free_dmatrix(diff_C68, 1, num_bp, 1, 3);
    free_dmatrix(N_C1_1, 1, num_bp, 1, 3);
    free_dmatrix(N_C1_2, 1, num_bp, 1, 3);
}

void groove_width(long *P, long midnum, long *O4, double **xyz, long num_bp,
                  char **bp_seq, FILE * fo)
     /* GROOVE_WIDTH get the reduced P-P(O4'-O4') distances */
{
    long i, j, k;
    double **xyz1, **xyz2, **adist, *tmp_vec;
    long ncircle, inum, icol, ist, max_col = 12;

    prt_sep(fo, '*', 76);
    fprintf(fo, "Groove width parameters\n\n");
    fprintf(fo, "Note: Strand I up to down in 5' to 3' direction\n");
    fprintf(fo, "      Strand II left to right in 3' to 5' direction\n");

    /* Part 1: Phosphorus-Phosphorus distances */

    /* Get xyz coordinates for each strand */
    xyz1 = dmatrix(1, num_bp, 1, 3);
    xyz2 = dmatrix(1, num_bp, 1, 3);

    for (i = 1; i <= midnum; i++) {
        k = 2 * midnum - i + 1;  /* Matching pair: i <===> 2*midnum-i+1 */
        for (j = 1; j <= 3; j++) {
            xyz1[i][j] = xyz[P[i]][j];
            xyz2[i][j] = xyz[P[k]][j];
        }
    }

    /* Get the distance */
    adist = dmatrix(1, num_bp, 1, num_bp);
    tmp_vec = dvector(1, 3);

    for (i = 1; i <= midnum; i++) {
        for (j = 1; j <= midnum; j++) {
            for (k = 1; k <= 3; k++)
                tmp_vec[k] = xyz1[i][k] - xyz2[j][k];
            adist[i][j] = mag_dvector(tmp_vec, 3) - 5.8;
        }
    }

    /* Print out the adist matrix */
    fprintf(fo, "\nPhosphorus-phosphorus distance reduced by 5.8A ");
    fprintf(fo, "for phosphate vdW radii\n\n");

    /* Do not Print P5' end associated distance if any */
    if (midnum == num_bp) {
        for (i = 1; i < midnum; i++)
            for (j = 1; j < midnum; j++)
                adist[i][j] = adist[i + 1][j];
        fprintf(fo, "(Note: 5' terminal P atoms NOT included!)\n\n");
    }

    /* max_col (12) dist per line maximum */
    ncircle = (long) ceil((num_bp - 1) / (double) max_col);

    for (inum = 1; inum <= ncircle; inum++) {
        if (inum == ncircle) {
            icol = (num_bp - 1) % max_col;  /* Remainder */
            if (!icol)
                icol = max_col;
        } else
            icol = max_col;

        ist = (inum - 1) * max_col;  /* Starting column */

        fprintf(fo, "      ");  /* The sequential of strand II P atoms */
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%6ld", i);
        fprintf(fo, "\n");

        fprintf(fo, "        ");  /* Base step */
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "  %c-%c ", bp_seq[i][2], bp_seq[i + 1][2]);
        fprintf(fo, "\n");

        for (i = 1; i <= num_bp - 1; i++) {
            fprintf(fo, "%2ld %c-%c  ", i, bp_seq[i][1], bp_seq[i + 1][1]);
            for (j = ist + 1; j <= ist + icol; j++)
                fprintf(fo, "%6.2f", adist[i][j]);
            fprintf(fo, "\n");
        }

        fprintf(fo, "\n");
        prt_sep(fo, '-', 80);
    }

    /* Part 2: O4'-O4' distances */
    /* Get xyz coordinates for each strand */
    for (i = 1; i <= num_bp; i++) {
        k = 2 * num_bp - i + 1;  /* Matching pair: i <===> 2*num_bp-i+1 */
        for (j = 1; j <= 3; j++) {
            xyz1[i][j] = xyz[O4[i]][j];
            xyz2[i][j] = xyz[O4[k]][j];
        }
    }

    /* Get the distance */
    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= num_bp; j++) {
            for (k = 1; k <= 3; k++)
                tmp_vec[k] = xyz1[i][k] - xyz2[j][k];
            adist[i][j] = mag_dvector(tmp_vec, 3) - 2.8;
        }
    }

    /* Print out the adist matrix */
    fprintf(fo, "O4'-O4' distance reduced by 2.8A ");
    fprintf(fo, "for oxygen vdW radii\n\n");

    /* max_col (12) dist per line maximum */
    ncircle = (long) ceil(num_bp / (double) max_col);

    for (inum = 1; inum <= ncircle; inum++) {
        if (inum == ncircle) {
            icol = num_bp % max_col;  /* Remainder */
            if (!icol)
                icol = max_col;
        } else
            icol = max_col;

        ist = (inum - 1) * max_col;  /* Starting column */

        fprintf(fo, "      ");  /* The sequential of strand II O4' atoms */
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%6ld", i);
        fprintf(fo, "\n");

        fprintf(fo, "      ");  /* Base */
        for (i = ist + 1; i <= ist + icol; i++)
            fprintf(fo, "%6c", bp_seq[i][2]);
        fprintf(fo, "\n");

        for (i = 1; i <= num_bp; i++) {
            fprintf(fo, "%2ld  %c   ", i, bp_seq[i][1]);
            for (j = ist + 1; j <= ist + icol; j++)
                fprintf(fo, "%6.2f", adist[i][j]);
            fprintf(fo, "\n");
        }

        fprintf(fo, "\n");
        prt_sep(fo, '-', 80);
    }

    free_dvector(tmp_vec, 1, 3);
    free_dmatrix(xyz1, 1, num_bp, 1, 3);
    free_dmatrix(xyz2, 1, num_bp, 1, 3);
    free_dmatrix(adist, 1, num_bp, 1, num_bp);
}

void get_frame(double **bx, double **by, double **bz, double **borg, long *xdir,
               long *RY_idx, long *C68, long num_bp, long *base, long n_base,
               long *base_num, char **atom_id, double **xyz, char **bp_seq, FILE * fo)
     /* GET_FRAME get the reference frame for each base & base-pair
        bx,by,bz--(x,y,z)-axes of the reference triads
        borg--xyz coordinates of the reference triads
        xdir--the +x-axis direction */
{
    long i, j, k, lbn;
    long *tmp1, *tmp2, *idx;
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
    idx = lvector(1, n_base);

    /* Good for handling bases with continuous sequential # */
    for (i = 1; i <= n_base; i++)
        idx[i] = base[i] + 100 * base_num[base[i]];

    i_diff(tmp1, idx, n_base);
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

    *xdir = sum_sign;

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

void get_bp_par(double *bp_par, double **xyz_rb, double *org_rb,
                double **xyz_lb, double *org_lb)
     /* GET_BP_PAR get the 6 local CEHS base-pair parameters */
{
    long i;
    double *hinge, buckleopening, phi;
    double propeller, opening, buckle, stretch = 0, stagger = 0, shear = 0;
    double *tv1, *tv2, *tv3;
    double **para_left, **para_right, **rotmat;
    double *mbpx, *mbpy, *mbpz;

    hinge = dvector(1, 3);
    tv1 = dvector(1, 3);
    tv2 = dvector(1, 3);
    tv3 = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        tv1[i] = xyz_rb[i][2];
        tv2[i] = xyz_lb[i][2];
    }

    /* Get the hinge by vector product of the y-axes of the two bases */
    cross_dvector(hinge, tv1, tv2);

    /* Get the buckleopening angle.
       It is a positive angle between the y-axes of the two base */
    buckleopening = magang(tv1, tv2);

    para_left = dmatrix(1, 3, 1, 3);
    para_right = dmatrix(1, 3, 1, 3);
    rotmat = dmatrix(1, 3, 1, 3);
    /* % We rotate the left base by -buckleopening/2
       %               right base by +buckleopening/2
       %   Note: The two bases are now co-planar */
    arbrot(rotmat, hinge, -0.5 * buckleopening);
    mul_dmatrix(para_left, rotmat, xyz_lb, 3, 3, 3, 3);

    arbrot(rotmat, hinge, 0.5 * buckleopening);
    mul_dmatrix(para_right, rotmat, xyz_rb, 3, 3, 3, 3);

    /* The propeller angle is that between the two x-axes of the two bases */
    for (i = 1; i <= 3; i++) {
        tv1[i] = para_right[i][1];
        tv2[i] = para_left[i][1];
        tv3[i] = para_left[i][2];
    }
    propeller = vec_ang(tv1, tv2, tv3);

    /* % We get the MBP
       %   Its y axis is the same as the two bases given above
       %   Its x axis is the mean of the two bases given above
       %   Its z axis can be got by the vector product of its x & y-axes */
    mbpx = dvector(1, 3);
    mbpy = dvector(1, 3);
    mbpz = dvector(1, 3);

    for (i = 1; i <= 3; i++)
        mbpy[i] = para_left[i][2];
    for (i = 1; i <= 3; i++)
        mbpx[i] = 0.5 * (para_left[i][1] + para_right[i][1]);
    norm_dvector(mbpx, mbpx, 3);
    cross_dvector(mbpz, mbpx, mbpy);

    /* We then need to get the phi angle which is defined by hinge & mbpx */
    phi = vec_ang(hinge, mbpx, mbpy);

    /* We then get buckle & opening angles */
    buckle = buckleopening * cos(deg2rad(phi));
    opening = buckleopening * sin(deg2rad(phi));

    /* Finally we get the xyz displacement parameters */
    for (i = 1; i <= 3; i++) {
        tv1[i] = org_lb[i] - org_rb[i];
        shear += tv1[i] * mbpx[i];
        stretch += tv1[i] * mbpy[i];
        stagger += tv1[i] * mbpz[i];
    }

    /* Now the whole set of parameters */
    bp_par[1] = shear;
    bp_par[2] = stretch;
    bp_par[3] = stagger;
    bp_par[4] = buckle;
    bp_par[5] = propeller;
    bp_par[6] = opening;

    free_dvector(hinge, 1, 3);
    free_dvector(tv1, 1, 3);
    free_dvector(tv2, 1, 3);
    free_dvector(tv3, 1, 3);
    free_dvector(mbpx, 1, 3);
    free_dvector(mbpy, 1, 3);
    free_dvector(mbpz, 1, 3);
    free_dmatrix(para_left, 1, 3, 1, 3);
    free_dmatrix(para_right, 1, 3, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
}

void get_step_par(double *step_par, double **mstxyz, double **xyz_lbp,
                  double *org_lbp, double **xyz_ubp, double *org_ubp)
     /* GET_BP_PAR get the 6 local CEHS base-pair parameters */
{
    long i;
    double *hinge, rolltilt, phi;
    double twist, tilt, roll, rise = 0, shift = 0, slide = 0;
    double *tv1, *tv2, *tv3;
    double **para_upper, **para_lower, **rotmat;
    double *mstx, *msty, *mstz;

    hinge = dvector(1, 3);
    tv1 = dvector(1, 3);
    tv2 = dvector(1, 3);
    tv3 = dvector(1, 3);

    for (i = 1; i <= 3; i++) {
        tv1[i] = xyz_lbp[i][3];
        tv2[i] = xyz_ubp[i][3];
    }

    /* Get the hinge by vector product of the z-axes of the two base-pairs */
    cross_dvector(hinge, tv1, tv2);

    /* Get the rolltilt angle.
       It is a positive angle between the z-axes of the two base */
    rolltilt = magang(tv1, tv2);

    para_upper = dmatrix(1, 3, 1, 3);
    para_lower = dmatrix(1, 3, 1, 3);
    rotmat = dmatrix(1, 3, 1, 3);
    /* % We rotate the upper base-pair by -rolltilt/2
       %               lower base-pair by +rolltilt/2
       %   Note: The two base-pairs are now parallel */
    arbrot(rotmat, hinge, -0.5 * rolltilt);
    mul_dmatrix(para_upper, rotmat, xyz_ubp, 3, 3, 3, 3);

    arbrot(rotmat, hinge, 0.5 * rolltilt);
    mul_dmatrix(para_lower, rotmat, xyz_lbp, 3, 3, 3, 3);

    /* The twist angle is that between the two y-axes of the two base-pairs */
    for (i = 1; i <= 3; i++) {
        tv1[i] = para_lower[i][2];
        tv2[i] = para_upper[i][2];
        tv3[i] = para_upper[i][3];
    }
    twist = vec_ang(tv1, tv2, tv3);

    /* We get the MST
       %   Its z axis is the same as the two base_pairs as given above
       %   Its y axis is the mean of the two base_pairs as given above
       %   Its x axis can be got by the vector product of its y & z-axes */
    mstx = dvector(1, 3);
    msty = dvector(1, 3);
    mstz = dvector(1, 3);

    for (i = 1; i <= 3; i++)
        mstz[i] = para_upper[i][3];
    if (fabs(twist - 180) < EPS)  /* Special case: +180 */
        for (i = 1; i <= 3; i++)
            msty[i] = para_upper[i][1];
    else if (fabs(twist + 180) < EPS)  /* Special case: -180 */
        for (i = 1; i <= 3; i++)
            msty[i] = -para_upper[i][1];
    else {
        for (i = 1; i <= 3; i++)
            msty[i] = 0.5 * (para_upper[i][2] + para_lower[i][2]);
        norm_dvector(msty, msty, 3);
    }
    cross_dvector(mstx, msty, mstz);

    for (i = 1; i <= 3; i++) {
        mstxyz[i][1] = mstx[i];
        mstxyz[i][2] = msty[i];
        mstxyz[i][3] = mstz[i];
    }

    /* We then need to get the phi angle which is defined by hinge & msty */
    phi = vec_ang(hinge, msty, mstz);

    /* We then get roll & tilt angles */
    roll = rolltilt * cos(deg2rad(phi));
    tilt = rolltilt * sin(deg2rad(phi));

    /* Finally we get the xyz displacement parameters */
    for (i = 1; i <= 3; i++) {
        tv1[i] = org_ubp[i] - org_lbp[i];
        shift += tv1[i] * mstx[i];
        slide += tv1[i] * msty[i];
        rise += tv1[i] * mstz[i];
    }

    /* Now the whole set of parameters */
    step_par[1] = shift;
    step_par[2] = slide;
    step_par[3] = rise;
    step_par[4] = tilt;
    step_par[5] = roll;
    step_par[6] = twist;

    free_dvector(hinge, 1, 3);
    free_dvector(tv1, 1, 3);
    free_dvector(tv2, 1, 3);
    free_dvector(tv3, 1, 3);
    free_dvector(mstx, 1, 3);
    free_dvector(msty, 1, 3);
    free_dvector(mstz, 1, 3);
    free_dmatrix(para_upper, 1, 3, 1, 3);
    free_dmatrix(para_lower, 1, 3, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
}

void CEHS_par(double **mxyz, double **morg, double **bp_par, double **P_mst,
              double *mtwist, double **bx, double **by, double **bz, double **borg,
              long *P, long numP1, long num_bp, char **bp_seq,
              double **xyz, char *filstr, FILE * fo)
     /* CEHS_PAR get the local CEHS base-pair and step parameters
        %
        % Output: mxyz--(x,y,z)-axes of the MST
        %         morg--origins of the MST
        %         bp_par--local CEHS base-pair parameters
        %         P_mst--phosphorus atom's coordinates with reference to the MST */
{
    long i, j, k, ip1;
    double **xyz_P1, **xyz_P2;

    double **xyz_lb, **xyz_rb, *org_lb, *org_rb;
    double **xyz_ubp, **xyz_lbp, *org_ubp, *org_lbp;
    double **step_par, **mstxyz;  /* step parameters */

    double *dP1, *dP2, *tmp_vec;

    double *ave_par, *std_par;

    char filnam[BUF512];
    FILE *fc;

    ip1 = num_bp - 1;
    xyz_P1 = dmatrix(1, ip1, 1, 3);  /* Terminal P does not consider */
    xyz_P2 = dmatrix(1, ip1, 1, 3);

    if (numP1 == num_bp) {  /* With terminal P */
        for (i = 1; i <= ip1; i++) {
            k = 2 * num_bp - i + 1;
            for (j = 1; j <= 3; j++) {
                xyz_P1[i][j] = xyz[P[i + 1]][j];
                xyz_P2[i][j] = xyz[P[k]][j];
            }
        }
    } else {  /* Normal case */
        for (i = 1; i <= ip1; i++) {
            k = 2 * numP1 - i + 1;  /* numP1 = num_bp-1 */
            for (j = 1; j <= 3; j++) {
                xyz_P1[i][j] = xyz[P[i]][j];
                xyz_P2[i][j] = xyz[P[k]][j];
            }
        }
    }

    /* Get the base-pair parameters */
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

    dP1 = dvector(1, 3);
    dP2 = dvector(1, 3);
    tmp_vec = dvector(1, 3);

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

        for (j = 1; j <= 3; j++) {
            morg[i][j] = 0.5 * (borg[i + 1][j + 6] + borg[i][j + 6]);
            dP1[j] = xyz_P1[i][j] - morg[i][j];
            dP2[j] = xyz_P2[i][j] - morg[i][j];
        }

        dvec_x_dmtx(tmp_vec, dP1, mstxyz, 3, 3, 3);
        for (j = 1; j <= 3; j++)
            P_mst[i][j] = tmp_vec[j];
        dvec_x_dmtx(tmp_vec, dP2, mstxyz, 3, 3, 3);
        /* Reverse y & z coordinates for strand II */
        tmp_vec[2] = -tmp_vec[2];
        tmp_vec[3] = -tmp_vec[3];
        for (j = 1; j <= 3; j++)
            P_mst[i][j + 3] = tmp_vec[j];

        for (j = 1; j <= 3; j++) {
            mxyz[i][j] = mstxyz[j][1];
            mxyz[i][j + 3] = mstxyz[j][2];
            mxyz[i][j + 6] = mstxyz[j][3];
        }
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
        *mtwist = ave_par[6];
        std_dmatrix(std_par, step_par, num_bp - 1, 6);
        fprintf(fo, "    ave.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%10.2f", ave_par[i]);
        fprintf(fo, "\n    s.d.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%10.2f", std_par[i]);
        fprintf(fo, "\n\n");
    } else
        *mtwist = step_par[1][6];  /* Only one step */

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

    free_dvector(org_lb, 1, 3);
    free_dvector(org_rb, 1, 3);
    free_dvector(org_ubp, 1, 3);
    free_dvector(org_lbp, 1, 3);
    free_dvector(dP1, 1, 3);
    free_dvector(dP2, 1, 3);
    free_dvector(tmp_vec, 1, 3);
    free_dvector(ave_par, 1, 6);
    free_dvector(std_par, 1, 6);
    free_dmatrix(xyz_P1, 1, ip1, 1, 3);
    free_dmatrix(xyz_P2, 1, ip1, 1, 3);
    free_dmatrix(xyz_lb, 1, 3, 1, 3);
    free_dmatrix(xyz_rb, 1, 3, 1, 3);
    free_dmatrix(step_par, 1, num_bp - 1, 1, 6);
    free_dmatrix(xyz_ubp, 1, 3, 1, 3);
    free_dmatrix(xyz_lbp, 1, 3, 1, 3);
    free_dmatrix(mstxyz, 1, 3, 1, 3);
}

void wrt_mst(char **atom_id, char *base_type, char *strand, long *base_num,
             double **xyz, char **bp_seq, long num_atom, long num_bp,
             double **mxyz, double **morg, char *filstr)
     /* WRT_MST write out each step in its MST */
{
    long i, j, k, lbn, ns = 0, ne = 0;
    long *bidx, *bsidx, *beidx, *tmp1, *tmp2;
    long *bsidx1, *beidx1, *bsidx2, *beidx2;

    double **mxyzi, **step_mst, **tmp_mtx;
    long step_num, *step_idx, nb1, nb2;

    char filnam[BUF512];
    FILE *fmst;

    lbn = 2 * num_bp;

    /* Get the index for each residue: same code as in <get_data> */
    bidx = lvector(1, num_atom);
    bsidx = lvector(1, lbn);
    beidx = lvector(1, lbn);
    tmp1 = lvector(1, num_atom - 1);
    tmp2 = lvector(1, num_atom + 1);

    for (i = 1; i <= num_atom; i++)
        bidx[i] = base_num[i] + 100 * base_type[i];  /* char as integer & use bidx temp */

    i_diff(tmp1, bidx, num_atom);
    logical_not(tmp1, tmp1, num_atom - 1);

    tmp2[1] = 0;
    tmp2[num_atom + 1] = 0;
    for (i = 1; i <= num_atom - 1; i++)
        tmp2[i + 1] = tmp1[i];

    i_diff(bidx, tmp2, num_atom + 1);

    for (i = 1; i <= num_atom; i++) {
        if (bidx[i] == 1) {
            ns++;
            if (ns > lbn)
                nrerror("Something wrong with naming in you PDB file!");
            bsidx[ns] = i;
        }
        if (bidx[i] == -1) {
            ne++;
            if (ne > lbn)
                nrerror("Something wrong with naming in you PDB file!");
            beidx[ne] = i;
        }
    }  /* ns=ne=lbn: 2*num_bp */

    if (ns != lbn)
        nrerror("Something wrong with naming in you PDB file!");

    bsidx1 = lvector(1, num_bp);
    beidx1 = lvector(1, num_bp);
    bsidx2 = lvector(1, num_bp);
    beidx2 = lvector(1, num_bp);

    for (i = 1; i <= num_bp; i++) {
        bsidx1[i] = bsidx[i];
        beidx1[i] = beidx[i];
        bsidx2[i] = bsidx[num_bp + i];
        beidx2[i] = beidx[num_bp + i];
    }

    strcpy(filnam, filstr);
    fmst = open_file(strcat(filnam, "mst"), "w");

    mxyzi = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_bp - 1; i++) {
        for (j = 1; j <= 3; j++)
            mxyzi[j][1] = mxyz[i][j];
        for (j = 1; j <= 3; j++)
            mxyzi[j][2] = mxyz[i][j + 3];
        for (j = 1; j <= 3; j++)
            mxyzi[j][3] = mxyz[i][j + 6];

        nb1 = beidx1[i + 1] - bsidx1[i] + 1;
        nb2 = beidx2[num_bp - i + 1] - bsidx2[num_bp - i] + 1;
        step_num = nb1 + nb2;
        step_idx = lvector(1, step_num);

        for (j = 1; j <= nb1; j++)
            step_idx[j] = bsidx1[i] + j - 1;
        for (j = 1; j <= nb2; j++)
            step_idx[j + nb1] = bsidx2[num_bp - i] + j - 1;

        tmp_mtx = dmatrix(1, step_num, 1, 3);
        step_mst = dmatrix(1, step_num, 1, 3);

        for (j = 1; j <= step_num; j++)
            for (k = 1; k <= 3; k++)
                tmp_mtx[j][k] = xyz[step_idx[j]][k] - morg[i][k];
        mul_dmatrix(step_mst, tmp_mtx, mxyzi, step_num, 3, 3, 3);

        fprintf(fmst, "HEADER    Base-Pair step #%2.2ld with reference to the MST\n", i);
        fprintf(fmst, "COMPND    Strand  I in 5'-->3': %c%c\n", bp_seq[i][1],
                bp_seq[i + 1][1]);
        fprintf(fmst, "COMPND    Strand II in 3'-->5': %c%c\n", bp_seq[i][2],
                bp_seq[i + 1][2]);
        fprintf(fmst, "AUTHOR    Program written by Xiang-Jun Lu (1996--1997)\n");

        for (j = 1; j <= step_num; j++) {
            k = step_idx[j];
            fprintf(fmst, "ATOM  %5ld  %3s %3c %c%4ld    %8.3f%8.3f%8.3f\n",
                    j, atom_id[k], base_type[k], strand[k], base_num[k],
                    step_mst[j][1], step_mst[j][2], step_mst[j][3]);
        }
        fprintf(fmst, "END\n");

        free_lvector(step_idx, 1, step_num);
        free_dmatrix(tmp_mtx, 1, step_num, 1, 3);
        free_dmatrix(step_mst, 1, step_num, 1, 3);
    }
    fclose(fmst);

    free_lvector(bidx, 1, num_atom);
    free_lvector(bsidx, 1, lbn);
    free_lvector(beidx, 1, lbn);
    free_lvector(tmp1, 1, num_atom - 1);
    free_lvector(tmp2, 1, num_atom + 1);
    free_lvector(bsidx1, 1, num_bp);
    free_lvector(beidx1, 1, num_bp);
    free_lvector(bsidx2, 1, num_bp);
    free_lvector(beidx2, 1, num_bp);
    free_dmatrix(mxyzi, 1, 3, 1, 3);
}

void GLH_par(double **bx, double **by, double **bz, double **borg,
             long num_bp, char **bp_seq, double **bp_par, char *filstr, FILE * fo)
     /*
        % GLH_PAR calculate the global helical parameters which are consistent
        %         with the Cambridge Convention and are completely reversible */
{
    long i, j;
    double *gy, *gz, *hinge, *bpz, *TipInc;
    double **bpxyz_i, **rotmat, **tmpdat;
    double *z_rot, *z_disp, *tip, *inc, *x_disp, *y_disp;
    double *twist, *rise, *bpx_p, *bpy_p, *bporg, theta_angle;
    double mrise, mtwist, nbp_turn, repeat_dist;

    char filnam[BUF512], bstr[BUF512];
    FILE *fg;

    gy = dvector(1, 3);
    gz = dvector(1, 3);
    hinge = dvector(1, 3);
    bpz = dvector(1, 3);
    bpx_p = dvector(1, 3);
    bpy_p = dvector(1, 3);
    bporg = dvector(1, 3);
    TipInc = dvector(1, num_bp);

    z_rot = dvector(1, num_bp);
    z_disp = dvector(1, num_bp);
    tip = dvector(1, num_bp);
    inc = dvector(1, num_bp);
    x_disp = dvector(1, num_bp);
    y_disp = dvector(1, num_bp);
    twist = dvector(1, num_bp);
    rise = dvector(1, num_bp);

    bpxyz_i = dmatrix(1, 3, 1, 3);
    rotmat = dmatrix(1, 3, 1, 3);
    tmpdat = dmatrix(1, 3, 1, 3);

    gy[1] = 0.0;
    gy[2] = 1.0;
    gy[3] = 0.0;
    gz[1] = 0.0;
    gz[2] = 0.0;
    gz[3] = 1.0;

    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++)
            bpz[j] = bz[i][j + 6];
        cross_dvector(hinge, gz, bpz);
        TipInc[i] = magang(gz, bpz);

        for (j = 1; j <= 3; j++)
            bpxyz_i[j][1] = bx[i][j + 6];
        for (j = 1; j <= 3; j++)
            bpxyz_i[j][2] = by[i][j + 6];
        for (j = 1; j <= 3; j++)
            bpxyz_i[j][3] = bz[i][j + 6];

        arbrot(rotmat, hinge, -TipInc[i]);
        mul_dmatrix(tmpdat, rotmat, bpxyz_i, 3, 3, 3, 3);
        for (j = 1; j <= 3; j++)
            bpx_p[j] = tmpdat[j][1];
        for (j = 1; j <= 3; j++)
            bpy_p[j] = tmpdat[j][2];

        for (j = 1; j <= 3; j++)
            bporg[j] = borg[i][j + 6];

        /* Get z_rot & z_disp with reference to the global reference triad
           Note: by defintion, z_rot[1]=0.0 & z_disp[1]=0.0 */
        z_rot[i] = vec_ang(gy, bpy_p, gz);
        z_disp[i] = bporg[3];  /* dot_dvector(bporg, gz, 3); */

        /* TIP and INCLINATION */
        theta_angle = deg2rad(vec_ang(hinge, bpy_p, gz));
        tip[i] = TipInc[i] * cos(theta_angle);
        inc[i] = TipInc[i] * sin(theta_angle);

        /* X & Y-displacement: the projections of the bporg's x/y components
           %                   along x and y-axes */
        x_disp[i] = dot_dvector(bporg, bpx_p, 3);
        y_disp[i] = dot_dvector(bporg, bpy_p, 3);
    }

    /* TWIST and RISE: Twist is defined as the angle between successive y-axes
       %               Rise is the difference between bporg's projection along
       %                 the helical axis.
       % See also the MATLAB version */
    for (i = 1; i <= num_bp - 1; i++) {
        twist[i] = z_rot[i + 1] - z_rot[i];
        if (twist[i] < -180.0)
            twist[i] += 360.0;
        if (twist[i] > +180.0)
            twist[i] += -360.0;
        rise[i] = z_disp[i + 1] - z_disp[i];
    }

    /* Mean rise & twist */
    mrise = mean_dvector(rise, num_bp - 1);
    mtwist = mean_dvector(twist, num_bp - 1);
    if ((mtwist < 0) && ((num_bp - 1) % 2)) {  /* Z- and W-DNA: 2-bp a unit */
        mrise = mean_dvector(rise, num_bp - 2);
        mtwist = mean_dvector(twist, num_bp - 2);
    }

    /* Write out the result */
    prt_sep(fo, '*', 76);
    fprintf(fo, "Global helical parameters by Lu & Hunter\n");
    fprintf(fo, "    Note: Would be MEANINGLESS for strongly curved structures!\n\n");
    fprintf(fo, "    bp       X-disp    Y-disp     Rise");
    fprintf(fo, "     Incl.      Tip      Twist\n");
    for (i = 1; i <= num_bp - 1; i++) {
        fprintf(fo, " %2ld %c-%c ", i, bp_seq[i][1], bp_seq[i][2]);
        fprintf(fo, "%10.2f%10.2f%10.2f", x_disp[i], y_disp[i], rise[i]);
        fprintf(fo, "%10.2f%10.2f%10.2f\n", inc[i], tip[i], twist[i]);
    }

    i = num_bp;
    strcpy(bstr, "       -- ");
    fprintf(fo, " %2ld %c-%c ", i, bp_seq[i][1], bp_seq[i][2]);
    fprintf(fo, "%10.2f%10.2f%s", x_disp[i], y_disp[i], bstr);
    fprintf(fo, "%10.2f%10.2f%s\n", inc[i], tip[i], bstr);
    fprintf(fo, "          ");
    prt_sep(fo, '~', 58);

    if (num_bp > 2) {
        fprintf(fo, "    ave.");
        fprintf(fo, "%10.2f", mean_dvector(x_disp, num_bp));
        fprintf(fo, "%10.2f", mean_dvector(y_disp, num_bp));
        fprintf(fo, "%10.2f", mean_dvector(rise, num_bp - 1));
        fprintf(fo, "%10.2f", mean_dvector(inc, num_bp));
        fprintf(fo, "%10.2f", mean_dvector(tip, num_bp));
        fprintf(fo, "%10.2f\n", mean_dvector(twist, num_bp - 1));

        fprintf(fo, "    s.d.");
        fprintf(fo, "%10.2f", std_dvector(x_disp, num_bp));
        fprintf(fo, "%10.2f", std_dvector(y_disp, num_bp));
        fprintf(fo, "%10.2f", std_dvector(rise, num_bp - 1));
        fprintf(fo, "%10.2f", std_dvector(inc, num_bp));
        fprintf(fo, "%10.2f", std_dvector(tip, num_bp));
        fprintf(fo, "%10.2f\n", std_dvector(twist, num_bp - 1));
    }

    /* Print out overall helical information */
    nbp_turn = 360 / fabs(mtwist);
    repeat_dist = nbp_turn * mrise;
    prt_sep(fo, '-', 76);
    fprintf(fo, "Overall helical parameters:\n\n");
    fprintf(fo, "Mean helical rise per residue %6.2f Angstrom\n", mrise);
    fprintf(fo, "Mean helical twist per residue %6.2f degrees\n", mtwist);
    fprintf(fo, "There are %5.2f base-pairs per turn\n", nbp_turn);
    fprintf(fo, "The repeat distance is %6.2f Angstrom\n\n", repeat_dist);

    /* Write bp & global parameter data file for rebuilding process! */
    strcpy(filnam, filstr);
    fg = open_file(strcat(filnam, "glh"), "w");
    fprintf(fg, "%4ld base-pairs\n", num_bp);
    fprintf(fg, "       Shear     Str    Stag    Buck    Prop    Open");
    fprintf(fg, "   X-dsp   Y-dsp    Rise    Incl    Tip    Twist\n");

    i = 1;
    fprintf(fg, "%c-%c ", bp_seq[i][1], bp_seq[i][2]);
    for (j = 1; j <= 6; j++)
        fprintf(fg, "%8.2f", bp_par[i][j]);
    fprintf(fg, "%8.2f%8.2f%8.2f", x_disp[i], y_disp[i], 0.0);
    fprintf(fg, "%8.2f%8.2f%8.2f\n", inc[i], tip[i], 0.0);

    for (i = 2; i <= num_bp; i++) {
        fprintf(fg, "%c-%c ", bp_seq[i][1], bp_seq[i][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fg, "%8.2f", bp_par[i][j]);
        fprintf(fg, "%8.2f%8.2f%8.2f", x_disp[i], y_disp[i], rise[i - 1]);
        fprintf(fg, "%8.2f%8.2f%8.2f\n", inc[i], tip[i], twist[i - 1]);
    }
    fclose(fg);

    free_dvector(gy, 1, 3);
    free_dvector(gz, 1, 3);
    free_dvector(hinge, 1, 3);
    free_dvector(bpz, 1, 3);
    free_dvector(bpx_p, 1, 3);
    free_dvector(bpy_p, 1, 3);
    free_dvector(bporg, 1, 3);
    free_dvector(TipInc, 1, num_bp);
    free_dvector(z_rot, 1, num_bp);
    free_dvector(z_disp, 1, num_bp);
    free_dvector(tip, 1, num_bp);
    free_dvector(inc, 1, num_bp);
    free_dvector(x_disp, 1, num_bp);
    free_dvector(y_disp, 1, num_bp);
    free_dvector(twist, 1, num_bp);
    free_dvector(rise, 1, num_bp);
    free_dmatrix(bpxyz_i, 1, 3, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
    free_dmatrix(tmpdat, 1, 3, 1, 3);
}

void classify(double **P_mst, long num_bp, char **bp_seq, long xdir,
              double mtwist, FILE * fo)
     /* % CLASSIFY write out phosphorus atom's xyz coordinates with regard
        %          to the MST, and classify the structure into A, B, Z, W
        %          or R-form */
{
    long i, j;
    double *ave_P, *std_P, Zp;

    /* xyz coordinates of phosphorus atoms with regard to the MST */
    prt_sep(fo, '*', 76);
    fprintf(fo, "The xyz coordinates of phosphorus atoms with respect to");
    fprintf(fo, " the middle-step-triad");
    if (mtwist > 0)
        fprintf(fo,
                "\n\n  step      xI      yI      zI      xII     yII     zII     Form\n");
    else
        fprintf(fo, "\n\n  step      xI      yI      zI      xII     yII     zII\n");

    for (i = 1; i <= num_bp - 1; i++) {
        fprintf(fo, "%2ld %c%c/%c%c", i, bp_seq[i][1], bp_seq[i + 1][1],
                bp_seq[i + 1][2], bp_seq[i][2]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%8.2f", P_mst[i][j]);
        if (mtwist > 0) {
            Zp = (P_mst[i][3] + P_mst[i][6]) / 2;  /* Average for the step */
            if (Zp > 1.5)
                fprintf(fo, "  |  A\n");
            else if (Zp < 0.5)
                fprintf(fo, "  |  B\n");
            else
                fprintf(fo, "  |  N\n");
        } else
            fprintf(fo, "\n");
    }
    ave_P = dvector(1, 6);
    std_P = dvector(1, 6);

    if (num_bp > 2) {
        mean_dmatrix(ave_P, P_mst, num_bp - 1, 6);
        std_dmatrix(std_P, P_mst, num_bp - 1, 6);
        fprintf(fo, "        ");
        prt_sep(fo, '~', 56);
        fprintf(fo, "    ave.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%8.2f", ave_P[i]);
        fprintf(fo, "\n    s.d.");
        for (i = 1; i <= 6; i++)
            fprintf(fo, "%8.2f", std_P[i]);
        fprintf(fo, "\n\n");
        Zp = (ave_P[3] + ave_P[6]) / 2.0;
    } else  /* num_bp=2, only one middle P on each strand */
        Zp = (P_mst[1][3] + P_mst[1][6]) / 2.0;

    prt_sep(fo, '*', 76);
    fprintf(fo, "Structure classification: \n\n");

    if (labs(xdir) != num_bp - 1)
        fprintf(fo, "This structure has bases flipped over.\n");

    if (xdir > 0)
        fprintf(fo, "Positive x-axis points from MINOR groove to MAJOR groove\n\n");
    else
        fprintf(fo, "Positive x-axis points from MAJOR groove to MINOR groove\n\n");

    if (mtwist < 0) {  /* Left handed */
        if (xdir > 0)
            fprintf(fo, "This is a left-handed W-form DNA\n");
        else
            fprintf(fo, "This is a left-handed Z-form DNA\n");
    } else {  /* Right handed */
        fprintf(fo, "This is a right-handed double helix\n");
        if (xdir > 0) {  /* Either A or B */
            if (Zp > 1.5)
                fprintf(fo, "  Since its mean Zp (%6.2f) > 1.5A, it is in A-form\n", Zp);
            else if (Zp < 0.5)
                fprintf(fo, "  Since its mean Zp (%6.2f) < 0.5A, it is in B-form\n", Zp);
            else
                fprintf(fo,
                        "  Mean Zp (%6.2f) in [0.5--1.5]A, it is between A- and B-form\n",
                        Zp);
        } else  /* Non-exist form! But logically possible */
            fprintf(fo, "This is an R-form DNA (non-experimental evidence yet)!\n");
    }

    fprintf(fo, "\n");

    free_dvector(ave_P, 1, 6);
    free_dvector(std_P, 1, 6);
}
