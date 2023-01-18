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

/* %%%%%%%%%%%%%%%%%%% Common functions BEGIN %%%%%%%%%%%%%%%%%%% */

/* ALIGN_HELIX align the helix using an axis defined by <org_xyz> */
void align_helix(double **oxyz, double **xyz, long natom, double **org_xyz, long nbp)
{
    char str[BUF512];
    long ich, i, j;
    double *hel_axis, *ave_xyz, *pnts_dist, **diff_org_xyz, dvar;
    double **rotmat, **rotmat_T, **tmpxyz;

/* ----- Commented to probe for efficiency ----------- */

ich = 1;

/*

    printf("\n1. View perpendicular to the 1st base-pair\n");
    printf("2. View parallel to the 1st base-pair\n");
    printf("3. View perpendicular to the 'best-fit' line through mid-bps\n");
    printf("4. View perpendicular to the 'best-fit' helical axis\n");
    printf("Your choice (1-4, Dft 1): ");
    read_stdin_str(str);
    ich = atoi(str);

*/

    if ((ich < 1) || (ich > 4))
        ich = 1;

    hel_axis = dvector(1, 3);
    ave_xyz = dvector(1, 3);
    pnts_dist = dvector(1, nbp);
    diff_org_xyz = dmatrix(1, nbp - 1, 1, 3);

    if (ich == 1) {
        for (i = 1; i <= natom; i++)
            for (j = 1; j <= 3; j++)
                oxyz[i][j] = xyz[i][j];
        free_dvector(hel_axis, 1, 3);
        free_dvector(ave_xyz, 1, 3);
        free_dvector(pnts_dist, 1, nbp);
        free_dmatrix(diff_org_xyz, 1, nbp - 1, 1, 3);
        return;
    } else if (ich == 2) {
        hel_axis[1] = 0.0;
        hel_axis[2] = 0.0;
        hel_axis[3] = 1.0;
    } else if (ich == 3)
        lsline(hel_axis, ave_xyz, pnts_dist, org_xyz, nbp);
    else {
        for (i = 1; i <= nbp - 1; i++)
            for (j = 1; j <= 3; j++)
                diff_org_xyz[i][j] = org_xyz[i + 1][j] - org_xyz[i][j];
        lsplane(hel_axis, &dvar, pnts_dist, diff_org_xyz, nbp - 1);
    }

    rotmat = dmatrix(1, 3, 1, 3);
    rotmat_T = dmatrix(1, 3, 1, 3);
    tmpxyz = dmatrix(1, natom, 1, 3);

    alignz(rotmat, hel_axis);  /* rotation matrix */
    tr_dmatrix(rotmat_T, rotmat, 3, 3);  /* its transpose */
    mul_dmatrix(oxyz, xyz, rotmat_T, natom, 3, 3, 3);

    rotz(rotmat, -90.0);
    mul_dmatrix(tmpxyz, oxyz, rotmat, natom, 3, 3, 3);

    rotx(rotmat, 90.0);
    mul_dmatrix(oxyz, tmpxyz, rotmat, natom, 3, 3, 3);

    /* matlab code: xyz=xyz*alignz(hel_axis)'*rotz(90)'*rotx(-90)'; */

    free_dvector(hel_axis, 1, 3);
    free_dvector(ave_xyz, 1, 3);
    free_dvector(pnts_dist, 1, nbp);
    free_dmatrix(diff_org_xyz, 1, nbp - 1, 1, 3);
    free_dmatrix(rotmat, 1, 3, 1, 3);
    free_dmatrix(rotmat_T, 1, 3, 1, 3);
    free_dmatrix(tmpxyz, 1, natom, 1, 3);
}

/* xdir--0 for minor-->major as in A/B/W; 1 otherwise as in Z-DNA */
long get_xdir(double twist)
{
    char yn, str[BUF512];
    long xdir;

/* ----- Commented to probe for efficiency ----------- */

xdir = 0;

/*
    if (twist < 0.0) {  // Left-handed DNA
        printf("\nPositive x-axis from MAJOR-->minor groove (y/n, Dft: y): ");
        read_stdin_str(str);
        yn = toupper(str[0]);
        if (yn != 'N')
            xdir = 1;
        else
            xdir = 0;
    } else {  // Right-handed DNA
        printf("\nPositive x-axis from minor-->MAJOR groove (y/n, Dft: y): ");
        read_stdin_str(str);
        yn = toupper(str[0]);
        if (yn != 'N')
            xdir = 0;
        else
            xdir = 1;
    }

*/
    return xdir;
}

/* STR_DISPLAY graphical display of the generated structure using RasMol
   %
   % Input: imdl--representation model, 1 for atomic, and 2 for schematic */
void str_display(long imdl, char *filnam)
{
    char *str, *cmd_str, *env_add, yn;

    str = cvector(1, BUF512);
    cmd_str = cvector(1, BUF512);

    printf("\nDisplay of the generated structure [Y(Dft--RasMol)/N]: \n");
    read_stdin_str(str);
    yn = toupper(str[0]);

    if (yn == 'N') {
        free_cvector(str, 1, BUF512);
        free_cvector(cmd_str, 1, BUF512);
        return;
    }

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(str, env_add);
    strcat(str, "/Utilities/");

    if (imdl == 1) {  /* Atomic model: PDB format */
        strcpy(cmd_str, "rasmol ");
        strcat(cmd_str, filnam);
        strcat(cmd_str, " -script ");
        strcat(cmd_str, str);
        strcat(cmd_str, "ATMrasmol.scr");
        system(cmd_str);
    } else {  /* Rectangular block model: ALCHEMY format */
        strcpy(cmd_str, "rasmol -alchemy ");
        strcat(cmd_str, filnam);
        strcat(cmd_str, " -script ");
        strcat(cmd_str, str);
        strcat(cmd_str, "BLKrasmol.scr");
        system(cmd_str);
    }
    free_cvector(str, 1, BUF512);
    free_cvector(cmd_str, 1, BUF512);
}

void set_bppar(double **bppar, long nbp, char **bpseq)
     /* Set proper base-pair parameters for A-T & G-C */
{
    long i, j, ibp;
    char *str, *env_add, bpfile[BUF512], bchar;
    char A_T[4], G_C[4];
    double *A_Tpar, *G_Cpar, *T_Apar, *C_Gpar;
    FILE *fbp;

    str = cvector(1, BUF512);


/* ----- Commented to probe for efficiency ----------- */

ibp = 1;

/*

    printf("\nWhich base-pair geometry to use:\n");
    printf("1. Flat base-pair\n");
    printf("2. Non-flat base-pair (overall X-ray average)\n");
    printf("3. X-ray average for A-T and G-C separately\n");
    printf("4. User defined geometry\n");
    printf("Your choice (1-4, Dft 1): ");
    read_stdin_str(str);
    ibp = atoi(str);
*/

    if ((ibp < 1) || (ibp > 4))
        ibp = 1;

    if (ibp == 1)
        strcpy(bpfile, "BPflat.dat");
    else if (ibp == 2)
        strcpy(bpfile, "BPnflat.dat");
    else if (ibp == 3)
        strcpy(bpfile, "BPxray.dat");
    else {
        printf("Name of the base-pair geometry file: ");
        read_stdin_str(bpfile);
    }

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(str, env_add);
    strcat(str, "/BaseGeo/");
    strcat(str, bpfile);

    fbp = open_file(str, "r");

    A_Tpar = dvector(1, 6);
    G_Cpar = dvector(1, 6);
    T_Apar = dvector(1, 6);
    C_Gpar = dvector(1, 6);

    fgets(str, BUF512 - 1, fbp);  /* Skip the first title line */

    /* For A-T & T-A base-pairs */
    fscanf(fbp, "%s", A_T);
    upper(A_T, strlen(A_T));
    if (strcmp(A_T, "A-T"))
        nrerror("The first bp should be A-T!");
    for (i = 1; i <= 6; i++) {
        fscanf(fbp, "%lf", &A_Tpar[i]);
        if (i == 1 || i == 4)
            T_Apar[i] = -A_Tpar[i];
        else
            T_Apar[i] = A_Tpar[i];
    }

    /* For G-C & C-G base-pairs */
    fscanf(fbp, "%s", G_C);
    upper(G_C, strlen(G_C));
    if (strcmp(G_C, "G-C"))
        nrerror("The second bp should be G-C!");
    for (i = 1; i <= 6; i++) {
        fscanf(fbp, "%lf", &G_Cpar[i]);
        if (i == 1 || i == 4)
            C_Gpar[i] = -G_Cpar[i];
        else
            C_Gpar[i] = G_Cpar[i];
    }
    fclose(fbp);

    /* Set base-pair parameters for each pair in the sequence */
    for (i = 1; i <= nbp; i++) {
        bchar = bpseq[1][i - 1];
        if (bchar == 'A')
            for (j = 1; j <= 6; j++)
                bppar[i][j] = A_Tpar[j];
        else if (bchar == 'T')
            for (j = 1; j <= 6; j++)
                bppar[i][j] = T_Apar[j];
        else if (bchar == 'G')
            for (j = 1; j <= 6; j++)
                bppar[i][j] = G_Cpar[j];
        else
            for (j = 1; j <= 6; j++)
                bppar[i][j] = C_Gpar[j];
    }

    free_cvector(str, 1, BUF512);
    free_dvector(A_Tpar, 1, 6);
    free_dvector(G_Cpar, 1, 6);
    free_dvector(T_Apar, 1, 6);
    free_dvector(C_Gpar, 1, 6);
}

/* Connect successive blocks of ALCHEMY data file */
/* This is a function called directly by SCHNArP, but
   there is also a stand alone command version */
void alc_connect(char *filinp, char *filout)
{
    long i, num, nbond, **lkg;
    char **asym;
    double **xyz;
    long n_O1, n_O2, n_N1, n_N2;
    long *O1, *O2, *N1, *N2;

    lkg = lmatrix(1, MAX_ATOM, 1, 2);
    asym = cmatrix(1, MAX_ATOM, 1, 3);  /* 2-characters long */
    xyz = dmatrix(1, MAX_ATOM, 1, 3);

    rdalc(&num, &nbond, asym, lkg, xyz, filinp);

    O1 = lvector(1, num);
    O2 = lvector(1, num);
    N1 = lvector(1, num);
    N2 = lvector(1, num);

    /* Note: No "O1" atom type in ALCHEMY */
    str_idx(&n_O1, O1, asym, "O3", num);
    str_idx(&n_O2, O2, asym, "O2", num);
    str_idx(&n_N1, N1, asym, "N1", num);
    str_idx(&n_N2, N2, asym, "N2", num);

    if (n_O1 == 0) {  /* No C1' atoms, so connect successive N atoms */
        /* Strand I */
        for (i = 1; i <= (n_N1 - 1); i++) {
            nbond++;
            lkg[nbond][1] = N1[i];
            lkg[nbond][2] = N1[i + 1];
        }
        /* Strand II */
        for (i = 1; i <= (n_N2 - 1); i++) {
            nbond++;
            lkg[nbond][1] = N2[i];
            lkg[nbond][2] = N2[i + 1];
        }
    } else {  /* There are O atoms */
        /* Strand I */
        for (i = 1; i <= (n_O1 - 1); i++) {
            nbond++;
            lkg[nbond][1] = O1[i];
            lkg[nbond][2] = O1[i + 1];
        }
        /* Strand II */
        for (i = 1; i <= (n_O2 - 1); i++) {
            nbond++;
            lkg[nbond][1] = O2[i];
            lkg[nbond][2] = O2[i + 1];
        }
    }

    for (i = 1; i <= num; i++)
        if (!strcmp(asym[i], "BR"))
            strcpy(asym[i], "O2");  /* O3 will also do */

    wrtalc(num, nbond, asym, lkg, xyz, filout);

    free_lvector(O1, 1, num);
    free_lvector(O2, 1, num);
    free_lvector(N1, 1, num);
    free_lvector(N2, 1, num);
    free_lmatrix(lkg, 1, MAX_ATOM, 1, 2);
    free_cmatrix(asym, 1, MAX_ATOM, 1, 3);
    free_dmatrix(xyz, 1, MAX_ATOM, 1, 3);
}

/* %%%%%%%%%%%%%%%%%%% Common functions END %%%%%%%%%%%%%%%%%%% */

/* %%%%%%%%%%%%%%%%%%% CEHS Building functions BEGIN %%%%%%%%%%%%%%%%%%% */

/* WRT_SEQ_CEHS write the sequence and its associated CEHS step parameters
   %             This file can be used for error check or fined-turned to
   %             get the desired result.
   %
   % Input: nbp--number of base-pairs
   %        bpseq--[2 by nbp] base-pair sequence
   %        steppar--[nbp by 6] CEHS step parameters */
void wrt_seq_CEHS(long nbp, char **bpseq, double **steppar)
{
    long i, j;
    char fstep[BUF512];
    FILE *fo;

    printf("\nOutput sequence CEHS parameter data file name (Dft CEHS_seq.dat): ");
    read_stdin_str(fstep);
    if (strlen(fstep) == 0)
        strcpy(fstep, "CEHS_seq.dat");

    fo = open_file(fstep, "w");

    fprintf(fo, "%5ld base-pairs\n", nbp);
    fprintf(fo, "       Shift   Slide    Rise    Tilt    Roll   Twist\n");
    for (i = 1; i <= nbp; i++) {
        fprintf(fo, "%c-%c ", bpseq[1][i - 1], bpseq[2][i - 1]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%8.2f", steppar[i][j]);
        fprintf(fo, "\n");
    }
    fclose(fo);
}

/*  STD_CEHS build A, B, C, Z and other regular double helical
    structure using CEHS step parameters */
void std_CEHS(long *nbp, char **bpseq, double **steppar)
{
    long i, j, ich;
    char str[BUF512];
    double *stepblk, *Z_Cg, *Z_Gc;
    double d1, d2, d3, d4, d5, d6;

    printf("1. Standard A-form DNA (Twist=32, Roll=12, Slide=-1.5A)\n");
    printf("2. Standard B-form DNA (Twist=36, Roll= 0, Slide=   0A)\n");
    printf("3. Standard C-form DNA (Twist=40, Roll=-6, Slide=   1A)\n");
    printf("4. Standard Z-form DNA (Twist=-10, Slide=5A for CG step)\n");
    printf("                        Twist=-50, Slide=0A for GC step)\n");
    printf("5. Other regular double helical structures\n");
    printf("Input your choice (1-5, Dft 2): ");
    read_stdin_str(str);
    ich = atoi(str);
    if ((ich < 1) || (ich > 5))
        ich = 2;

    stepblk = dvector(1, 6);
    Z_Cg = dvector(1, 6);
    Z_Gc = dvector(1, 6);

    /* Get the step parameters for building blocks */
    switch (ich) {
    case 1:  /* A-form: [0 -1.5 3.23 0 12 32] */
        stepblk[1] = 0.0;
        stepblk[2] = -1.5;
        stepblk[3] = 3.23;
        stepblk[4] = 0.0;
        stepblk[5] = 12.0;
        stepblk[6] = 32.0;
        break;
    case 2:  /* B-form: [0 0 3.38 0 0 36] */
        stepblk[1] = 0.0;
        stepblk[2] = 0.0;
        stepblk[3] = 3.38;
        stepblk[4] = 0.0;
        stepblk[5] = 0.0;
        stepblk[6] = 36.0;
        break;
    case 3:  /* C-form: [0 1 3.36 0 -6 40] */
        stepblk[1] = 0.0;
        stepblk[2] = 1.0;
        stepblk[3] = 3.36;
        stepblk[4] = 0.0;
        stepblk[5] = -6.0;
        stepblk[6] = 40.0;
        break;
    case 4:  /* Z-form: Two kinds of steps */
        /* CG/CG [0 5 3.19 0 0 -10] */
        Z_Cg[1] = 0.0;
        Z_Cg[2] = 5.0;
        Z_Cg[3] = 3.19;
        Z_Cg[4] = 0.0;
        Z_Cg[5] = 0.0;
        Z_Cg[6] = -10.0;
        /* GC/GC [0 0 3.85 0 0 -50] */
        Z_Gc[1] = 0.0;
        Z_Gc[2] = 0.0;
        Z_Gc[3] = 3.85;
        Z_Gc[4] = 0.0;
        Z_Gc[5] = 0.0;
        Z_Gc[6] = -50.0;
        break;
    default:  /* Other regular form */
        printf("\nInput 6 step parameters in the order of: \n");
        printf("Shift  Slide  Rise  Tilt  Roll  Twist\n");
        while (1) {  /* Infinite loop */
            read_stdin_str(str);
            j = sscanf(str, "%lf%lf%lf%lf%lf%lf", &d1, &d2, &d3, &d4, &d5, &d6);
            if (j == 6)
                break;
            printf("Six step-parameters needed. Try again!\n");
        }
        stepblk[1] = d1;
        stepblk[2] = d2;
        stepblk[3] = d3;
        stepblk[4] = d4;
        stepblk[5] = d5;
        stepblk[6] = d6;
    }

    /* Read in base sequence */
    get_sequence(bpseq, nbp);

    for (i = 1; i <= 6; i++)
        steppar[1][i] = 0.0;

    /* Get the whole set of step parameters */
    if (ich != 4)  /* Not Z-DNA */
        for (i = 2; i <= (*nbp); i++)
            for (j = 1; j <= 6; j++)
                steppar[i][j] = stepblk[j];
    else {
        for (i = 2; i <= (*nbp); i += 2)  /* CG/CG steps */
            for (j = 1; j <= 6; j++)
                steppar[i][j] = Z_Cg[j];
        for (i = 3; i <= (*nbp); i += 2)  /* CG/CG steps */
            for (j = 1; j <= 6; j++)
                steppar[i][j] = Z_Gc[j];
    }

    /* Write out base-pair sequence & step parameters */
    wrt_seq_CEHS(*nbp, bpseq, steppar);

    free_dvector(stepblk, 1, 6);
    free_dvector(Z_Cg, 1, 6);
    free_dvector(Z_Gc, 1, 6);
}

/* % DIMER_CEHS build arbitrary sequence DNA structure using various
   %           dimer models
   %
   % Note: All existing models can ONLY handle bend/curvature as a function
   %       of "Twist, Roll and/or Tilt" of the ten dinucleotide steps.
   %
   %       The program given here is universal in that it accepts all six
   %       step parameters as in the X-ray model */
void dimer_CEHS(long *nbp, char **bpseq, double **steppar)
{
    long i, j, k, ich;
    char *str, *env_add, mfilnam[BUF512];
    double **dmstep, *tmppar;
    FILE *fm;

    char bpstep[3];

    /* --------------------------------------------------------------------
       This initialization methods works for ANSI-C. But for unknown
       reasons, it does NOT work on HP computers ......
       See also <cehs.c> & <analysis.c>

       Order of the 16 dinucleotide steps

       char *dncl[17] = {"XX", "AA", "GA", "AG", "GG", "AT", "AC", "GC", "TA",
       "CA", "CG", "TT", "TC", "CT", "CC", "GT", "TG"};

       -------------------------------------------------------------------- */

    /* The following stupid method should always work */
    char *dncl[17];
    dncl[0] = "XX";
    dncl[1] = "AA";
    dncl[2] = "GA";
    dncl[3] = "AG";
    dncl[4] = "GG";
    dncl[5] = "AT";
    dncl[6] = "AC";
    dncl[7] = "GC";
    dncl[8] = "TA";
    dncl[9] = "CA";
    dncl[10] = "CG";
    dncl[11] = "TT";
    dncl[12] = "TC";
    dncl[13] = "CT";
    dncl[14] = "CC";
    dncl[15] = "GT";
    dncl[16] = "TG";
    /* Stop here ====================================== */

    str = cvector(1, BUF512);

    printf("Choose from one of the following dimer models\n");
    printf("1. Calladine\n");
    printf("2. Bolshoy\n");
    printf("3. De Santis\n");
    printf("4. Average step parameters from X-ray structures\n");
    printf("5. User defined model\n");
    printf("Your choice (1-5, Dft 4): ");
    read_stdin_str(str);
    ich = atoi(str);

    if ((ich < 1) || (ich > 5))
        ich = 4;

    if (ich == 1)
        strcpy(mfilnam, "dimer_CDM");
    else if (ich == 2)
        strcpy(mfilnam, "dimer_BMHT");
    else if (ich == 3)
        strcpy(mfilnam, "dimer_SPSS");
    else if (ich == 4)
        strcpy(mfilnam, "xray_AB");
    else {  /* User defined model */
        printf("Name of your model data file: ");
        read_stdin_str(mfilnam);
    }

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(str, env_add);
    strcat(str, "/BendModel/");
    strcat(str, mfilnam);

    fm = open_file(str, "r");

    /* Read in the model data file. The 10 steps must be in order! */
    dmstep = dmatrix(1, 16, 1, 6);
    tmppar = dvector(1, 6);

    zero_dmatrix(dmstep, 16, 6);

    /* Skip the first 3 title lines */
    for (i = 1; i <= 3; i++)
        fgets(str, BUF512 - 1, fm);

    for (i = 1; i <= 10; i++) {
        fscanf(fm, "%s", bpstep);
        upper(bpstep, strlen(bpstep));
        if (strcmp(bpstep, dncl[i]))
            nrerror("Base step in wrong order!!!");
        for (j = 1; j <= 6; j++)
            fscanf(fm, "%lf", &dmstep[i][j]);
    }
    fclose(fm);

    /* Get the step parameters for the complementary steps
       Keep the sign of "Slide, Rise, Roll & Twist" unchanged
       Reverse the sign of "Shift & Tilt" */
    for (i = 11; i <= 16; i++) {
        if (i <= 14)
            j = i - 10;
        else if (i == 15)
            j = 6;
        else
            j = 9;
        dmstep[i][2] = dmstep[j][2];
        dmstep[i][3] = dmstep[j][3];
        dmstep[i][5] = dmstep[j][5];
        dmstep[i][6] = dmstep[j][6];

        dmstep[i][1] = -dmstep[j][1];
        dmstep[i][4] = -dmstep[j][4];
    }

    /* Read in the base sequence */
    get_sequence(bpseq, nbp);

    for (i = 1; i <= 6; i++)
        steppar[1][i] = 0.0;  /* Null operation */

    for (i = 2; i <= (*nbp); i++) {
        strncpy(bpstep, (bpseq[1] + (i - 2)), 2);  /* Note bpseq start from 0 not 1 */
        for (j = 1; j <= 16; j++)
            if (!strcmp(dncl[j], bpstep))
                break;
        for (k = 1; k <= 6; k++)
            steppar[i][k] = dmstep[j][k];
    }

    /* Write out base-sequence and step parameters */
    wrt_seq_CEHS(*nbp, bpseq, steppar);

    free_cvector(str, 1, BUF512);
    free_dvector(tmppar, 1, 6);
    free_dmatrix(dmstep, 1, 16, 1, 6);
}

/* TRIMER_CEHS build arbitrary sequence DNA structure using
   %            (a)  Satchwell, Drew & Travers's trimer model
   %            (b)  Brukner et al's trinucleotide model
   %
   % Note: The bending parameter correlates with "Roll" angle */
void trimer_CEHS(long *nbp, char **bpseq, double **steppar)
{
    long i, j, ich;
    char *str, *env_add, mfilnam[BUF512], **tbp;
    char tbp1[4], tbp2[4], ctbp2[4], tmp_tbp[8];
    FILE *fid;
    double *troll, low_roll, up_roll, max_roll, min_roll;
    double dval;

    str = cvector(1, BUF512);

    printf("Choose from one of the following trimer models\n");
    printf("1. Satchwell et al's model\n");
    printf("2. Brukner et al's model\n");
    printf("3. User defined model\n");
    printf("Your choice (1-3, Dft 1): ");
    read_stdin_str(str);
    ich = atoi(str);

    if ((ich < 1) || (ich > 3))
        ich = 1;

    if (ich == 1)
        strcpy(mfilnam, "trimer_SDT");
    else if (ich == 2)
        strcpy(mfilnam, "trimer_BSSP");
    else {  /* User defined model */
        printf("Name of your model data file: ");
        read_stdin_str(mfilnam);
    }

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(str, env_add);
    strcat(str, "/BendModel/");
    strcat(str, mfilnam);

    fid = open_file(str, "r");

    /* Read in the model data file */
    for (i = 1; i <= 3; i++)
        fgets(str, BUF512 - 1, fid);

    tbp = cmatrix(1, 32, 1, 8);
    troll = dvector(1, 32);

    for (i = 1; i <= 32; i++) {
        fscanf(fid, "%s", tmp_tbp);
        upper(tmp_tbp, strlen(tmp_tbp));
        for (j = 0; j <= 2; j++) {
            tbp1[j] = tmp_tbp[j];
            tbp2[j] = tmp_tbp[6 - j];
        }
        comp_base(ctbp2, tbp2, 3);
        if (strcmp(tbp1, ctbp2))
            nrerror("There are mismatched bases!");
        strcpy(tbp[i], tmp_tbp);

        fscanf(fid, "%lf", &troll[i]);
    }

    fclose(fid);

    /* Change bending parameters to roll angle */
    printf("Minimum roll (Dft 0.0): ");
    read_stdin_str(str);
    if (strlen(str) == 0)
        low_roll = 0.0;
    else
        sscanf(str, "%lf", &low_roll);

    printf("Maximum roll (Dft: 10.0): ");
    read_stdin_str(str);
    if (strlen(str) == 0)
        up_roll = 10.0;
    else
        sscanf(str, "%lf", &up_roll);

    if (up_roll < low_roll)
        nrerror("Error roll angle range!");

    max_roll = max_dvector(troll, 32);
    min_roll = min_dvector(troll, 32);

    /* Scale the roll into [low_roll...up_roll] */
    for (i = 1; i <= 32; i++) {
        dval = (troll[i] - min_roll) / (max_roll - min_roll);
        troll[i] = (up_roll - low_roll) * dval + low_roll;
    }

    /* Read in the base sequence */
    get_sequence(bpseq, nbp);

    /* % The way to get a trinucleotide unit: e.g.
       %    1 2 3 4
       %  A-G-T-C-A-A
       % (1): A-G-T; (2): G-T-C; (3): T-C-A; (4): C-A-A
       % So there are nbp-2 trinucleotide units! */

    /* Deduce the trinucleotide step parameters from base sequence */
    (*nbp)--;
    zero_dmatrix(steppar, *nbp, 6);
    for (i = 2; i <= (*nbp); i++) {
        steppar[i][3] = 3.34;  /* Rise */
        steppar[i][6] = 34.3;  /* Twist */
    }

    for (i = 1; i <= (*nbp) - 1; i++) {
        strncpy(tbp1, (bpseq[1] + (i - 1)), 3);
        for (j = 1; j <= 32; j++)
            if (strstr(tbp[j], tbp1) != NULL)
                break;
        steppar[i + 1][5] = troll[j];
    }

    /* Write out base-sequence and step parameters */
    wrt_seq_CEHS(*nbp, bpseq, steppar);

    free_cvector(str, 1, BUF512);
    free_dvector(troll, 1, 32);
    free_cmatrix(tbp, 1, 32, 1, 8);
}

/* STEP_CEHS build arbitrary sequence DNA structure using the
   %          CEHS step parameters supplied.
   %
   %         The input file can be most easily generated by firstly
   %         running a proper dimer/trimer model, and then modifying
   %         the output parameter file as desired. */
void step_CEHS(long *nbp, char **bpseq, double **steppar, int itot, char **ist, double **xconf, int phi, char **output) //char **output put this in argument not to go over output file output_coordinates.dat
{
    long i, j;
    char str[BUF512], tmpstr[4];
    char seq1;	
//    FILE *fp;

/* ----- Commented to probe for efficiency ----------- */



//sprintf(str, "output_dnaflex/output_coordinates_%.6d.dat", phi);

//strcpy(str, "output_coordinates.dat");

/*
    // Read in the step parameters
    printf("Sequence-CEHS step parameter file name (Dft output_coordinates.dat): ");
    read_stdin_str(str);

    if (strlen(str) == 0)
        strcpy(str, "output_coordinates.dat");
*/

*nbp = (long)itot + 1;

    for (i = 0; i < (*nbp); i++) {
        
        sscanf(output[i+2], "%s", tmpstr);
        bpseq[1][i] = toupper(tmpstr[0]);
        bpseq[2][i] = toupper(tmpstr[2]);

        for (j = 0; j < 6; j++){
            steppar[i+1][j+1] = xconf[j][i-1];	// hier wird schon was reingeschrieben obwohl es noch nicht aufs array xconf zugreift
				}

				  }


        sscanf(output[*nbp+2], "%s", tmpstr);
        bpseq[1][*nbp-1] = toupper(tmpstr[0]);
        bpseq[2][*nbp-1] = toupper(tmpstr[2]);

/* This is the code if we want to read from the file

      fp = open_file(str, "r");

      fgets(str, BUF512 - 1, fp);
      sscanf(str, "%ld", nbp);
      fgets(str, BUF512 - 1, fp);  // skip one line

    for (i = 1; i <= (*nbp); i++) {
        fscanf(fp, "%s", tmpstr);
        bpseq[1][i - 1] = toupper(tmpstr[0]);
        bpseq[2][i - 1] = toupper(tmpstr[2]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &steppar[i][j]);
    }
*/

/*
printf ("%ld \n", *nbp);
   for (i = 0; i < (*nbp); i++) {
printf("%c \n", bpseq[2][i]);
}
printf("%lf \n", steppar[2][1]);
printf("%s \n", output[3]);
*/





// DO NOT DO THIS BECAUSE IT GIVES AN ERROR WHEN EXECUTING THE WHOLE MC STUFF 1000 TIMES

    /* Base-pair consistency check! */
/*
    comp_base(str, bpseq[1], *nbp);
    if (strcmp(str, bpseq[2]))
        nrerror("There are mismatched bases!");
*/
//   fclose(fp);
}

void exact_CEHS(long *nbp, char **bpseq, double **bppar, double **steppar)
{
    long i, j;
    char str[BUF512], tmpstr[4];
    FILE *fp;

    printf("Name of the local CEHS parameter data file: ");
    read_stdin_str(str);
    fp = open_file(str, "r");

    fgets(str, BUF512 - 1, fp);
    sscanf(str, "%ld", nbp);
    fgets(str, BUF512 - 1, fp);  /* skip one line */

    for (i = 1; i <= (*nbp); i++) {
        fscanf(fp, "%s", tmpstr);
        bpseq[1][i - 1] = toupper(tmpstr[0]);
        bpseq[2][i - 1] = toupper(tmpstr[2]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &bppar[i][j]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &steppar[i][j]);
    }

    fclose(fp);
}

/* % ATOMIC_CEHS build DNA atomic structure using CEHS parameters
   %
   % Input: nbp--number of base-pairs
   %        bpseq--base-pair sequence
   %        bppar--[nbp by 6] CEHS base-pair parameters
   %        steppar--[nbp by 6] CEHS step parameters
   %        xdir--+x-axis direction */
void atomic_CEHS(long nbp, char **bpseq, double **bppar, double **steppar,
                 long xdir, char *filnam, double **xyz_float, char **pdb_data, char **pdb_tot, char **n7_data, int *tot_num_at)
{
    long i, j, k, m, i2, num1, num2, tn1 = 0, tn2 = 0;
    char base_set[BUF512], *dir_name, *env_add;
    double **orien_next, **orien_next_T, *pos_next, **org_xyz;
    char **asym1, **asym2, **t_asym1, **t_asym2;
    double **xyz1, **xyz2, **t_xyz1, **t_xyz2, **xyztmp1, **xyztmp2;
    char *t_btype1, *t_btype2;
    long *t_bnum1, *t_bnum2, ATM_NUM = 100;
    char bpname[3];
    long **ts2;  /* keep track of strand II base atoms */

    double **orien, **mst, *pos, *tmp_vec;

    long *base_num, tot_num;
    char **tot_asym, *base_type, *strand;
    double **tot_xyz, **tot_xyz2;

    /* Read in atomic base geometry */
    dir_name = cvector(1, BUF512);


/* ----- Commented to probe for efficiency ----------- */

strcpy(base_set, "NDB96");

/*
    printf("\nBase geometry set to use (Dft NDB96): ");
    read_stdin_str(base_set);

    if (strlen(base_set) == 0)
        strcpy(base_set, "NDB96");
*/

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(dir_name, env_add);
    strcat(dir_name, "/BaseGeo/");
    strcat(dir_name, base_set);
    strcat(dir_name, "_");

    /* Variable initializations */
    orien_next = dmatrix(1, 3, 1, 3);
    orien_next_T = dmatrix(1, 3, 1, 3);
    pos_next = dvector(1, 3);
    org_xyz = dmatrix(1, nbp, 1, 3);

    identity_dmatrix(orien_next, 3);  /* for setting the first bp */
    zero_dvector(pos_next, 3);

    /* For strand I */
    t_asym1 = cmatrix(1, MAX_ATOM, 1, 4);
    t_btype1 = cvector(1, MAX_ATOM);
    t_bnum1 = lvector(1, MAX_ATOM);
    t_xyz1 = dmatrix(1, MAX_ATOM, 1, 3);

    /* For strand II */
    t_asym2 = cmatrix(1, MAX_ATOM, 1, 4);
    t_btype2 = cvector(1, MAX_ATOM);
    t_bnum2 = lvector(1, MAX_ATOM);
    t_xyz2 = dmatrix(1, MAX_ATOM, 1, 3);

    orien = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    pos = dvector(1, 3);
    tmp_vec = dvector(1, 3);

    /* column 1: starting # for each base */
    /* column 2: total atom # for each base */
    ts2 = lmatrix(1, nbp, 1, 2);

    /* Get the whole data set by iterations */
    for (i = 1; i <= nbp; i++) {
        asym1 = cmatrix(1, ATM_NUM, 1, 4);
        asym2 = cmatrix(1, ATM_NUM, 1, 4);
        xyz1 = dmatrix(1, ATM_NUM, 1, 3);
        xyz2 = dmatrix(1, ATM_NUM, 1, 3);
        xyztmp1 = dmatrix(1, ATM_NUM, 1, 3);
        xyztmp2 = dmatrix(1, ATM_NUM, 1, 3);

        bpname[0] = bpseq[1][i - 1];
        bpname[1] = bpseq[2][i - 1];
        bpname[2] = '\0';

        set_bp(&num1, &num2, asym1, asym2, xyz1, xyz2, bpname, bppar[i], xdir, dir_name);

        stepfunc(orien, mst, pos, steppar[i]);
        tr_dmatrix(orien_next_T, orien_next, 3, 3);
        dvec_x_dmtx(tmp_vec, pos, orien_next_T, 3, 3, 3);
        for (j = 1; j <= 3; j++)
            pos_next[j] += tmp_vec[j];

        duplicate_dmtx(mst, orien_next, 3, 3);
        mul_dmatrix(orien_next, mst, orien, 3, 3, 3, 3);

        for (j = 1; j <= 3; j++)
            org_xyz[i][j] = pos_next[j];

        tr_dmatrix(orien_next_T, orien_next, 3, 3);  /* New transpose!!! */
        mul_dmatrix(xyztmp1, xyz1, orien_next_T, num1, 3, 3, 3);
        for (j = 1; j <= num1; j++) {
            k = tn1 + j;
            for (m = 1; m <= 3; m++)
                t_xyz1[k][m] = xyztmp1[j][m] + pos_next[m];
            strcpy(t_asym1[k], asym1[j]);
            t_btype1[k] = bpname[0];
            t_bnum1[k] = i;
        }

        ts2[i][1] = tn2 + 1;
        ts2[i][2] = num2;

        mul_dmatrix(xyztmp2, xyz2, orien_next_T, num2, 3, 3, 3);
        for (j = 1; j <= num2; j++) {
            k = tn2 + j;
            for (m = 1; m <= 3; m++)
                t_xyz2[k][m] = xyztmp2[j][m] + pos_next[m];
            strcpy(t_asym2[k], asym2[j]);
            t_btype2[k] = bpname[1];
            t_bnum2[k] = 2 * nbp - i + 1;
        }
        tn1 += num1;
        tn2 += num2;

        free_cmatrix(asym1, 1, ATM_NUM, 1, 4);
        free_cmatrix(asym2, 1, ATM_NUM, 1, 4);
        free_dmatrix(xyz1, 1, ATM_NUM, 1, 3);
        free_dmatrix(xyz2, 1, ATM_NUM, 1, 3);
        free_dmatrix(xyztmp1, 1, ATM_NUM, 1, 3);
        free_dmatrix(xyztmp2, 1, ATM_NUM, 1, 3);
    }

    tot_num = tn1 + tn2;

    tot_xyz = dmatrix(1, tot_num, 1, 3);
    tot_xyz2 = dmatrix(1, tot_num, 1, 3);
    tot_asym = cmatrix(1, tot_num, 1, 4);
    base_type = cvector(1, tot_num);
    strand = cvector(1, tot_num);
    base_num = lvector(1, tot_num);

    /* Combine two strands together! */
    for (i = 1; i <= tn1; i++) {
        strcpy(tot_asym[i], t_asym1[i]);
        base_type[i] = t_btype1[i];
        strand[i] = 'A';
        base_num[i] = t_bnum1[i];
        for (j = 1; j <= 3; j++)
            tot_xyz[i][j] = t_xyz1[i][j];
    }
    k = tn1;  /* total atom number in strand I */

    /* The order of strand II is now from 3'--->5', this need
       to be reversed according to information in <ts2> */
    for (i = nbp; i >= 1; i--) {
        for (j = ts2[i][1]; j <= ts2[i][1] + ts2[i][2] - 1; j++) {
            m = k + j - ts2[i][1] + 1;
            strcpy(tot_asym[m], t_asym2[j]);
            base_type[m] = t_btype2[j];
            strand[m] = 'B';
            base_num[m] = t_bnum2[j];
            for (i2 = 1; i2 <= 3; i2++)
                tot_xyz[m][i2] = t_xyz2[j][i2];
        }
        k += ts2[i][2];
    }

    /* Align the helix */
    align_helix(tot_xyz2, tot_xyz, tot_num, org_xyz, nbp);

    /* Write out the result */

	for (i = 1; i <= 2*nbp; i++)
	{
	  for (j = 1; j <= tot_num; j++)
	  {
		if(base_num[j] == i){
			if(base_type[j] == 'C' || base_type[j] == 'T'){	
				if(strcmp(tot_asym[j],"C6 ") == 0){ 
							
							sprintf(pdb_data[i], "%3s   %c %c %4ld %8.3f%8.3f%8.3f\n", tot_asym[j],base_type[j], strand[j], base_num[j], tot_xyz2[j][1], tot_xyz2[j][2], tot_xyz2[j][3]);
							     }
									}

			else if(base_type[j] == 'A' || base_type[j] == 'G'){	
				      if(strcmp(tot_asym[j],"C8 ") == 0){ 
							sprintf(pdb_data[i], "%3s   %c %c %4ld %8.3f%8.3f%8.3f\n", tot_asym[j],base_type[j], strand[j], base_num[j], tot_xyz2[j][1], tot_xyz2[j][2], tot_xyz2[j][3]);
								     }
			     						     }
				    }
				
	   }
	}
		

j=0;	
for (i = 1; i <= tot_num; i++){
		if(strcmp(tot_asym[i],"N7 ") == 0){ 
						sprintf(n7_data[j], "%4ld %8.3f%8.3f%8.3f\n", base_num[i], tot_xyz2[i][1], tot_xyz2[i][2], tot_xyz2[i][3]);
						  j++;
						  }
			       }


for (i = 1; i <= tot_num; i++){
        	sprintf(pdb_tot[i], "%3s  D%c %c%4ld    %8.3f%8.3f%8.3f\n", tot_asym[i],base_type[i], strand[i], base_num[i], tot_xyz2[i][1], tot_xyz2[i][2], tot_xyz2[i][3]);
				}

	*tot_num_at = (int) tot_num;

	/*for(i=2*nbp-10;i<=2*nbp;i++){
		printf("%s",pdb_tot[i]);
					}*/
	//"ATOM  %5ld  %3s  D%c %c%4ld    %8.3f%8.3f%8.3f\n"
    //wrt_my_pdb(tot_num, tot_asym, base_type, strand, base_num, tot_xyz2, filnam, xyz_float);

    /* free memory */
    free_cvector(dir_name, 1, BUF512);
    free_dvector(pos_next, 1, 3);
    free_cvector(t_btype1, 1, MAX_ATOM);
    free_lvector(t_bnum1, 1, MAX_ATOM);
    free_cvector(t_btype2, 1, MAX_ATOM);
    free_lvector(t_bnum2, 1, MAX_ATOM);
    free_dvector(pos, 1, 3);
    free_dvector(tmp_vec, 1, 3);
    free_cvector(base_type, 1, tot_num);
    free_cvector(strand, 1, tot_num);
    free_lvector(base_num, 1, tot_num);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_dmatrix(orien_next_T, 1, 3, 1, 3);
    free_dmatrix(org_xyz, 1, nbp, 1, 3);
    free_cmatrix(t_asym1, 1, MAX_ATOM, 1, 4);
    free_dmatrix(t_xyz1, 1, MAX_ATOM, 1, 3);
    free_cmatrix(t_asym2, 1, MAX_ATOM, 1, 4);
    free_dmatrix(t_xyz2, 1, MAX_ATOM, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_lmatrix(ts2, 1, nbp, 1, 2);
    free_dmatrix(tot_xyz, 1, tot_num, 1, 3);
    free_dmatrix(tot_xyz2, 1, tot_num, 1, 3);
    free_cmatrix(tot_asym, 1, tot_num, 1, 4);
}

/* % BLOCK_CEHS build DNA schematic structure using CEHS parameters
   %
   % Input: nblk--number of regular blocks
   %        steppar--[nblk by 6] CEHS step parameters
   %        xdir--+x-axis direction  */
void block_CEHS(long nblk, double **steppar, long xdir, char *filnam)
{
    long i, j, k, m, vnum, nbond, **lkg, NBLK = 100;
    char blknam[BUF512], *dir_name, *env_add, **asym, **tot_asym;
    double **blkxyz, **orien_next, *pos_next, **org_xyz;
    double **tot_xyz, **tot_xyz2;
    long tot_vnum, tot_nbond, **tot_lkg;
    long num_vtx = 0, num_lkg = 0;

    double **orien_next_T, **orien, **mst, *pos;
    double *tmp_vec, **blkxyz2;

    /* Read in rigid block geometry */
    dir_name = cvector(1, BUF512);

    printf("\nGeometry data file name for the block (Dft BLOCK.alc): ");
    read_stdin_str(blknam);

    if (strlen(blknam) == 0)
        strcpy(blknam, "BLOCK.alc");

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(dir_name, env_add);
    strcat(dir_name, "/BaseGeo/");
    strcat(dir_name, blknam);

    asym = cmatrix(1, NBLK, 1, 3);  /* 2 characters long */
    lkg = lmatrix(1, NBLK, 1, 2);
    blkxyz = dmatrix(1, NBLK, 1, 3);

    rdalc(&vnum, &nbond, asym, lkg, blkxyz, dir_name);

    /* Set block orientation properly */
    if (xdir == 1)
        for (i = 1; i <= vnum; i++) {
            blkxyz[i][1] = -blkxyz[i][1];
            blkxyz[i][3] = -blkxyz[i][3];
        }

    tot_vnum = nblk * vnum;
    tot_nbond = nblk * nbond;

    orien_next = dmatrix(1, 3, 1, 3);
    orien_next_T = dmatrix(1, 3, 1, 3);
    orien = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);
    pos_next = dvector(1, 3);
    pos = dvector(1, 3);
    tmp_vec = dvector(1, 3);
    org_xyz = dmatrix(1, nblk, 1, 3);
    tot_lkg = lmatrix(1, tot_nbond, 1, 2);
    tot_asym = cmatrix(1, tot_vnum, 1, 3);
    tot_xyz = dmatrix(1, tot_vnum, 1, 3);
    tot_xyz2 = dmatrix(1, tot_vnum, 1, 3);

    identity_dmatrix(orien_next, 3);
    zero_dvector(pos_next, 3);
    blkxyz2 = dmatrix(1, NBLK, 1, 3);

    /* Get the whole data by iterations */
    for (i = 1; i <= nblk; i++) {
        stepfunc(orien, mst, pos, steppar[i]);
        tr_dmatrix(orien_next_T, orien_next, 3, 3);
        dvec_x_dmtx(tmp_vec, pos, orien_next_T, 3, 3, 3);
        for (j = 1; j <= 3; j++)
            pos_next[j] += tmp_vec[j];

        duplicate_dmtx(mst, orien_next, 3, 3);
        mul_dmatrix(orien_next, mst, orien, 3, 3, 3, 3);

        for (j = 1; j <= 3; j++)
            org_xyz[i][j] = pos_next[j];

        tr_dmatrix(orien_next_T, orien_next, 3, 3);  /* New transpose!!! */
        mul_dmatrix(blkxyz2, blkxyz, orien_next_T, vnum, 3, 3, 3);
        for (j = 1; j <= vnum; j++) {
            k = num_vtx + j;
            strcpy(tot_asym[k], asym[j]);
            for (m = 1; m <= 3; m++)
                tot_xyz[k][m] = blkxyz2[j][m] + pos_next[m];
        }

        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tot_lkg[num_lkg + j][k] = lkg[j][k] + (i - 1) * vnum;

        num_vtx += vnum;
        num_lkg += nbond;
    }

    /* Align the helix */
    align_helix(tot_xyz2, tot_xyz, tot_vnum, org_xyz, nblk);

    /* Write out the result */
    wrtalc(tot_vnum, tot_nbond, tot_asym, tot_lkg, tot_xyz2, filnam);

    /* free memory */
    free_cvector(dir_name, 1, BUF512);
    free_dvector(pos_next, 1, 3);
    free_dvector(pos, 1, 3);
    free_dvector(tmp_vec, 1, 3);
    free_cmatrix(asym, 1, NBLK, 1, 3);
    free_lmatrix(lkg, 1, NBLK, 1, 2);
    free_dmatrix(blkxyz, 1, NBLK, 1, 3);
    free_dmatrix(orien_next, 1, 3, 1, 3);
    free_dmatrix(orien_next_T, 1, 3, 1, 3);
    free_dmatrix(orien, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(org_xyz, 1, nblk, 1, 3);
    free_cmatrix(tot_asym, 1, tot_vnum, 1, 3);
    free_lmatrix(tot_lkg, 1, tot_nbond, 1, 2);
    free_dmatrix(tot_xyz, 1, tot_vnum, 1, 3);
    free_dmatrix(tot_xyz2, 1, tot_vnum, 1, 3);
    free_dmatrix(blkxyz2, 1, NBLK, 1, 3);
}

/* CEHS_BUILD is a double helical structure building program using local
   CEHS parameters. Program written in C by Xiang-Jun Lu (1996--1997) */
void CEHS_build(int itot, char **ist, double **xconf, int phi, char **output, double **xyz_float, char **pdb_data, char **pdb_tot, char **n7_data, int *tot_num_at)
{
    long ich, imdl, xdir, nbp, ipos;
    char str[BUF512], filnam[BUF512], **bpseq, filcnt[BUF512], *ptr;
    double *ave_par, **bppar, **steppar;
    

/* ----- Commented to probe for efficiency ----------- */
ich = 4;
/*
    system("clear");
    printf("Structure rebuilding based on CEHS parameters\n\n");
    printf("1. Uniform regular helix\n");
    printf("2. Step parameters from dinucleotide bending models\n");
    printf("3. Parameters from trinucleotide bending models\n");
    printf("4. Step parameters from a data file\n");
    printf("5. Step and base-pair parameters from a data file\n");
    printf("6. Quit\n");
    printf("Your choice(1-6, Dft 5): ");
    read_stdin_str(str);
    ich = atoi(str);
    if ((ich < 1) || (ich > 6))
        ich = 5;
*/



    bpseq = cmatrix(1, 2, 1, MAX_BP);  /* NOTE: bpseq is [2 by MAX_BP] */
    bppar = dmatrix(1, MAX_BP, 1, 6);
    steppar = dmatrix(1, MAX_BP, 1, 6);

    if (ich == 1)
        std_CEHS(&nbp, bpseq, steppar);
    else if (ich == 2)
        dimer_CEHS(&nbp, bpseq, steppar);
    else if (ich == 3)
        trimer_CEHS(&nbp, bpseq, steppar);
    else if (ich == 4)
        step_CEHS(&nbp, bpseq, steppar, itot, ist, xconf, phi, output);
    else if (ich == 5)
        exact_CEHS(&nbp, bpseq, bppar, steppar);
    else
        exit(0);

    ave_par = dvector(1, 6);
    /* The first row of <twist> is all ZEROs */
    mean_dmatrix(ave_par, steppar, nbp, 6);
    xdir = get_xdir(ave_par[6]);  /* Set +x-axis direction */






/* -------- Commented to probe for efficiency ----------- */

sprintf(filnam, "output_schnarp/structure_%.6d.pdb", phi);

//strcpy(filnam, "CEHS_gen.out");

/*
    printf("\nOutput data file name (Dft CEHS_gen.out): ");
    read_stdin_str(filnam);
    if (strlen(filnam) == 0)
        strcpy(filnam, "CEHS_gen.out");
*/





/* ----- Commented to probe for efficiency ----------- */

imdl = 1;

/*
    printf("\nChoose from one of the following representations:\n");
    printf("1. Real atomic model (in PDB format)\n");
    printf("2. Rectangular block model (in ALCHEMY format)\n");
    printf("Your choice (1 or 2, Dft 1): ");
    read_stdin_str(str);
    imdl = atoi(str);
*/
    if ((imdl < 1) || (imdl > 2))
        imdl = 1;

    if (imdl == 1) {  /* Use atomic model */
        if (ich != 5)
            set_bppar(bppar, nbp, bpseq);  /* Get base-pair parameters */
        atomic_CEHS(nbp, bpseq, bppar, steppar, xdir, filnam, xyz_float,pdb_data,pdb_tot,n7_data,tot_num_at);
    } else {  /* Use rectangular block model */
        block_CEHS(nbp, steppar, xdir, filnam);
        printf("\n\aFile %s has NO backbone connection information", filnam);

        ptr = strchr(filnam, '.');
        if (ptr != NULL) {
            ipos = ptr - filnam + 1;
            strncpy(filcnt, filnam, (size_t) ipos);
            filcnt[ipos] = '\0';  /* Just to make sure .... */
        } else {
            strcpy(filcnt, filnam);
            strcat(filcnt, ".\0");
        }
        strcat(filcnt, "cnt");
        alc_connect(filnam, filcnt);
        printf("\n\aFile %s HAS backbone connection information\n", filcnt);

        strcpy(filnam, filcnt);  /* For display! */
    }

/* ----- Commented to probe for efficiency ----------- */

/* We will never need to display the structure with RasMol*/

/*
   str_display(imdl, filnam); // To display the structure with RasMol
*/

    free_dvector(ave_par, 1, 6);
    free_cmatrix(bpseq, 1, 2, 1, MAX_BP);
    free_dmatrix(bppar, 1, MAX_BP, 1, 6);
    free_dmatrix(steppar, 1, MAX_BP, 1, 6);
}

/* %%%%%%%%%%%%%%%%%%% CEHS Building functions END %%%%%%%%%%%%%%%%%%% */

/* %%%%%%%%%%%%%%%%%%% GLH Building functions BEGIN %%%%%%%%%%%%%%%%%%% */

/*  WRT_SEQ_GLH write out the sequence and its associated GLOBAL parameters
    %             This file can be used for error check or fined-turned
    %             to get the desired result.
    %
    % Input: nbp--number of base-pairs
    %        bpseq--[2 by nbp] base-pair sequence
    %        glhpar--[nbp by 6] GLOBAL parameters         */
void wrt_seq_GLH(long nbp, char **bpseq, double **glhpar)
{
    long i, j;
    char fstep[BUF512];
    FILE *fo;

    printf("\nOutput sequence GLH parameter data file name (Dft GLH_seq.dat): ");
    read_stdin_str(fstep);

    if (strlen(fstep) == 0)
        strcpy(fstep, "GLH_seq.dat");

    fo = open_file(fstep, "w");

    fprintf(fo, "%5ld base-pairs\n", nbp);
    fprintf(fo, "       X-dsp   Y-dsp    Rise   Incl.    Tip    Twist\n");
    for (i = 1; i <= nbp; i++) {
        fprintf(fo, "%c-%c ", bpseq[1][i - 1], bpseq[2][i - 1]);
        for (j = 1; j <= 6; j++)
            fprintf(fo, "%8.2f", glhpar[i][j]);
        fprintf(fo, "\n");
    }
    fclose(fo);
}

/* STD_GLH build A, B, C, Z and other regular double helical structures
   using a set of 6 GLOBAL parameters */
void std_GLH(long *nbp, char **bpseq, double **glhpar)
{
    long i, j, ich;
    char str[BUF512];
    double *glhblk, *Z_CG, *Z_GC;
    double d1, d2, d3, d4, d5, d6;

    printf("1. Standard A-form DNA (Twist=33, Rise=2.6, Incl=19, X_disp=-4.5)\n");
    printf("2. Standard B-form DNA (Twist=36, Rise=3.4, Incl=-6)\n");
    printf("3. Standard C-form DNA (Twist=39, Rise=3.3, Incl=-8  X_disp=+0.9)\n");
    printf("4. Standard Z-form DNA (Twist=-11, Y_disp=-2.2 for CG step)\n");
    printf("                        Twist=-49, Y_disp=+2.2 for GC step)\n");
    printf("               Rise=3.7, Incl=6 and X_disp=-3 for all steps\n");
    printf("5. Other regular double helical structures\n");
    printf("Input your choice (1-5, Dft 2): ");
    read_stdin_str(str);
    ich = atoi(str);
    if ((ich < 1) || (ich > 5))
        ich = 2;

    glhblk = dvector(1, 6);
    Z_CG = dvector(1, 6);
    Z_GC = dvector(1, 6);

    /* Get the GLOBAL parameters for building blocks, in the order of
       X-displacement, Y-displacement, Rise, Inclination, Tip & Twist */
    switch (ich) {
    case 1:  /* A-form: [-4.5 0 2.6 19 0 33] */
        glhblk[1] = -4.5;
        glhblk[2] = 0.0;
        glhblk[3] = 2.6;
        glhblk[4] = 19.0;
        glhblk[5] = 0.0;
        glhblk[6] = 33.0;
        break;
    case 2:  /* B-form: [0 0 3.4 -6 0 36] */
        glhblk[1] = 0.0;
        glhblk[2] = 0.0;
        glhblk[3] = 3.4;
        glhblk[4] = -6.0;
        glhblk[5] = 0.0;
        glhblk[6] = 36.0;
        break;
    case 3:  /* C-form: [0.9 0 3.3 -8 0 39] */
        glhblk[1] = 0.9;
        glhblk[2] = 0.0;
        glhblk[3] = 3.3;
        glhblk[4] = -8.0;
        glhblk[5] = 0.0;
        glhblk[6] = 39.0;
        break;
    case 4:
        /* Z-form: Note its specialty (sign for x-displacement) */
        /* C-G base-pair [-3 -2.2 3.7 6 0 -11] */
        Z_CG[1] = -3.0;
        Z_CG[2] = -2.2;
        Z_CG[3] = 3.7;
        Z_CG[4] = 6.0;
        Z_CG[5] = 0.0;
        Z_CG[6] = -11.0;
        /* G-C base-pair [-3  2.2 3.7 6 0 -49] */
        Z_GC[1] = -3.0;
        Z_GC[2] = 2.2;
        Z_GC[3] = 3.7;
        Z_GC[4] = 6.0;
        Z_GC[5] = 0.0;
        Z_GC[6] = -49.0;
        break;
    default:  /* Other regular form */
        printf("\nInput 6 GLOBAL parameters in the order of: \n");
        printf("X-disp  Y-disp  Rise  Inclination  Tip  Twist\n");
        while (1) {  /* Infinite loop */
            read_stdin_str(str);
            j = sscanf(str, "%lf%lf%lf%lf%lf%lf", &d1, &d2, &d3, &d4, &d5, &d6);
            if (j == 6)
                break;
            printf("Six Global-parameters needed. Try again!\n");
        }
        glhblk[1] = d1;
        glhblk[2] = d2;
        glhblk[3] = d3;
        glhblk[4] = d4;
        glhblk[5] = d5;
        glhblk[6] = d6;
    }

    /* Read in base sequence */
    get_sequence(bpseq, nbp);

    /* Get the whole set of GLOBAL parameters */
    if (ich != 4)  /* Not Z-DNA */
        for (i = 1; i <= (*nbp); i++)
            for (j = 1; j <= 6; j++)
                glhpar[i][j] = glhblk[j];
    else {
        for (i = 1; i <= (*nbp); i += 2)  /* for C-G base-pair */
            for (j = 1; j <= 6; j++)
                glhpar[i][j] = Z_CG[j];
        for (i = 2; i <= (*nbp); i += 2)  /* for G-C base-pair */
            for (j = 1; j <= 6; j++)
                glhpar[i][j] = Z_GC[j];
    }

    /* Adjust rise & twist */
    for (i = (*nbp); i >= 2; i--) {
        glhpar[i][3] = glhpar[i - 1][3];  /* Rise */
        glhpar[i][6] = glhpar[i - 1][6];  /* Twist */
    }

    /* Rise & twist for the 1st base-pair */
    glhpar[1][3] = 0.0;
    glhpar[1][6] = 0.0;

    /* Write out base-pair sequence & global parameters */
    wrt_seq_GLH(*nbp, bpseq, glhpar);

    free_dvector(glhblk, 1, 6);
    free_dvector(Z_CG, 1, 6);
    free_dvector(Z_GC, 1, 6);
}

/*  STEP_GLH build arbitrary sequence DNA structure using the
    GLOBAL parameters supplied. */
void step_GLH(long *nbp, char **bpseq, double **glhpar)
{
    long i, j;
    char str[BUF512], tmpstr[4];
    FILE *fp;

    /* Read in the GLOBAL parameters */
    printf("Sequence-GLH parameter file name (Dft GLH_seq.dat): ");
    read_stdin_str(str);

    if (strlen(str) == 0)
        strcpy(str, "GLH_seq.dat");

    fp = open_file(str, "r");

    fgets(str, BUF512 - 1, fp);
    sscanf(str, "%ld", nbp);
    fgets(str, BUF512 - 1, fp);  /* skip one line */

    for (i = 1; i <= (*nbp); i++) {
        fscanf(fp, "%s", tmpstr);
        bpseq[1][i - 1] = toupper(tmpstr[0]);
        bpseq[2][i - 1] = toupper(tmpstr[2]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &glhpar[i][j]);
    }

    /* Base-pair consistency check! */
    comp_base(str, bpseq[1], *nbp);
    if (strcmp(str, bpseq[2]))
        nrerror("There are mismatched bases!");

    fclose(fp);
}

void exact_GLH(long *nbp, char **bpseq, double **bppar, double **glhpar)
{
    long i, j;
    char str[BUF512], tmpstr[4];
    FILE *fp;

    printf("Name of the GLH parameter data file: ");
    read_stdin_str(str);
    fp = open_file(str, "r");

    fgets(str, BUF512 - 1, fp);
    sscanf(str, "%ld", nbp);
    fgets(str, BUF512 - 1, fp);  /* skip one line */

    for (i = 1; i <= (*nbp); i++) {
        fscanf(fp, "%s", tmpstr);
        bpseq[1][i - 1] = toupper(tmpstr[0]);
        bpseq[2][i - 1] = toupper(tmpstr[2]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &bppar[i][j]);
        for (j = 1; j <= 6; j++)
            fscanf(fp, "%lf", &glhpar[i][j]);
    }

    fclose(fp);
}

/* % ATOMIC_GLH build DNA atomic structure using GLOBAL parameters
   %
   % Input: nbp--number of base-pairs
   %        npseq--base-pair sequence
   %        bppar--[nbp by 6] CEHS base-pair parameters
   %        glhpar--[nbp by 6] GLOBAL parameters
   %        xdir--+x-axis direction */
void atomic_GLH(long nbp, char **bpseq, double **bppar, double **glhpar,
                long xdir, char *filnam)
{
    long i, j, k, m, i2, num1, num2, tn1 = 0, tn2 = 0;
    char base_set[BUF512], *dir_name, *env_add;
    double **org_xyz;
    char **asym1, **asym2, **t_asym1, **t_asym2;
    double **xyz1, **xyz2, **t_xyz1, **t_xyz2, **xyztmp1, **xyztmp2;
    char *t_btype1, *t_btype2;
    long *t_bnum1, *t_bnum2, ATM_NUM = 100;
    char bpname[3];
    long **ts2;  /* keep track of strand II base atoms */

    long *base_num, tot_num;
    char **tot_asym, *base_type, *strand;
    double **tot_xyz, **tot_xyz2;

    double twist_g = 0.0, rise_g = 0.0, tip, incl, y_disp, x_disp;
    double **v_orien, **v_orien_T, *bp_pos, *tmp_vec;
    double *hinge, TipInc, **bp_orien, **bp_orien_T;

    /* Read in atomic base geometry */
    dir_name = cvector(1, BUF512);

    printf("\nBase geometry set to use (Dft NDB96): ");
    read_stdin_str(base_set);

    if (strlen(base_set) == 0)
        strcpy(base_set, "NDB96");

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(dir_name, env_add);
    strcat(dir_name, "/BaseGeo/");
    strcat(dir_name, base_set);
    strcat(dir_name, "_");

    /* Variable initializations */
    org_xyz = dmatrix(1, nbp, 1, 3);

    /* For strand I */
    t_asym1 = cmatrix(1, MAX_ATOM, 1, 4);
    t_btype1 = cvector(1, MAX_ATOM);
    t_bnum1 = lvector(1, MAX_ATOM);
    t_xyz1 = dmatrix(1, MAX_ATOM, 1, 3);

    /* For strand II */
    t_asym2 = cmatrix(1, MAX_ATOM, 1, 4);
    t_btype2 = cvector(1, MAX_ATOM);
    t_bnum2 = lvector(1, MAX_ATOM);
    t_xyz2 = dmatrix(1, MAX_ATOM, 1, 3);

    /* column 1: starting # for each base */
    /* column 2: total atom # for each base */
    ts2 = lmatrix(1, nbp, 1, 2);

    v_orien = dmatrix(1, 3, 1, 3);
    v_orien_T = dmatrix(1, 3, 1, 3);
    bp_orien = dmatrix(1, 3, 1, 3);
    bp_orien_T = dmatrix(1, 3, 1, 3);
    bp_pos = dvector(1, 3);
    tmp_vec = dvector(1, 3);
    hinge = dvector(1, 3);

    /* Get the whole data set by iterations */
    for (i = 1; i <= nbp; i++) {
        asym1 = cmatrix(1, ATM_NUM, 1, 4);
        asym2 = cmatrix(1, ATM_NUM, 1, 4);
        xyz1 = dmatrix(1, ATM_NUM, 1, 3);
        xyz2 = dmatrix(1, ATM_NUM, 1, 3);
        xyztmp1 = dmatrix(1, ATM_NUM, 1, 3);
        xyztmp2 = dmatrix(1, ATM_NUM, 1, 3);

        bpname[0] = bpseq[1][i - 1];
        bpname[1] = bpseq[2][i - 1];
        bpname[2] = '\0';

        set_bp(&num1, &num2, asym1, asym2, xyz1, xyz2, bpname, bppar[i], xdir, dir_name);

        /* Decompose the 6 GLOBAL parameters */
        x_disp = glhpar[i][1];
        y_disp = glhpar[i][2];
        rise_g += glhpar[i][3];
        incl = glhpar[i][4];
        tip = glhpar[i][5];
        twist_g += glhpar[i][6];

        rotz(v_orien, twist_g);  /* bp orientation after "twisting" */
        tr_dmatrix(v_orien_T, v_orien, 3, 3);
        tmp_vec[1] = x_disp;
        tmp_vec[2] = y_disp;
        tmp_vec[3] = rise_g;
        dvec_x_dmtx(bp_pos, tmp_vec, v_orien_T, 3, 3, 3);

        /* Each bp's origin in global reference frame */
        for (j = 1; j <= 3; j++)
            org_xyz[i][j] = bp_pos[j];

        for (j = 1; j <= 3; j++)
            hinge[j] = incl * v_orien[j][1] + tip * v_orien[j][2];  /* axis */
        TipInc = dist_ab(incl, tip);  /* angle */
        arbrot(v_orien_T, hinge, TipInc);  /* v_orien_T as rotmat */
        /* Final bp orientation */
        mul_dmatrix(bp_orien, v_orien_T, v_orien, 3, 3, 3, 3);
        tr_dmatrix(bp_orien_T, bp_orien, 3, 3);

        mul_dmatrix(xyztmp1, xyz1, bp_orien_T, num1, 3, 3, 3);
        for (j = 1; j <= num1; j++) {
            k = tn1 + j;
            for (m = 1; m <= 3; m++)
                t_xyz1[k][m] = xyztmp1[j][m] + bp_pos[m];
            strcpy(t_asym1[k], asym1[j]);
            t_btype1[k] = bpname[0];
            t_bnum1[k] = i;
        }

        ts2[i][1] = tn2 + 1;
        ts2[i][2] = num2;

        mul_dmatrix(xyztmp2, xyz2, bp_orien_T, num2, 3, 3, 3);
        for (j = 1; j <= num2; j++) {
            k = tn2 + j;
            for (m = 1; m <= 3; m++)
                t_xyz2[k][m] = xyztmp2[j][m] + bp_pos[m];
            strcpy(t_asym2[k], asym2[j]);
            t_btype2[k] = bpname[1];
            t_bnum2[k] = 2 * nbp - i + 1;
        }
        tn1 += num1;
        tn2 += num2;

        free_cmatrix(asym1, 1, ATM_NUM, 1, 4);
        free_cmatrix(asym2, 1, ATM_NUM, 1, 4);
        free_dmatrix(xyz1, 1, ATM_NUM, 1, 4);
        free_dmatrix(xyz2, 1, ATM_NUM, 1, 4);
        free_dmatrix(xyztmp1, 1, ATM_NUM, 1, 4);
        free_dmatrix(xyztmp2, 1, ATM_NUM, 1, 4);
    }

    tot_num = tn1 + tn2;

    tot_xyz = dmatrix(1, tot_num, 1, 3);
    tot_xyz2 = dmatrix(1, tot_num, 1, 3);
    tot_asym = cmatrix(1, tot_num, 1, 4);
    base_type = cvector(1, tot_num);
    strand = cvector(1, tot_num);
    base_num = lvector(1, tot_num);

    /* Combine two strands together! */
    for (i = 1; i <= tn1; i++) {
        strcpy(tot_asym[i], t_asym1[i]);
        base_type[i] = t_btype1[i];
        strand[i] = 'A';
        base_num[i] = t_bnum1[i];
        for (j = 1; j <= 3; j++)
            tot_xyz[i][j] = t_xyz1[i][j];
    }
    k = tn1;  /* total atom number in strand I */

    /* The order of strand II is now from 3'--->5', this need
       to be reversed according to information in <ts2> */
    for (i = nbp; i >= 1; i--) {
        for (j = ts2[i][1]; j <= ts2[i][1] + ts2[i][2] - 1; j++) {
            m = k + j - ts2[i][1] + 1;
            strcpy(tot_asym[m], t_asym2[j]);
            base_type[m] = t_btype2[j];
            strand[m] = 'B';
            base_num[m] = t_bnum2[j];
            for (i2 = 1; i2 <= 3; i2++)
                tot_xyz[m][i2] = t_xyz2[j][i2];
        }
        k += ts2[i][2];
    }

    /* Align the helix */
    align_helix(tot_xyz2, tot_xyz, tot_num, org_xyz, nbp);

    /* Write out the result */
    wrtpdb(tot_num, tot_asym, base_type, strand, base_num, tot_xyz2, filnam);

    /* free memory */
    free_cvector(dir_name, 1, BUF512);
    free_cvector(t_btype1, 1, MAX_ATOM);
    free_lvector(t_bnum1, 1, MAX_ATOM);
    free_cvector(t_btype2, 1, MAX_ATOM);
    free_lvector(t_bnum2, 1, MAX_ATOM);
    free_dvector(bp_pos, 1, 3);
    free_dvector(tmp_vec, 1, 3);
    free_dvector(hinge, 1, 3);
    free_cvector(base_type, 1, tot_num);
    free_cvector(strand, 1, tot_num);
    free_lvector(base_num, 1, tot_num);
    free_dmatrix(org_xyz, 1, nbp, 1, 3);
    free_cmatrix(t_asym1, 1, MAX_ATOM, 1, 4);
    free_dmatrix(t_xyz1, 1, MAX_ATOM, 1, 3);
    free_cmatrix(t_asym2, 1, MAX_ATOM, 1, 4);
    free_dmatrix(t_xyz2, 1, MAX_ATOM, 1, 3);
    free_lmatrix(ts2, 1, nbp, 1, 2);
    free_dmatrix(v_orien, 1, 3, 1, 3);
    free_dmatrix(v_orien_T, 1, 3, 1, 3);
    free_dmatrix(bp_orien, 1, 3, 1, 3);
    free_dmatrix(bp_orien_T, 1, 3, 1, 3);
    free_dmatrix(tot_xyz, 1, tot_num, 1, 3);
    free_dmatrix(tot_xyz2, 1, tot_num, 1, 3);
    free_cmatrix(tot_asym, 1, tot_num, 1, 4);
}

/* % BLOCK_GLH build DNA schematic structure using GLOBAL parameters
   %
   % Input: nblk--number of regular blocks
   %        glhpar--[nblk by 6] GLOBAL parameters
   %        xdir--+x-axis direction */
void block_GLH(long nblk, double **glhpar, long xdir, char *filnam)
{
    long i, j, k, m, vnum, nbond, **lkg, NBLK = 100;
    char blknam[BUF512], *dir_name, *env_add, **asym, **tot_asym;
    double twist_g = 0.0, rise_g = 0.0, tip, incl, y_disp, x_disp;

    double **blkxyz, **blkxyz2, **org_xyz;
    double **tot_xyz, **tot_xyz2;
    long tot_vnum, tot_nbond, **tot_lkg;
    long num_vtx = 0, num_lkg = 0;

    double **v_orien, **v_orien_T, *tmp_vec, *hinge, TipInc;
    double **bp_orien, **bp_orien_T, *bp_pos;

    /* Read in rigid block geometry */
    dir_name = cvector(1, BUF512);

    printf("\nGeometry data file name for the block (Dft BLOCK.alc): ");

    read_stdin_str(blknam);

    if (strlen(blknam) == 0)
        strcpy(blknam, "BLOCK.alc");

    if ((env_add = getenv("SCHNA_AR_P")) == NULL)
        nrerror("Environment variable SCHNA_AR_P undefined!!!\n");
    strcpy(dir_name, env_add);
    strcat(dir_name, "/BaseGeo/");
    strcat(dir_name, blknam);

    asym = cmatrix(1, NBLK, 1, 3);  /* 2 characters long */
    lkg = lmatrix(1, NBLK, 1, 2);
    blkxyz = dmatrix(1, NBLK, 1, 3);
    blkxyz2 = dmatrix(1, NBLK, 1, 3);

    rdalc(&vnum, &nbond, asym, lkg, blkxyz, dir_name);

    /* Set block orientation properly */
    if (xdir == 1)
        for (i = 1; i <= vnum; i++) {
            blkxyz[i][1] = -blkxyz[i][1];
            blkxyz[i][3] = -blkxyz[i][3];
        }

    tot_vnum = nblk * vnum;
    tot_nbond = nblk * nbond;

    org_xyz = dmatrix(1, nblk, 1, 3);
    tot_lkg = lmatrix(1, tot_nbond, 1, 2);
    tot_asym = cmatrix(1, tot_vnum, 1, 3);
    tot_xyz = dmatrix(1, tot_vnum, 1, 3);
    tot_xyz2 = dmatrix(1, tot_vnum, 1, 3);

    v_orien = dmatrix(1, 3, 1, 3);
    v_orien_T = dmatrix(1, 3, 1, 3);
    tmp_vec = dvector(1, 3);
    hinge = dvector(1, 3);

    bp_orien = dmatrix(1, 3, 1, 3);
    bp_orien_T = dmatrix(1, 3, 1, 3);
    bp_pos = dvector(1, 3);

    /* Get the whole data by iterations */
    for (i = 1; i <= nblk; i++) {
        x_disp = glhpar[i][1];
        y_disp = glhpar[i][2];
        rise_g += glhpar[i][3];
        incl = glhpar[i][4];
        tip = glhpar[i][5];
        twist_g += glhpar[i][6];

        rotz(v_orien, twist_g);  /* bp orientation after "twisting" */
        tr_dmatrix(v_orien_T, v_orien, 3, 3);
        tmp_vec[1] = x_disp;
        tmp_vec[2] = y_disp;
        tmp_vec[3] = rise_g;
        dvec_x_dmtx(bp_pos, tmp_vec, v_orien_T, 3, 3, 3);

        /* Each bp's origin in global reference frame */
        for (j = 1; j <= 3; j++)
            org_xyz[i][j] = bp_pos[j];

        for (j = 1; j <= 3; j++)
            hinge[j] = incl * v_orien[j][1] + tip * v_orien[j][2];  /* axis */
        TipInc = dist_ab(incl, tip);  /* angle */
        arbrot(v_orien_T, hinge, TipInc);  /* v_orien_T as rotmat */
        /* Final bp orientation */
        mul_dmatrix(bp_orien, v_orien_T, v_orien, 3, 3, 3, 3);
        tr_dmatrix(bp_orien_T, bp_orien, 3, 3);

        mul_dmatrix(blkxyz2, blkxyz, bp_orien_T, vnum, 3, 3, 3);
        for (j = 1; j <= vnum; j++) {
            k = num_vtx + j;
            strcpy(tot_asym[k], asym[j]);
            for (m = 1; m <= 3; m++)
                tot_xyz[k][m] = blkxyz2[j][m] + bp_pos[m];
        }

        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tot_lkg[num_lkg + j][k] = lkg[j][k] + (i - 1) * vnum;

        num_vtx += vnum;
        num_lkg += nbond;
    }

    /* Align the helix */
    align_helix(tot_xyz2, tot_xyz, tot_vnum, org_xyz, nblk);

    /* Write out the result */
    wrtalc(tot_vnum, tot_nbond, tot_asym, tot_lkg, tot_xyz2, filnam);

    free_cvector(dir_name, 1, BUF512);
    free_dvector(tmp_vec, 1, 3);
    free_dvector(hinge, 1, 3);
    free_dvector(bp_pos, 1, 3);
    free_cmatrix(asym, 1, NBLK, 1, 3);
    free_lmatrix(lkg, 1, NBLK, 1, 2);
    free_dmatrix(blkxyz, 1, NBLK, 1, 3);
    free_dmatrix(blkxyz2, 1, NBLK, 1, 3);
    free_dmatrix(org_xyz, 1, nblk, 1, 3);
    free_lmatrix(tot_lkg, 1, tot_nbond, 1, 2);
    free_cmatrix(tot_asym, 1, tot_vnum, 1, 3);
    free_dmatrix(tot_xyz, 1, tot_vnum, 1, 3);
    free_dmatrix(tot_xyz2, 1, tot_vnum, 1, 3);
    free_dmatrix(v_orien, 1, 3, 1, 3);
    free_dmatrix(v_orien_T, 1, 3, 1, 3);
    free_dmatrix(bp_orien, 1, 3, 1, 3);
    free_dmatrix(bp_orien_T, 1, 3, 1, 3);
}

void GLH_build(void)
     /* GLH_BUILD is a double helical structure building program using
        GLOBAL parameters. Program written in C by Xiang-Jun Lu (1996--1997) */
{
    char filnam[BUF512], filcnt[BUF512], str[BUF512], **bpseq, *ptr;
    long ich, imdl, xdir, nbp, ipos;
    double *ave_par, **bppar, **glhpar;

    system("clear");
    printf("Structure rebuilding based on GLOBAL parameters\n\n");
    printf("1. Uniform regular helix\n");
    printf("2. Use a full set of GLOBAL parameters\n");
    printf("3. Exact structure based on CEHS base-pair & GLOBAL parameters\n");
    printf("4. Quit\n");
    printf("Your choice(1-4, Dft 3): ");
    read_stdin_str(str);
    ich = atoi(str);
    if ((ich < 1) || (ich > 4))
        ich = 3;

    bpseq = cmatrix(1, 2, 1, MAX_BP);  /* NOTE: bpseq is [2 by MAX_BP] */
    bppar = dmatrix(1, MAX_BP, 1, 6);
    glhpar = dmatrix(1, MAX_BP, 1, 6);

    if (ich == 1)
        std_GLH(&nbp, bpseq, glhpar);
    else if (ich == 2)
        step_GLH(&nbp, bpseq, glhpar);
    else if (ich == 3)
        exact_GLH(&nbp, bpseq, bppar, glhpar);
    else
        exit(0);

    ave_par = dvector(1, 6);
    /* Twist & rise are ZEROs for the 1st bp */
    mean_dmatrix(ave_par, glhpar, nbp, 6);
    xdir = get_xdir(ave_par[6]);  /* Set +x-axis direction */

    printf("\nOutput data file name (Dft GLH_gen.out): ");
    read_stdin_str(filnam);
    if (strlen(filnam) == 0)
        strcpy(filnam, "GLH_gen.out");

    printf("\nChoose from one of the following representations:\n");
    printf("1. Real atomic model (in PDB format)\n");
    printf("2. Rectangular block model (in ALCHEMY format)\n");
    printf("Your choice (1 or 2, Dft 1): ");
    read_stdin_str(str);
    imdl = atoi(str);
    if ((imdl < 1) || (imdl > 2))
        imdl = 1;

    if (imdl == 1) {  /* Use atomic model */
        if (ich != 3)
            set_bppar(bppar, nbp, bpseq);  /* Get base-pair parameters */
        atomic_GLH(nbp, bpseq, bppar, glhpar, xdir, filnam);
    } else {  /* Use rectangular block model */
        block_GLH(nbp, glhpar, xdir, filnam);
        printf("\n\aFile %s has NO backbone connection information", filnam);

        ptr = strchr(filnam, '.');
        if (ptr != NULL) {
            ipos = ptr - filnam + 1;
            strncpy(filcnt, filnam, (size_t) ipos);
            filcnt[ipos] = '\0';  /* Just to make sure .... */
        } else {
            strcpy(filcnt, filnam);
            strcat(filcnt, ".\0");
        }
        strcat(filcnt, "cnt");
        alc_connect(filnam, filcnt);
        printf("\n\aFile %s HAS backbone connection information\n", filcnt);

        strcpy(filnam, filcnt);
    }
    str_display(imdl, filnam);

    free_dvector(ave_par, 1, 6);
    free_cmatrix(bpseq, 1, 2, 1, MAX_BP);
    free_dmatrix(bppar, 1, MAX_BP, 1, 6);
    free_dmatrix(glhpar, 1, MAX_BP, 1, 6);
}

/* %%%%%%%%%%%%%%%%%%% GLH Building functions END %%%%%%%%%%%%%%%%%%% */
