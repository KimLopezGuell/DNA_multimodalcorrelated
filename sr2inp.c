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

/*
   % To extract base paired residues and put them in SCHNAaP input format
   %
   % The input file should in the following format:
   % Line 1: Name of original PDB file
   % Line 2: Name of the new PDB file with only base-paired residues
   % Line 3: Name of SCHNAaP input file (with only base-paired residues)
   % Line 4: Number of base-pairs
   % From Line 5 up to Number-of-base-pairs Lines:
   % sI-name sI-residue-number  sII-name sII-residue-number
   %
   % where  sI: the leading strand I (1-char),
   %       sII: the complementary strand II (1-char).
   %      The four fields are separated by white space
   %
   % An example for PDT015 is as follows
   %     (See file <pdt015.sr> in directory "Examples" of the package)
   % ====================================================================
   % pdt015.pdb
   % pdt015_bp.pdb
   % pdt015_bp.inp
   % 10
   % C 2  D 11
   % C 3  D 10
   % C 4  D  9
   % C 5  D  8
   % C 6  D  7
   % C 7  D  6
   % C 8  D  5
   % C 9  D  4
   % C 10 D  3
   % C 11 D  2
   % ====================================================================
   %
   % NB: It is always a good idea to check your PDB file with RASMOL
   % --------------------------------------------------------------------
   %
   % Program written by Xiang-Jun Lu <xiangjun@rutchem.rutgers.edu>
   % Date: Dec. 6, 1997
   %
   % Please feel free to modify this utility program as necessary.
 */

#include "schna_ar_p.h"

int main(long argc, char *argv[])
{
    char filnam[80], orgpdb[80], outpdb[80], inpana[80];
    char **asym, **asymA, **asymB;
    char *btype, *btypeA, *btypeB;
    char *strand, *strandA, *strandB;
    char *sA, *sB;

    long i, j, k, itmp, num_bp;
    long num, numA = 0, numB = 0, numAB;
    long *bnum, *bnumA, *bnumB;
    long *bA, *bB, *idx, nidx;

    double **xyz, **xyzA, **xyzB;

    FILE *fid;

    /* Read in the base-pairing information */
    if (argc == 1) {
        printf("Name of the [strand-name, residue-number] file: ");
        read_stdin_str(filnam);
    } else
        strcpy(filnam, argv[1]);

    if ((fid = fopen(filnam, "r")) == NULL) {
        printf("Failed to open file: %s\n", filnam);
        exit(1);
    }

    fscanf(fid, "%s", orgpdb);
    fscanf(fid, "%s", outpdb);
    fscanf(fid, "%s", inpana);
    fscanf(fid, "%ld", &num_bp);

    sA = cvector(1, num_bp);
    sB = cvector(1, num_bp);
    bA = lvector(1, num_bp);
    bB = lvector(1, num_bp);

    for (i = 1; i <= num_bp; i++) {
        j = num_bp - i + 1;
        fscanf(fid, " %c %ld %c %ld", &sA[i], &bA[i], &sB[j], &bB[j]);
    }

    fclose(fid);

    /* Read in the original PDB data file */
    asym = cmatrix(1, MAX_ATOM, 1, 4);
    btype = cvector(1, MAX_ATOM);
    strand = cvector(1, MAX_ATOM);
    bnum = lvector(1, MAX_ATOM);
    xyz = dmatrix(1, MAX_ATOM, 1, 3);

    rdpdb(&num, asym, btype, strand, bnum, xyz, orgpdb);

    /* Process each residue and put the sorted information in I and II */
    asymA = cmatrix(1, num, 1, 4);
    btypeA = cvector(1, num);
    strandA = cvector(1, num);
    bnumA = lvector(1, num);
    xyzA = dmatrix(1, num, 1, 3);

    asymB = cmatrix(1, num, 1, 4);
    btypeB = cvector(1, num);
    strandB = cvector(1, num);
    bnumB = lvector(1, num);
    xyzB = dmatrix(1, num, 1, 3);

    idx = lvector(1, 100);  /* Assume maximum atom per residue <= 100 */

    for (i = 1; i <= num_bp; i++) {
        /* Base in leading strand I */
        nidx = 0;
        for (j = 1; j <= num; j++)
            if (strand[j] == sA[i] && bnum[j] == bA[i])
                idx[++nidx] = j;

        if (!nidx) {
            printf("No residue <%4ld> on strand <%c> in this PDB file!\n", bA[i], sA[i]);
            exit(1);
        }
        for (j = 1; j <= nidx; j++) {
            itmp = numA + j;
            strcpy(asymA[itmp], asym[idx[j]]);
            btypeA[itmp] = btype[idx[j]];
            strandA[itmp] = strand[idx[j]];
            bnumA[itmp] = bnum[idx[j]];
            for (k = 1; k <= 3; k++)
                xyzA[itmp][k] = xyz[idx[j]][k];
        }
        numA += nidx;

        /* Base in strand II */
        nidx = 0;
        for (j = 1; j <= num; j++)
            if (strand[j] == sB[i] && bnum[j] == bB[i])
                idx[++nidx] = j;
        if (!nidx) {
            printf("No residue <%4ld> on strand <%c> in this PDB file!\n", bB[i], sB[i]);
            exit(1);
        }
        for (j = 1; j <= nidx; j++) {
            itmp = numB + j;
            strcpy(asymB[itmp], asym[idx[j]]);
            btypeB[itmp] = btype[idx[j]];
            strandB[itmp] = strand[idx[j]];
            bnumB[itmp] = bnum[idx[j]];
            for (k = 1; k <= 3; k++)
                xyzB[itmp][k] = xyz[idx[j]][k];
        }
        numB += nidx;
    }

    /* Combine II with I, both along 5'--->3' direction */
    numAB = numA + numB;
    for (i = 1; i <= numB; i++) {
        j = numA + i;
        strcpy(asymA[j], asymB[i]);
        btypeA[j] = btypeB[i];
        strandA[j] = strandB[i];
        bnumA[j] = bnumB[i];
        for (k = 1; k <= 3; k++)
            xyzA[j][k] = xyzB[i][k];
    }

    /* Write out the base-paired PDB file */
    wrtpdb(numAB, asymA, btypeA, strandA, bnumA, xyzA, outpdb);

    /* Write out simplified data file as input to SCHNAaP */
    if ((fid = fopen(inpana, "w")) == NULL) {
        printf("Failed to open file: %s\n", filnam);
        exit(1);
    }
    for (i = 1; i <= numAB; i++)
        fprintf(fid, "%3s %c %3ld %8.3f %8.3f %8.3f\n", asymA[i], btypeA[i],
                bnumA[i], xyzA[i][1], xyzA[i][2], xyzA[i][3]);

    fclose(fid);

    free_cvector(sA, 1, num_bp);
    free_cvector(sB, 1, num_bp);
    free_lvector(bA, 1, num_bp);
    free_lvector(bB, 1, num_bp);

    free_cmatrix(asym, 1, MAX_ATOM, 1, 4);
    free_cvector(btype, 1, MAX_ATOM);
    free_cvector(strand, 1, MAX_ATOM);
    free_lvector(bnum, 1, MAX_ATOM);
    free_dmatrix(xyz, 1, MAX_ATOM, 1, 3);

    free_cmatrix(asymA, 1, num, 1, 4);
    free_cvector(btypeA, 1, num);
    free_cvector(strandA, 1, num);
    free_lvector(bnumA, 1, num);
    free_dmatrix(xyzA, 1, num, 1, 3);

    free_cmatrix(asymB, 1, num, 1, 4);
    free_cvector(btypeB, 1, num);
    free_cvector(strandB, 1, num);
    free_lvector(bnumB, 1, num);
    free_dmatrix(xyzB, 1, num, 1, 3);

    free_lvector(idx, 1, 100);

    return 0;
}
