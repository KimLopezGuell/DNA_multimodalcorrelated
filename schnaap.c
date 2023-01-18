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
    long ipos, i, j, numP1;

    long *RY_idx, *base_num, num_atom, num_bp, nmc, n_base, xdir;
    char **atom_id, **bp_seq, *base_type, *strand;
    double **xyz, mtwist;

    /* Get atomic listing */
    long *C1, *N, *P, *O4, *C68, *chi, *mc, *sugar, *base;

    /* Reference triads */
    double **bx, **by, **bz, **borg;

    /* MST triads */
    double **mxyz, **morg, **bp_par, **P_mst;

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
    strand = cvector(1, MAX_ATOM);
    base_num = lvector(1, MAX_ATOM);

    /* Read in the data file and get the information about each atom */
    get_data(RY_idx, atom_id, bp_seq, xyz, base_type, strand, base_num,
             &num_atom, &num_bp, inpfil, fo);

    /* Get atomic listing */
    C1 = lvector(1, 2 * num_bp);
    N = lvector(1, 2 * num_bp);
    P = lvector(1, 2 * num_bp);
    O4 = lvector(1, 2 * num_bp);
    C68 = lvector(1, 2 * num_bp);
    chi = lvector(1, 8 * num_bp);
    mc = lvector(1, 12 * num_bp);
    sugar = lvector(1, 10 * num_bp);
    base = lvector(1, num_atom);

    get_list(C1, N, P, O4, C68, chi, mc, &nmc, sugar, base, &n_base,
             RY_idx, atom_id, num_bp, num_atom, xyz);

    if (((nmc / 2) % 2) == 1)
        numP1 = num_bp - 1;  /* Normal case */
    else
        numP1 = num_bp;  /* P is the first atom */

    /* Define the global reference frame */
    set_str(xyz, C1, N, C68, base, base_num, n_base, num_atom, num_bp, fo);

    /* Get reference triad of each base and base-pair */
    bx = dmatrix(1, num_bp, 1, 9);
    by = dmatrix(1, num_bp, 1, 9);
    bz = dmatrix(1, num_bp, 1, 9);
    borg = dmatrix(1, num_bp, 1, 9);
    get_frame(bx, by, bz, borg, &xdir, RY_idx, C68, num_bp, base, n_base,
              base_num, atom_id, xyz, bp_seq, fo);

    /* Get local CEHS parameters, and MST */
    mxyz = dmatrix(1, num_bp - 1, 1, 9);
    morg = dmatrix(1, num_bp - 1, 1, 3);
    bp_par = dmatrix(1, num_bp, 1, 6);
    P_mst = dmatrix(1, num_bp - 1, 1, 6);
    CEHS_par(mxyz, morg, bp_par, P_mst, &mtwist,
             bx, by, bz, borg, P, numP1, num_bp, bp_seq, xyz, filstr, fo);

    /* Set each step with reference to its MST */
    wrt_mst(atom_id, base_type, strand, base_num, xyz, bp_seq,
            num_atom, num_bp, mxyz, morg, filstr);

    /* Deduce a set of global helical parameters */
    GLH_par(bx, by, bz, borg, num_bp, bp_seq, bp_par, filstr, fo);

    /* Reset the structure w.r.t the global reference frame */
    wrt_hel_dat(num_atom, atom_id, base_type, strand, base_num, xyz,
                base, n_base, filstr);

    /* Write out P_mst coordinates, and overall classifications */
    classify(P_mst, num_bp, bp_seq, xdir, mtwist, fo);

    /* Get backbone torsion angle and glycosyl angle chi */
    main_chain(mc, nmc, chi, num_bp, xyz, bp_seq, fo);

    /* Get the conformational parameters of the sugar ring. */
    sugar_ana(sugar, num_bp, xyz, bp_seq, fo);

    /* Get lambda: C1'--N1/C1'--N9 with respect to C1'--C1' */
    get_lambda(C1, N, C68, num_bp, xyz, bp_seq, fo);

    /* Get groove width based on P-P (-5.8A) and O4'-O4'(-2.8A) distance */
    groove_width(P, numP1, O4, xyz, num_bp, bp_seq, fo);

    /* Get the polar coordinates etc of P/C1'/O4' atoms */
    if (numP1 > 2) {
        prt_sep(fo, '*', 76);
        fprintf(fo, "Backbone polar coordinates with respect to the ");
        fprintf(fo, "global reference frame\n");
        fprintf(fo, "    Note: Parameters other than `dxyz' (the interatomic distance)");
        fprintf(fo,
                "\n          would be MEANINGLESS for strongly curved structures!\n\n");
        polar_coord(P, numP1, "Phosphorus", xyz, bp_seq, fo);
        prt_sep(fo, '-', 76);
        polar_coord(C1, num_bp, "C1'", xyz, bp_seq, fo);
        prt_sep(fo, '-', 76);
        polar_coord(O4, num_bp, "O4'", xyz, bp_seq, fo);
    }

    fclose(fo);

    free_lvector(RY_idx, 1, MAX_ATOM);
    free_cvector(base_type, 1, MAX_ATOM);
    free_cvector(strand, 1, MAX_ATOM);
    free_lvector(base_num, 1, MAX_ATOM);
    free_lvector(C1, 1, 2 * num_bp);
    free_lvector(N, 1, 2 * num_bp);
    free_lvector(P, 1, 2 * num_bp);
    free_lvector(O4, 1, 2 * num_bp);
    free_lvector(C68, 1, 2 * num_bp);
    free_lvector(chi, 1, 8 * num_bp);
    free_lvector(mc, 1, 12 * num_bp);
    free_lvector(sugar, 1, 10 * num_bp);
    free_lvector(base, 1, num_atom);
    free_cmatrix(atom_id, 1, MAX_ATOM, 1, 4);
    free_cmatrix(bp_seq, 1, MAX_ATOM, 1, 3);
    free_dmatrix(xyz, 1, MAX_ATOM, 1, 3);
    free_dmatrix(bx, 1, num_bp, 1, 9);
    free_dmatrix(by, 1, num_bp, 1, 9);
    free_dmatrix(bz, 1, num_bp, 1, 9);
    free_dmatrix(borg, 1, num_bp, 1, 9);
    free_dmatrix(mxyz, 1, num_bp - 1, 1, 9);
    free_dmatrix(morg, 1, num_bp - 1, 1, 3);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(P_mst, 1, num_bp - 1, 1, 6);

    return 0;
}
