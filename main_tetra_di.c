/* --------------- CHANGE IN THIS PROGRAM TO ORIGINAL ---------------- */

/* 

1) The arrays seq1,stif and geom are adjusted so that to each single bp-step a stiffness matrix can be assigned and not just an average stiffness matrix as before 

2) Another output file "table_all_(helpar).dat" is provided which gives all (helpar) parameters for each MC optimization (each line of the table equals the (helpar) parameters of the examined sequence)

*/
//Some changes have been made to the old Jurgen version to include correlation effects in both the unimodal approach (selection number 5) and multimodal approach (selection number 6).



/* Headers from DNA Flex */
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "dna_flex.h"	// Functions used in DNA Flex
#include "energy.h"	// Functions for calculating int energies

/* Headers from SCHNArP */

#include "schna_ar_p.h"

#include <time.h>

/* End Headers */




/*--------------------------Beginning of DNAFlex-------------------------*/


int main(int argc, char *argv[])
{

char seq_file[200];
sscanf(argv[1], "%s", seq_file);	// get sequence

int buf_num_it;
buf_num_it = atoi (argv[2]);	// get number of structures to generate

char out_f[1000],out_folder[1000];
sscanf(argv[3], "%s", out_f);        // get output folder
sprintf(out_folder, "%s", out_f);
printf("%s\n",out_folder);

int buf_selection;
buf_selection = atoi (argv[4]);	// selection of form of code (0=individual,1=dimer,2=tetramer,3=tet/dimerend,4=tet/dimerend multimodal,5=tet/dimerend/correlations,6=tet/dimerend/multimodal/correlations)

/* Filenames for sequence and stiffness (only applicable for individual bp steps) */

char stif_file[] = "input_stiff_avg_1kx5nucl_with4qlc_ends1.dat";
char avg_file[] = "stif_bsc1_k_avg_miniabc.dat"; // for average values
char avg_file_mumo[] = "stif_bsc1_k_avg_miniabc_allclust_TAAA_TTAA_corrected.dat";
char avg_file_cor[] = "stif_bsc1_k_avg_miniabc_cor.dat"; 
char avg_file_cor_mumo[] = "stif_bsc1_k_avg_miniabc_cor_multimodal.dat"; //for miniabc
//char avg_file_cor_mumo[] = "stif_bsc1_k_avg_miniabc_cor_40mer.dat"; //when I want to do the 40mer

/* Initialize all variables at beginning of code!!*/ 

int nmc,i0,ibuff2,juega,last;
int iflag; 	//it should be a bool, but I dont know how to assign true and false to a bool in c
double e_dist, e_dist2;    
char **seq1, **ist; //seq1 saves the read sequence, ist all the possible tetramers
double ***stif, ***stiff, **geom; // stif contains the read data matrices, stiff the assigned ones and geom the read data equilibrium values
int ord[6];   // no one knows the reason of its existence
double **xconfi, **xconf0, **xconf; // xconfi saves the values of the initial helpars, xconf0 the equilibrium ones and xconf the current accepted ones
double **xconft, *ener, *enerH; // xconft saves the currently tested new values of the helpars. ener is the energy of the accepted chain, enerH its harmonically computed analogous one
double sc[6],dd[6],dd2[6],dd3[6],dd4[6],dd5[6]; //sc is the scaling matrix for the change in helpars. dd holds the difference of helpars between the tested values and the equilibrium values. We only use dd2,dd3,dd4 and dd5 for computation of energy taking into account correlation between i,i+1 and i,i-1
double **dx; //dx[6][1000];
char *names[6];   				// No maximum size of char in C
char line;
double boltz=0.00198717;  //boltzmann constant in kcal/mol/K
double enex,enexA,enexB; //values of the energy of the currently tested bp, bp+1 (After), bp-1 (Before)
int itot; //number of bp of the sequence. It starts at 0. The last one is itot-1.
int i,j,k,l;	// Iteration parameters
int it_mc;	// Test parameter for # MC moves
char **output;  // output[1000][100];			// Writing helical coordinates in this string
int it_prog,phi;  //Total number of structures to simulate and its iteration parameter
double acc_multicor; //For counting the acceptance of structures in select_state in the multimodal correlated model (Ising version)
double ind_energy;
double *ind_ener;
float sum_energy;
int selection;
char filename[1100];
double ***stif_tetra, ***stif_tetraS, ***stif_di, ***stif_tetra_mumo, ***prob_tetra_mumocor, ***prob_tetra_mumocor_all, **geom_tetra, **geom_di, **geom_tetra_mumo, *perc, **perc_seq; //stif_tetra saves the covariance matrices for tetramer (12x12 if correlated), stif_tetraS the 6x6 extra ones for correlation simulation, stif_di for dimer, stif_tetra_mumo for mumo, prob_tetra_mumocor saves the cluster probabilities for Ising. prob_tetra_mumocor_all reads the data, as well as geom_tetra, geom_di, geom_tetra_mumo. perc and perc_seq save the probabilities of the multimodal states.
//Two stiffness matrices arrays are needed for correlation because the banded matrix approach requires us to sum different blocks. The 12x12 are in stif_tetra, already computed. The 6x6 ones in stif_tetraS are for the bp itot-4, which does not have correlation with i+1.
int *num_mumo, *sel_state, *read_mumo; //num_mumo has the number of states of all the bp pf the sequence, sel_state the states selected. 
char **seq1_tetra, **seq2_tetra, **seq3_tetra, **seq4_tetra, **seq1_di, **seq1_tetra_mumo;  //sequence read and assigned (dimer or tetramer, some for stiffnes and some for probabilities)
int it_snapshots, it_mcsteps, alpha, it_mod;
char **out_part;
double endenergy;

double **xyz_float;
int msec,msec1,msec2,msec3;
double **center_bp;
double r_bp;
char **pdb_data;
char **tmp_at;
char *baset;
char *strand;
int *base_num;
double **tmp_xyz;
double *q_fl;
int num_bp;
double **d_bp;
double **d_bp_fl;
double **lj_bp, **lj_bp_fl;
double **mid_bp;
double **dir_bp;
double **dir_ph;
double **pos_ph;
double **dh_ph;
double **dh_ph_fl;
double q_ph;
double **d_ph1, **d_ph2;
double **d_ph_fl;
double tot_dh_lj, tot_lj, tot_dh;
int k1,k2,k3,k4;
char **pdb_tot;
int tot_num_at;
int num_fl;	// Number of floating patricles
double shift;
double shift_in;
int acc;
int acc_ph;
double bp_excl_vol;
double bp_fl_excl_vol;
double lj_k,lj_cut_up,lj_cut_low,lj_a_bp,lj_a_fl;
double eps,kappa,f_pi_eps0,dh_cut_low,dh_cut_up;
double **xconf_buf;
double tot_en,prev_tot_en;
double d_ph_midbp;
double **xyz_float_buf;
int count_no_excl;
int *base_num_ph;
char *baset_ph;
char *strand_ph;
char **n7_data;
int *base_num_n7;
double **n7_xyz;
double **n7_xyz_buf;
double **c6_xyz, **c8_xyz;
double **vec_c6_c8, **vec_n7_c8, *cos_alpha, **vec_z, **vec_n7_z, **x_axis_bp, **y_axis_bp;
double **mid_bp_x, **mid_bp_y;
double **ph_pos_ex;
char filnam[BUF512];
double **pos_store, ***orient_store, ***rot_store;
char **table_bp;

int len_seq;	// Variable to determine dynamic size of all arrays
int seq_pcs;	// Length of sequence pieces to save (i.e. 4 for tetramer seq)
int no_helpars;
int len_out_line;
int len_mumo;
len_mumo = 7000;
seq_pcs = 10;
no_helpars = 6;
len_out_line = 100;
num_fl = 0;			// Number of floating particles
int no_rot_pars = 3;

get_seq_len(&len_seq,seq_file);

int alc_space;
if(len_seq < 256){alc_space = 256;}else{alc_space = len_seq;}

/* Allocate memory for arrays */
seq1 = char_aloc2d(alc_space,seq_pcs);
ist  = char_aloc2d(len_seq,seq_pcs);
stif = dbl_aloc3d(alc_space,7*no_helpars,7*no_helpars);
stiff = dbl_aloc3d(len_seq,7*2*no_helpars,7*2*no_helpars);
geom = dbl_aloc2d(alc_space,no_helpars);
xconfi = dbl_aloc2d(no_helpars,len_seq);
xconf0 = dbl_aloc2d(7*no_helpars,len_seq);
xconf = dbl_aloc2d(no_helpars,len_seq);
xconft = dbl_aloc2d(no_helpars,len_seq);
ener = malloc(len_seq*sizeof(double));
enerH = malloc(len_seq*sizeof(double));
dx = dbl_aloc2d(no_helpars,len_seq);
output  = char_aloc2d(len_seq,len_out_line);

seq1_tetra = char_aloc2d(65536,seq_pcs);
seq2_tetra = char_aloc2d(65536,seq_pcs);
seq3_tetra = char_aloc2d(65536,seq_pcs);
seq4_tetra = char_aloc2d(65536,seq_pcs);
seq1_di = char_aloc2d(16,seq_pcs);
seq1_tetra_mumo = char_aloc2d(len_mumo,seq_pcs);
stif_tetra = dbl_aloc3d(65536,2*no_helpars,2*no_helpars); 
stif_tetraS = dbl_aloc3d(65536,no_helpars,no_helpars);
stif_di = dbl_aloc3d(16,no_helpars,no_helpars);
stif_tetra_mumo = dbl_aloc3d(len_mumo,no_helpars,no_helpars);
prob_tetra_mumocor = dbl_aloc3d(7*len_seq,7,7);
prob_tetra_mumocor_all = dbl_aloc3d(65536,7,7);  
geom_tetra = dbl_aloc2d(65536,no_helpars);
geom_di = dbl_aloc2d(16,no_helpars);
geom_tetra_mumo = dbl_aloc2d(len_mumo,no_helpars);
perc = malloc(len_mumo*sizeof(double));
perc_seq = dbl_aloc2d(len_seq,2*no_helpars);
num_mumo = malloc(len_mumo*sizeof(int));
read_mumo = malloc(len_mumo*sizeof(int));
sel_state = malloc(len_seq*sizeof(int));

out_part  = char_aloc2d(len_seq,len_out_line);
xyz_float = dbl_aloc2d(num_fl,no_helpars);
q_fl = malloc(len_seq*sizeof(double));
center_bp = dbl_aloc2d(len_seq,no_helpars);
pdb_data = char_aloc2d(2*len_seq,len_out_line);
tmp_at = char_aloc2d(2*len_seq,no_helpars);
baset = malloc(2*len_seq*sizeof(char));
strand = malloc(2*len_seq*sizeof(char));
base_num = malloc(2*len_seq*sizeof(int));
tmp_xyz = dbl_aloc2d(2*len_seq,no_helpars);
d_bp = dbl_aloc2d(len_seq,len_seq);
d_bp_fl = dbl_aloc2d(len_seq,len_seq);
lj_bp = dbl_aloc2d(2*len_seq,no_helpars);
lj_bp_fl = dbl_aloc2d(2*len_seq,no_helpars);
mid_bp = dbl_aloc2d(len_seq,no_helpars);
dir_bp = dbl_aloc2d(len_seq,no_helpars);
dir_ph = dbl_aloc2d(len_seq,no_helpars);
pos_ph = dbl_aloc2d(2*len_seq,no_helpars);
d_ph1 = dbl_aloc2d(2*len_seq,2*len_seq);
d_ph2 = dbl_aloc2d(2*len_seq,2*len_seq);
d_ph_fl = dbl_aloc2d(2*len_seq,len_seq);
dh_ph = dbl_aloc2d(4*len_seq,no_helpars);
dh_ph_fl = dbl_aloc2d(2*len_seq,no_helpars);

xconf_buf = dbl_aloc2d(no_helpars,len_seq);
xyz_float_buf = dbl_aloc2d(num_fl,no_helpars);

pdb_tot = char_aloc2d(20*len_seq,len_out_line);
base_num_ph = malloc(2*len_seq*sizeof(int));
baset_ph = malloc(2*len_seq*sizeof(char));
strand_ph = malloc(2*len_seq*sizeof(char));
n7_data = char_aloc2d(2*len_seq,len_out_line);
base_num_n7 = malloc(2*len_seq*sizeof(int));
n7_xyz = dbl_aloc2d(2*len_seq,no_helpars);
n7_xyz_buf = dbl_aloc2d(2*len_seq,no_helpars);
c6_xyz = dbl_aloc2d(len_seq,no_helpars);
c8_xyz = dbl_aloc2d(len_seq,no_helpars);
vec_c6_c8 = dbl_aloc2d(len_seq,no_helpars);
vec_n7_c8 = dbl_aloc2d(len_seq,no_helpars);
cos_alpha = malloc(len_seq*sizeof(double));
vec_z = dbl_aloc2d(len_seq,no_helpars);
vec_n7_z = dbl_aloc2d(len_seq,no_helpars);
x_axis_bp = dbl_aloc2d(len_seq,no_helpars);
y_axis_bp = dbl_aloc2d(len_seq,no_helpars);
mid_bp_x = dbl_aloc2d(len_seq,no_helpars);
mid_bp_y = dbl_aloc2d(len_seq,no_helpars);
ph_pos_ex = dbl_aloc2d(2*len_seq,no_helpars);

orient_store = dbl_aloc3d(len_seq,no_rot_pars,no_rot_pars);
pos_store = dbl_aloc2d(len_seq,no_rot_pars);
rot_store = dbl_aloc3d(2*len_seq,no_rot_pars,no_rot_pars);
table_bp  = char_aloc2d(len_seq,len_out_line);

/* Initialize arrays and files for huge tables */

char **table_heli_shif;
char **table_heli_slid;
char **table_heli_rise;
char **table_heli_tilt;
char **table_heli_roll;
char **table_heli_twis;
char **table_state; // selected states for all pbs in multimodal approach
char **table_e_dist; // end 2 end distance 
char **table_e_dist2; // end 2 end distance for bo 2 to 13 (checking and analysis purposes, not real end 2 end distance)

table_heli_shif  = char_aloc2d(len_seq,len_out_line);
table_heli_slid  = char_aloc2d(len_seq,len_out_line);
table_heli_rise  = char_aloc2d(len_seq,len_out_line);
table_heli_tilt  = char_aloc2d(len_seq,len_out_line);
table_heli_roll  = char_aloc2d(len_seq,len_out_line);
table_heli_twis  = char_aloc2d(len_seq,len_out_line);
table_state  = char_aloc2d(len_seq,len_out_line);
table_e_dist= char_aloc2d(len_seq,len_out_line); 
table_e_dist2= char_aloc2d(len_seq,len_out_line);

sprintf(filename,"%s/output_tables_helpar/table_all_shif.dat",out_folder);
FILE *file_shif;
file_shif=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_slid.dat",out_folder);
FILE *file_slid;
file_slid=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_rise.dat",out_folder);
FILE *file_rise;
file_rise=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_tilt.dat",out_folder);
FILE *file_tilt;
file_tilt=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_roll.dat",out_folder);
FILE *file_roll;
file_roll=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_all_twis.dat",out_folder);
FILE *file_twis;
file_twis=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_states.dat",out_folder);
FILE *file_states;
file_states=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_e_dist.dat",out_folder);
FILE *file_dist;
file_dist=fopen(filename, "w");

sprintf(filename,"%s/output_tables_helpar/table_e_dist2.dat",out_folder);
FILE *file_dist2;
file_dist2=fopen(filename, "w");


/* ------- End Initialze for tables and files --------- */

srand (time(NULL));		// initialize random seed (the random seed is initialized to a value representing the current time (calling time) to generate a different value every time the program is run.)

clock_t start_all = clock(), diff_all; //Timer

/* icuan = 1;
icnt = 1; */

/* This is the part which should be read from a txt file (see fortran code) PROBLEM!!! */

nmc=100000;
double agr=1.0;	//agr=1.0 usually
int ibuff=10000;
ibuff2=100;
double temp=298.0,fact=0.3; //fact is a factor used in the MC algorithm to tune how strong we let the movements be
last=1000000;
sc[0]=5.148;
sc[1]=3.558;
sc[2]=2.052;
sc[3]=29.472;
sc[4]=34.266;
sc[5]=53.862;
names[0]="shif";
names[1]="slid";
names[2]="rise";
names[3]="tilt";
names[4]="roll";
names[5]="twis";

/* ----------------------------------End of part--------------------------------------- */



/* ------------------------------- Looping and selection variables ------------------------------- */


selection = buf_selection;			// 0 = individual, 1 = dimer, 2 = tetramer, 3 = tetramer with dimer end; 4 = tetramer with dimer end (multimodal); 5 = tetramer with dimer end (unimodal, with correlations between bp i and i+1, i and i-1); 6 = tetramer with dimer end (multimodal, with correlations whenever possible (bp 3 to itot-4), computed in the MC algorithm with 12x12 matrices and in selection process before MC with an Ising approach)


it_mcsteps = 1;		// number of MC steps in floating process before looking at position of particle
it_snapshots = 0;		// number of pdb's generated during the floating process
it_prog = buf_num_it;		// number of executing all the MC and reconstruction part (THE WHOLE PROGRAM)

it_mc = 5000*len_seq;		// number of MC steps to bring structure to minimum energy before particle inserted. Number tested by Jürgen.
//it_mc = 2;				
				/* 100000 -> 56mer
				   10000000 -> 1000mer
				   2000000 -> 450mer
				   5000000 -> 750mer
				*/



/* ------------------------------- END Looping and selection variables ------------------------------- */


/* ------------------------------- Parameters for dynamical part ------------------------------------- */

	
shift_in = 600.0;		// Volume occupied by floating particles
shift = 5.0;			// Range of Shift of floating particle per MC move

/* ---- Parameters for excluded volume ---- */
acc = 10;			// bp-steps for calculating distance
acc_ph = 20;			// phosphate = bp steps for calculating distance
bp_excl_vol = 20;		// minimum distance between two bp centers
bp_fl_excl_vol = 10;		// minimum distance between bp center and floating particle

d_ph_midbp = 10.25;		// Distance of Phosphate of helical axis



/* ---- Parameters for LJ ---- */	
lj_k = 4*0.0257; 		// in eV (0.0257 eV = 1kT, adapted from Arya 2014)
lj_cut_up = 1.5*bp_excl_vol;	
lj_cut_low = 0.1;

lj_a_bp = bp_excl_vol;		// = 10A as excluded volume (radius of excluded volume around center_bp pos)
lj_a_fl = bp_fl_excl_vol;	// 



/* ---- Parameters for DH ---- */
eps = 80.0;
kappa = 0.033; 			// unit: 1/A	taken from Arya paper
f_pi_eps0 = 1.1126*pow(10,-10);
q_ph = -1.0;

dh_cut_low = 0.1;
dh_cut_up = bp_excl_vol;	// 20.0


/* ------------------------------- END Parameters for dynamical part ----------------------------------- */


order(names,ord);


if(selection == 0){		// Select individual bsc1 stiffness matrices
readsq_ind(&itot,ist,seq_file);

readss_ind(stif,seq1,geom,ord,seq_file,stif_file);

assign_ind(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 1){		// Select average bsc0 dimer stiffness matrices
readsq_di(&itot,ist,seq_file);

readss_di(stif,seq1,geom,ord);

assign_di(itot,ist,seq1,stif,stiff,geom,xconf0);
}

if(selection == 2){		// Select average bsc0 tetramer stiffness matrices
readsq_tetra(&itot,ist,seq_file);

readss_tetra(stif,seq1,geom,ord,avg_file);

assign_tetra(itot,ist,seq1,stif,stiff,geom,xconf0);
}


if(selection == 3){		// Select average bsc1 tetramer stiffness matrices with bsc1 dimer at end
readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di(stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord,avg_file);

assign_tetra_di(itot,ist,stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0);

}

if(selection == 4){		// Select multimodal average bsc1 tetramer stiffness matrices with bsc1 dimer at end 
readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di(stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord,avg_file);

readss_tetra_mumo(stif_tetra_mumo,seq1_tetra_mumo,geom_tetra_mumo,perc,ord,avg_file_mumo,len_mumo); 

assign_tetra_di(itot,ist,stif_tetra,seq1_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0);

assign_tetra_di_mumo(itot,ist,stif_tetra_mumo,seq1_tetra_mumo,geom_tetra_mumo,stif_di,seq1_di,geom_di,stiff,xconf0,num_mumo,len_mumo,perc,perc_seq);

}

if(selection == 5){		// Select average bsc1 tetramer stiffness matrices with bsc1 dimer at end that contemplate correlations between bp at distance 1

readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di_cor(stif_tetra,stif_tetraS,seq1_tetra,seq2_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord,avg_file_cor);

assign_tetra_di_cor(itot,ist,stif_tetra,stif_tetraS,seq1_tetra,seq2_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0);

}

if(selection == 6){		// Select average bsc1 tetramer stiffness matrices with bsc1 dimer end but multimodal average xconf0 and probabilities for correlations between bp at distance 1

readsq_tetra_di(&itot,ist,seq_file);

readss_tetra_di_cor_mumo(stif_tetra,stif_tetraS,seq1_tetra,seq2_tetra,geom_tetra,stif_di,seq1_di,geom_di,ord,avg_file_cor_mumo,read_mumo);
//assign cor 12x12 matrices to central bps (except tails)
assign_tetra_di_cor_mumo(itot,ist,stif_tetra,stif_tetraS,seq1_tetra,seq2_tetra,geom_tetra,stif_di,seq1_di,geom_di,stiff,xconf0,read_mumo);

readss_tetra_mumo(stif_tetra_mumo,seq1_tetra_mumo,geom_tetra_mumo,perc,ord,avg_file_mumo,len_mumo); 
//assign 6x6 matrices to tails and equilibrium params to the whole sequence
assign_tetra_di_mumo_notK(itot,ist,stif_tetra_mumo,seq1_tetra_mumo,geom_tetra_mumo,stif_di,seq1_di,geom_di,stiff,xconf0,num_mumo,len_mumo,perc,perc_seq);

readss_mumocor_di(itot, prob_tetra_mumocor_all, read_mumo, seq3_tetra, seq4_tetra); //for miniabc
assign_mumocor_di(itot, prob_tetra_mumocor, prob_tetra_mumocor_all, read_mumo, num_mumo, seq3_tetra, seq4_tetra, ist); //for miniabc
//readss_mumocor_40mer(itot, prob_tetra_mumocor,num_mumo);  //Only if simulating the 40mer

}	

readini(itot,xconfi); //Reading the initial configuration and saving it in xconfi. Hardcoded file name in the function, in dna_flex.h.


/* ------------------------------- Looping starts ------------------------------- */

ind_ener = malloc(it_snapshots*sizeof(double));

for(phi=0;phi<it_prog;phi++){			// do optimization it_prog times to check if MC converges

	clock_t start = clock(), diff; //Timer

	//printf("0.5\n");
	// Compute initial energy. Beware that we assign a higher energy to each bp (the accumulated sum rather than the actual individual bp computed one) to help with optimization afterwards.
	if(selection == 5) {
	energy_cor(itot,xconfi,xconf0,dx,stiff,ener);
	}
	else if (selection == 6) {
	energy_cor(itot,xconfi,xconf0,dx,stiff,ener);
	energy(itot,xconfi,xconf0,dx,stiff,enerH);
	}
	else {
	energy(itot,xconfi,xconf0,dx,stiff,ener);
        }

	double ener0=0.0;
	double enertot;
	double enertotH; //enertotH is used in the multimodal correlated model because of the harmonic filter needed to prevent extreme energies (due to matrix divergence outside of clustered hp range)

	for(i=0; i<itot; i++){ener0 = ener0 + ener[i];}

	enertot = ener0;

	ener0 = 0.0;
	for(i=0; i<itot; i++){ener0 = ener0 + enerH[i];}

	enertotH = ener0;

	//printf("0.6\n");

	// Set working coordinates and energies equal to the starting ones

	transfer_hel_coord(itot,xconf,xconfi);

	/*for(i=0; i<6; i++){
	for(j=0; j<itot; j++){
				xconf[i][j] = xconfi[i][j];  
			      }	
			  }*/

	//printf("0.7\n");

	/*---------------------------------MC algorithm-----------------------------*/
	
	if(selection == 4){
	select_state(sel_state,perc_seq,itot,num_mumo);
	monte_carlo_mumo(it_mc,itot,enertot,ener,xconft,xconf,sc,xconf0,stiff,enex,dd,boltz,temp,iflag, &endenergy,num_mumo,perc_seq,sel_state);
			   }
	else if (selection == 5){
	monte_carlo_cor(it_mc,itot,enertot,ener,xconft,xconf,sc,xconf0,stiff,enex,enexA,enexB,dd,dd2,dd3,dd4,dd5,boltz,temp,iflag, &endenergy);
	}		
	else if (selection == 6) {
	select_state_cor(sel_state,perc_seq,itot,num_mumo,prob_tetra_mumocor,&acc_multicor);
	monte_carlo_cor_mumo(it_mc,itot,enertot,ener,enertotH,enerH,xconft,xconf,sc,xconf0,stiff,enex,enexA,enexB,dd,dd2,dd3,dd4,dd5,boltz,temp,iflag, &endenergy,sel_state);
	}    else{
	monte_carlo(it_mc,itot,enertot,ener,xconft,xconf,sc,xconf0,stiff,enex,dd,boltz,temp,iflag, &endenergy);
	 }
					
	enertot = endenergy/0.043424; // in desired units

	/*---------------------------------End of MC algorithm-----------------------------*/

	/*for(k=0;k<itot;k++){
		for(i=0;i<6;i++){printf("%lf ", xconf[i][k]);}
				printf("\n");
			}
				printf("\n");*/
			

	//printf("enertot %lf \n", enertot);
	//printf("endenergy %lf \n", endenergy);


	//for(j=0; j<itot; j++){for(i=0; i<6; i++){printf("%lf ",xconf0[i][j]);} printf("\n");}


	/* ------ Set end structure equal to equilibrum structure ------- */
	/*for(i=0; i<6; i++){
	for(j=0; j<itot; j++){
				xconf[i][j] = xconf0[i][j];  
			      }	
			  }*/
	/* ------ END Set end structure equal to equilibrum structure ------- */

	
	diff = clock() - start;
	//printf("1.1\n");
	msec = diff * 1000 / CLOCKS_PER_SEC;
	//printf("1.2\n");


		if(selection == 0 || selection == 1){
		output_str_di(output,itot,ist,xconf,alpha);	// Output hel coord in file, for input to SCHNArP
				  }
		if(selection == 2 || selection == 3 || selection == 4 || selection == 5 || selection == 6){
		output_str_tetra(output,itot,ist,xconf,alpha);
				  }


	//Here starts the backbone reconstruction part of the code

	//printf("1.3\n");
	/* ------------- Floating particle -------------- */

	
	double br_mot[3];		// Initial values floating particle
		br_mot[0] = 0.0;
		br_mot[1] = 0.0;
		br_mot[2] = 0.0;


	//printf("1.4\n");



	//printf("1\n");

	num_bp = itot+1;

	CEHS_build(itot,ist,xconf,phi,output,xyz_float,pdb_data,pdb_tot,n7_data,&tot_num_at);



	/* -------------- Calculate Phosphate position --------------*/
	

	for(i=0;i<2*num_bp;i++){
		sscanf(pdb_data[i+1], "%s %c %c %d %lf %lf %lf\n", tmp_at[i], &baset[i], &strand[i], &base_num[i], &tmp_xyz[i][0], &tmp_xyz[i][1], &tmp_xyz[i][2]);
				}

	// Get coordinates for N7 (specific to purines A,G)
	for(i=0;i<num_bp;i++){
        	sscanf(n7_data[i], "%d %lf %lf %lf\n", &base_num_n7[i], &n7_xyz[i][0], &n7_xyz[i][1], &n7_xyz[i][2]);
				}

	//printf("2\n");
	for(i=0;i<num_bp;i++){		// Sort base number
			if(base_num_n7[i] > num_bp){
						base_num_n7[i] = 2*num_bp + 1 - base_num_n7[i];
						    }
			      }

	//printf("3\n");
	for(i=0;i<num_bp;i++){		// Copy coordinates in buffer
				n7_xyz_buf[i][0] = n7_xyz[i][0];
	
			n7_xyz_buf[i][1] = n7_xyz[i][1];
				n7_xyz_buf[i][2] = n7_xyz[i][2];
			      }

	for(i=0;i<num_bp;i++){		// Sort coordinates according to base number sorting
	for(j=0;j<num_bp;j++){
			if(base_num_n7[j] == i+1){
						n7_xyz[i][0] = n7_xyz_buf[j][0];
						n7_xyz[i][1] = n7_xyz_buf[j][1];
						n7_xyz[i][2] = n7_xyz_buf[j][2];
						  }
			      }
			      }



	for(i=0;i<num_bp;i++){
			if(strcmp(tmp_at[i],"C8") == 0){
					for(j=0;j<3;j++){
					c8_xyz[i][j] = tmp_xyz[i][j];
					c6_xyz[i][j] = tmp_xyz[2*itot+1-i][j];
							}
						}
			else{
					for(j=0;j<3;j++){
					c6_xyz[i][j] = tmp_xyz[i][j];
					c8_xyz[i][j] = tmp_xyz[2*itot+1-i][j]; 	
							}
				}
			       }


	for(i=0;i<num_bp;i++){
	for(j=0;j<3;j++){
	vec_c6_c8[i][j] = c6_xyz[i][j] - c8_xyz[i][j];
				}
				}

	for(i=0;i<num_bp;i++){
	for(j=0;j<3;j++){
	vec_n7_c8[i][j] = n7_xyz[i][j] - c8_xyz[i][j];
				}
				}

	for(i=0;i<num_bp;i++){
	cos_alpha[i] = dot_product(vec_c6_c8[i],vec_n7_c8[i]) / (vec_len(vec_c6_c8[i],3) * vec_len(vec_n7_c8[i],3)); // angle between vec(C8C6) and vec(C8N7)
				}
	
	double len_az;
	for(i=0;i<num_bp;i++){
	len_az = cos_alpha[i] * vec_len(vec_n7_c8[i],3);
			for(j=0;j<3;j++){
	vec_z[i][j] = c8_xyz[i][j] + len_az * (vec_c6_c8[i][j])/vec_len(vec_c6_c8[i],3); // C8 -> C6 direction
				}	
					}

	for(i=0;i<num_bp;i++){
	for(j=0;j<3;j++){
	vec_n7_z[i][j] = n7_xyz[i][j] - vec_z[i][j];		// getting z -> N7 vector
			}
				}

	for(i=0;i<num_bp;i++){
	for(j=0;j<3;j++){
	x_axis_bp[i][j] = vec_n7_z[i][j] / vec_len(vec_n7_z[i],3);		// normalized
			}
				}
	
	for(i=0;i<num_bp;i++){
	for(j=0;j<3;j++){
	y_axis_bp[i][j] = -1 * vec_c6_c8[i][j] / vec_len(vec_c6_c8[i],3);		// (Idk if -1 or 1) normalized
			}
				}



	/* ---- Get Phosphate positions more accurate ---- */

	double cos_x,cos_y;	
	double buf_x, buf_y;	

	for(i=0;i<num_bp-1;i++){
	cos_x = dot_product(x_axis_bp[i],x_axis_bp[i+1]);
	cos_y = dot_product(y_axis_bp[i],y_axis_bp[i+1]);

	//printf("%lf %lf %lf %lf \n", cos_x, cos(acos(cos_x)), cos_y, cos(acos(cos_y)));

	

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = 0.5*cos_x*x_axis_bp[i][j] + 0.5*x_axis_bp[i+1][j];
	mid_bp_y[i][j] = 0.5*cos_y*y_axis_bp[i][j] + 0.5*y_axis_bp[i+1][j];
			}
	
	buf_x = vec_len(mid_bp_x[i],3);
	buf_y = vec_len(mid_bp_y[i],3);

	for(j=0;j<3;j++){
	mid_bp_x[i][j] = mid_bp_x[i][j] / buf_x;
	mid_bp_y[i][j] = mid_bp_y[i][j] / buf_y;
			}

				}

	for(i=0;i<num_bp;i++){	// Center bp reconstruction from the extraction from the pdb
				center_bp[i][0] = (tmp_xyz[i][0] + tmp_xyz[2*itot+1-i][0])/2;
				center_bp[i][1] = (tmp_xyz[i][1] + tmp_xyz[2*itot+1-i][1])/2;
				center_bp[i][2] = (tmp_xyz[i][2] + tmp_xyz[2*itot+1-i][2])/2;
				}

	for(i=0;i<num_bp-1;i++){
			mid_bp[i][0] = 0.5*(center_bp[i+1][0] + center_bp[i][0]);
			mid_bp[i][1] = 0.5*(center_bp[i+1][1] + center_bp[i][1]);
			mid_bp[i][2] = 0.5*(center_bp[i+1][2] + center_bp[i][2]);		
				}


	double a1,a2,a_adj;
	a_adj = 12;		// Adjust angle values
	a1 = 33+90-a_adj;	
	a2 = 147+90+a_adj;



	k=0;
	for(i=0;i<num_bp-1;i++){
	for(j=0;j<3;j++){
	ph_pos_ex[k][j] = mid_bp[i][j] + cos(a2*PI/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a2*PI/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;

	for(j=0;j<3;j++){
	ph_pos_ex[k][j] = mid_bp[i][j] + cos(a1*PI/180)*d_ph_midbp*mid_bp_x[i][j] + sin(a1*PI/180)*d_ph_midbp*mid_bp_y[i][j];
			}
	k++;
				}


	j=0;
	for(i=0;i<num_bp-1;i++){
			if((baset[i] == 'G' && baset[i+1] == 'C') || (baset[i] == 'C' && baset[i+1] == 'C') || (baset[i] == 'G' && baset[i+1] == 'T') || (baset[i] == 'G' && baset[i+1] == 'T') || (baset[i] == 'A' && baset[i+1] == 'T') || (baset[i] == 'C' && baset[i+1] == 'T') || (baset[i] == 'T' && baset[i+1] == 'C') || (baset[i] == 'T' && baset[i+1] == 'T') || (baset[i] == 'A' && baset[i+1] == 'C')){
				base_num_ph[j] = 2 + i;
				base_num_ph[j+1] = 2*num_bp - i;
				j=j+2;
				
							}

			else{
				base_num_ph[j] = 2*num_bp - i;
				base_num_ph[j+1] = 2 + i;
				j=j+2;
 			     }
				}

	
	for(i=0;i<2*(num_bp-1);i++){
	for(j=0;j<2*num_bp;j++){
				if(base_num_ph[i] == base_num[j]){	// do same for strand and incorporate in pdb
					baset_ph[i] = baset[j];
					strand_ph[i] = strand[j];
								}
				}
				}


	// Cartesian reconstruction
	cart_rec_old(itot,no_rot_pars,xconf,pos_store,orient_store,rot_store,&e_dist,&e_dist2);
	write_bp_pos_in_table(table_bp,itot,phi,pos_store,out_folder);

	/* ---- Write into pdb ---- */	// Addition to floating algorithm

	sprintf(filnam, "%s/output_schnarp/structure_%.6d.pdb", out_folder, phi);

	wrt_my_pdb_from_pdbtot_phosphate(tot_num_at,pdb_tot,filnam,ph_pos_ex,num_bp,baset_ph,strand_ph,base_num_ph);	// More exact Ph pos


	//printf("9\n");	


		/* ---- Write values for output tables ---- */

		int ind;
	
		ind=ord[0];						// assigns which hel par is printed in table
		write_hel_coord_in_table(table_heli_shif, itot, ord, xconf, phi, ind);
		
		ind=ord[1];
		write_hel_coord_in_table(table_heli_slid, itot, ord, xconf, phi, ind);
	
		ind=ord[2];
		write_hel_coord_in_table(table_heli_rise, itot, ord, xconf, phi, ind);
	
		ind=ord[3];
		write_hel_coord_in_table(table_heli_tilt, itot, ord, xconf, phi, ind);
	
		ind=ord[4];
		write_hel_coord_in_table(table_heli_roll, itot, ord, xconf, phi, ind);
	
		ind=ord[5];
		write_hel_coord_in_table(table_heli_twis, itot, ord, xconf, phi, ind);

		write_e_dist_in_table(table_e_dist, e_dist);

		write_e_dist_in_table(table_e_dist2, e_dist2);

		write_state_in_table(table_state,itot,sel_state);
	
	//printf("10\n");
	
	
		for(i=0;i<itot+1;i++){		//Put strings in file = one line in the file
		fputs(table_heli_shif[i],file_shif);}
	
		for(i=0;i<itot+1;i++){		
		fputs(table_heli_slid[i],file_slid);}
		
		for(i=0;i<itot+1;i++){		
		fputs(table_heli_rise[i],file_rise);}
	
		for(i=0;i<itot+1;i++){		
		fputs(table_heli_tilt[i],file_tilt);}
		
		for(i=0;i<itot+1;i++){		
		fputs(table_heli_roll[i],file_roll);}
	
		for(i=0;i<itot+1;i++){		
		fputs(table_heli_twis[i],file_twis);}

		for(i=0;i<itot+1;i++){		
		fputs(table_state[i],file_states);}

         	fputs(table_e_dist[0],file_dist);
		fputs(table_e_dist2[0],file_dist2);


//printf("12\n");
} 
/* ------------------- End of for loop to execute the whole program it_prog times ------------------- */

//printf("13\n");
//
//

if (selection == 6) { printf("acceptance rate select_state structures %lf \n", acc_multicor/it_prog); } //to know the selection rate in the Ising part of the mumocor model

fclose(file_shif);
fclose(file_slid);
fclose(file_rise);
fclose(file_tilt);
fclose(file_roll);
fclose(file_twis);
fclose(file_states);
fclose(file_dist);
fclose(file_dist2);



sum_energy=0.0;
for(i=0;i<it_snapshots;i++){
sum_energy = sum_energy + ind_ener[i];
//printf("%f \n", ind_ener[i]);
}
free (ind_ener);

/* Free memory of arrays */

free(seq1);
free(ist);
free_dbl3d(stif,alc_space);
free_dbl3d(stiff,len_seq);
free(geom);
free(xconfi);
free(xconf0);
free(xconf);
free(xconft);
free(ener);
free(enerH);
free(dx);
free(output);

free(seq1_tetra);
free(seq2_tetra);
free(seq3_tetra);
free(seq4_tetra);
free(seq1_di);
free(seq1_tetra_mumo);
free_dbl3d(stif_tetra,65536);
free_dbl3d(stif_tetraS,65536);
free_dbl3d(stif_di,16);
free_dbl3d(stif_tetra_mumo,len_mumo);
free_dbl3d(prob_tetra_mumocor,7*len_seq);
free_dbl3d(prob_tetra_mumocor_all,65536);
free(geom_tetra);
free(geom_di);
free(geom_tetra_mumo);
free(perc);
free(perc_seq);
free(num_mumo);
free(read_mumo);
free(sel_state);

free(table_heli_shif);
free(table_heli_slid);
free(table_heli_rise);
free(table_heli_tilt);
free(table_heli_roll);
free(table_heli_twis);
free(table_state);
free(table_e_dist);
free(table_e_dist2);

free(out_part);
free(xyz_float);
free(q_fl);
free(center_bp);
free(pdb_data);
free(tmp_at);
free(baset);
free(strand);
free(base_num);
free(tmp_xyz);
free(d_bp);
free(d_bp_fl);
free(lj_bp);
free(lj_bp_fl);
free(mid_bp);
free(dir_bp);
free(dir_ph);
free(pos_ph);
free(d_ph1);
free(d_ph2);
free(d_ph_fl);
free(dh_ph);
free(dh_ph_fl);

free(xconf_buf);
free(xyz_float_buf);

free(pdb_tot);
free(base_num_ph);
free(baset_ph);
free(strand_ph);
free(n7_data);
free(base_num_n7);
free(n7_xyz);
free(n7_xyz_buf);
free(c6_xyz);
free(c8_xyz);
free(vec_c6_c8);
free(vec_n7_c8);
free(cos_alpha);
free(vec_z);
free(vec_n7_z);
free(x_axis_bp);
free(y_axis_bp);
free(mid_bp_x);
free(mid_bp_y);
free(ph_pos_ex);

free_dbl3d(orient_store,len_seq);
free(pos_store);
free_dbl3d(rot_store,2*len_seq);
free(table_bp);

/*
free_char2d(seq1,len_seq);
free_char2d(ist,len_seq);
free_dbl3d(stif,len_seq,no_helpars);
free_dbl3d(stiff,len_seq,no_helpars);
free_dbl2d(geom,len_seq);
free_dbl2d(xconfi,no_helpars);
free_dbl2d(xconf0,no_helpars);
free_dbl2d(xconf,no_helpars);
free_dbl2d(xconft,no_helpars);
free(ener);
free_dbl2d(dx,no_helpars);
free_char2d(output,len_seq);

free_char2d(seq1_tetra,len_seq);
free_char2d(seq1_di,len_seq);
free_dbl3d(stif_tetra,len_seq,no_helpars);
free_dbl3d(stif_di,len_seq,no_helpars);
free_dbl2d(geom_tetra,len_seq);
free_dbl2d(geom_di,len_seq);

free_char2d(table_heli_shif,len_seq);
free_char2d(table_heli_slid,len_seq);
free_char2d(table_heli_rise,len_seq);
free_char2d(table_heli_tilt,len_seq);
free_char2d(table_heli_roll,len_seq);
free_char2d(table_heli_twis,len_seq);
*/


/* End free memory */

ind_energy = sum_energy/(double) it_snapshots;



printf("Time taken for stabilizing MC: %d seconds %d milliseconds\n", msec/1000, msec%1000);
/*printf("%d steps needed to get %d structures that are %d times energy-accepted by Metropolis = %lf efficiency\n", it_mod, count_no_excl, alpha ,(double) alpha/ (double) it_mod);
printf("Time taken for one of the %d variations with %d MC steps: %d seconds %d milliseconds\n", it_snapshots, it_mcsteps, msec2/1000, msec2%1000);
printf("Time taken for reconstruction: %d seconds %d milliseconds\n", msec1/1000, msec1%1000);
printf("Mean energy after %d stabilization steps and modifying %d times by %d MC steps: %f eV\n", it_mc, it_snapshots, it_mcsteps, ind_energy);*/



diff_all = clock() - start_all;

int msec_all = diff_all * 1000 / CLOCKS_PER_SEC;
printf("Time taken for whole algorithm: %d seconds %d milliseconds\n", msec_all/1000, msec_all%1000);

}
