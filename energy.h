/* ---- Functions for calculating interaction energies ---- */

double pow2(double x)
{
  return x*x;
}


double dist (double *at1, double *at2)
{
int i;
double dist,buf;

for(i=0;i<3;i++){
		buf = pow2(at1[i]-at2[i]);
		dist = dist + buf;
		}

		dist = sqrt(dist);
return dist;
}

double vec_len(double *vec, int n)
{
int i;
double len,res;
len = 0;
for(i=0;i<n;i++){len = len + vec[i]*vec[i];}
res = sqrt(len);

return(res);
}

double dot_product (double *vec1, double *vec2)
{
double dot_pr;
dot_pr = 0;
int i;

for(i=0;i<3;i++){dot_pr = dot_pr + vec1[i]*vec2[i];}

return dot_pr;
}


double lj (double r, double k, double a, double b)
{
double lj;

lj = k*(pow((a/r),12) - pow((b/r),6));

return lj;
}

double dh (double r, double q1, double q2, double eps, double kappa, double f_pi_eps0)
{
double dh;

dh = q1*1.6*pow(10,-19) * q2*1.6*pow(10,-19) / (f_pi_eps0*eps) * 1/r * pow(10,10) * exp(-kappa*r) /(1.6*pow(10,-19));

return dh;
}

void positions(int num_bp, int num_fl, int itot, int acc, double **center_bp, double **tmp_xyz, double **d_bp, double **d_bp_fl, double **xyz_float)
{
int i,j,k;
	

for(i=0;i<num_bp;i++){	// Center bp reconstruction from the extraction from the pdb
				center_bp[i][0] = (tmp_xyz[i][0] + tmp_xyz[2*itot+1-i][0])/2;
				center_bp[i][1] = (tmp_xyz[i][1] + tmp_xyz[2*itot+1-i][1])/2;
				center_bp[i][2] = (tmp_xyz[i][2] + tmp_xyz[2*itot+1-i][2])/2;
				}




	// Get distances between bp and between bp and floating particles

	k=1;
	for(i=0;i<num_bp;i++){
	for(j=0;j<num_bp;j++){ //num_bp-acc-i
				if(i+k*acc < num_bp){
					d_bp[i][k-1] = dist(center_bp[i],center_bp[i+k*acc]);
					//printf("%d %d %d\n",num_bp,i,i+k*acc);
					k++;
						     }
			     }
				k=1;
				}
	

	for(i=0;i<num_bp;i++){
	for(j=0;j<num_fl;j++){
				d_bp_fl[i][j] = dist(center_bp[i],xyz_float[j]);
			     }
				}

}

void ph_pos(int num_bp, int itot, double d_ph_midbp, double **mid_bp, double **center_bp, double **dir_bp, double **tmp_xyz, double **dir_ph, double **pos_ph)
{
int i,j,k;


	for(i=0;i<num_bp-1;i++){
			mid_bp[i][0] = 0.5*(center_bp[i+1][0] + center_bp[i][0]);
			mid_bp[i][1] = 0.5*(center_bp[i+1][1] + center_bp[i][1]);
			mid_bp[i][2] = 0.5*(center_bp[i+1][2] + center_bp[i][2]);		
				}

	for(i=0;i<num_bp;i++){
			dir_bp[i][0] = tmp_xyz[(i)][0] - tmp_xyz[2*itot+1-(i)][0];
			dir_bp[i][1] = tmp_xyz[(i)][1] - tmp_xyz[2*itot+1-(i)][1];
			dir_bp[i][2] = tmp_xyz[(i)][2] - tmp_xyz[2*itot+1-(i)][2];
				}

	// Direction of Phosphate position from mid bp
	
	double buf1,buf2,buf3;
	for(i=0;i<num_bp;i++){
			dir_ph[i][0] = dir_bp[i][0] + dir_bp[i+1][0];
			dir_ph[i][1] = dir_bp[i][1] + dir_bp[i+1][1];
			dir_ph[i][2] = dir_bp[i][2] + dir_bp[i+1][2];
	
			buf1 = dir_ph[i][0]/vec_len(dir_ph[i],3);
			buf2 = dir_ph[i][1]/vec_len(dir_ph[i],3);
			buf3 = dir_ph[i][2]/vec_len(dir_ph[i],3);

			dir_ph[i][0] = buf1;
			dir_ph[i][1] = buf2;
			dir_ph[i][2] = buf3;
			
				}

	// Phosphate position by translating 10.25A in both directions of dir_ph

	k=0;
	for(i=0;i<num_bp-1;i++){
				
				pos_ph[k][0] = mid_bp[i][0] + d_ph_midbp*dir_ph[i][0];
				pos_ph[k][1] = mid_bp[i][1] + d_ph_midbp*dir_ph[i][1];
				pos_ph[k][2] = mid_bp[i][2] + d_ph_midbp*dir_ph[i][2];
				k++;
				pos_ph[k][0] = mid_bp[i][0] - d_ph_midbp*dir_ph[i][0];
				pos_ph[k][1] = mid_bp[i][1] - d_ph_midbp*dir_ph[i][1];
				pos_ph[k][2] = mid_bp[i][2] - d_ph_midbp*dir_ph[i][2];
				k++;
				}
}


void check_olap(int num_bp, int num_fl, int itot, int acc, double **center_bp, double **tmp_xyz, double **d_bp, double **d_bp_fl, double **xyz_float, double bp_excl_vol, double bp_fl_excl_vol, int *vol_olap)
{
int i,j,k;

	/* ---- Calculate position and distances of bp's and floating particles ----*/
	positions(num_bp,num_fl,itot,acc,center_bp,tmp_xyz,d_bp,d_bp_fl,xyz_float);
		

	/* ---- Condition to start calculating LJ and DH potentials ---- */

	int overlap;
	overlap = 0;

	//double bp_excl_vol = 10;
	//double bp_fl_excl_vol = 5;	

	for(i=0;i<num_bp;i++){
	for(j=0;j<num_bp;j++){
				if(d_bp[i][j] > 0.1 && d_bp[i][j] < bp_excl_vol){overlap = 1;}
				if(d_bp_fl[i][j] > 0.1 && d_bp_fl[i][j] < bp_fl_excl_vol){overlap = 1;}
			      }
			      }

	*vol_olap = overlap;

}

void transfer_hel_coord(int itot, double **xconf, double **xconfi)
{
int i,j;

for(i=0; i<6; i++){
for(j=0; j<itot; j++){
			xconf[i][j] = xconfi[i][j];  
		      }	
		  }

}

void initial_fl_pos(int num_fl, double shift, double br_mot[3], double **xyz_float)
{
int i;
double rx,ry,rz;

for(i=0;i<num_fl;i++){
			rx = ((double) rand() / (RAND_MAX)); 
			ry = ((double) rand() / (RAND_MAX)); 
			rz = ((double) rand() / (RAND_MAX)); 			
			xyz_float[i][0] = br_mot[0] + shift*(rx - 0.5);	// Position particles between -200 and 200
			xyz_float[i][1] = br_mot[1] + shift*(ry - 0.5);
			xyz_float[i][2] = br_mot[2] + shift*(rz - 0.5);
		}

}


void dist_ph(int num_bp, int num_fl, int acc_ph, double **pos_ph, double **d_ph1, double **d_ph2, double **d_ph_fl, double **xyz_float)
{
int k,i,j;
k=1;
	for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){ //num_bp-acc-i
				if(i+1+k*acc_ph < 2*num_bp-2){
					d_ph1[i][k-1] = dist(pos_ph[i],pos_ph[i+k*acc_ph]);
					d_ph2[i][k-1] = dist(pos_ph[i],pos_ph[i+1+k*acc_ph]);
					//printf("%d %d %d %lf\n",num_bp,i,i+1+k*acc_ph,pos_ph[i+1+k*acc_ph][0]);
					k++;
				
						     }
			     }
				k=1;
				
				}


	for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<num_fl;j++){
				d_ph_fl[i][j] = dist(pos_ph[i],xyz_float[j]);
			     }
				}



	/*for(i=500;i<510;i++){	
	printf("%lf %lf %lf %lf\n",d_ph_fl[i][0],d_ph_fl[i][1],d_ph_fl[i][2],d_ph_fl[i][3]);
	}*/

}


void calc_lj(int *k1, int *k2, int num_bp, double **d_bp, double lj_cut_low, double lj_cut_up, double **lj_bp, double lj_k, double lj_a_bp, int acc, int num_fl, double **lj_bp_fl, double **d_bp_fl, double lj_a_fl)
{
int i,j,k1_buf,k2_buf;

	// Calculate inter-bp LJ pot (for lj_cut_low < r < lj_cut_up)


	k1_buf=0;
	for(i=0;i<num_bp;i++){
	for(j=0;j<num_bp;j++){		
				if(d_bp[i][j] > lj_cut_low && d_bp[i][j] < lj_cut_up){ //printf("yes1");
						lj_bp[k1_buf][0] = lj(d_bp[i][j],lj_k,lj_a_bp,lj_a_bp);
						lj_bp[k1_buf][1] = i;
						lj_bp[k1_buf][2] = i+j*acc;
						lj_bp[k1_buf][3] = d_bp[i][j];
						k1_buf++;
						    }
				  }
			     }

	

	
	// Calculate bp - floating patricle LJ pot (for lj_cut_low < r < lj_cut_up)
	k2_buf=0;
	for(i=0;i<num_bp;i++){
	for(j=0;j<num_fl;j++){		
				if(d_bp_fl[i][j] > lj_cut_low && d_bp_fl[i][j] < lj_cut_up){ //printf("yes2");
						lj_bp_fl[k2_buf][0] = lj(d_bp_fl[i][j],lj_k,lj_a_fl,lj_a_fl);
						lj_bp_fl[k2_buf][1] = i;
						lj_bp_fl[k2_buf][2] = j;
						lj_bp_fl[k2_buf][3] = d_bp_fl[i][j];
						k2_buf++;
						    }
				  }
			     }
	
	*k1 = k1_buf;
	*k2 = k2_buf;

}

void calc_dh(int *k3, int *k4, int num_bp, double **dh_ph, double **d_ph1, double **d_ph2, double dh_cut_low, double dh_cut_up, double q_ph, double eps, double kappa, double f_pi_eps0, double acc_ph, int num_fl, double **dh_ph_fl, double **d_ph_fl, double *q_fl)
{
int i,j,k3_buf,k4_buf;

	// Calculate inter-ph DH pot for dh_cut_low < r < dh_cut_up


	k3_buf=0;
	for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){		
				if(d_ph1[i][j] > dh_cut_low && d_ph1[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[k3_buf][0] = dh(d_ph1[i][j],q_ph,q_ph,eps,kappa,f_pi_eps0);  //in ev
						dh_ph[k3_buf][1] = i;
						dh_ph[k3_buf][2] = i+j*acc_ph;
						dh_ph[k3_buf][3] = d_ph1[i][j];
						k3_buf++;
						    }
				  }
			     }

	for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<2*num_bp-2;j++){		
				if(d_ph2[i][j] > dh_cut_low && d_ph2[i][j] < dh_cut_up){ //printf("yes3");
						dh_ph[k3_buf][0] = dh(d_ph2[i][j],q_ph,q_ph,eps,kappa,f_pi_eps0);
						dh_ph[k3_buf][1] = i;
						dh_ph[k3_buf][2] = i+1+j*acc_ph;
						dh_ph[k3_buf][3] = d_ph2[i][j];
						k3_buf++;
						    }
				  }
			     }
		
	

	// Calculate bp - floating patricle LJ pot for lj_cut_low < r < lj_cut_up 
	k4_buf=0;
	for(i=0;i<2*num_bp-2;i++){
	for(j=0;j<num_fl;j++){		
				if(d_ph_fl[i][j] > dh_cut_low && d_ph_fl[i][j] < dh_cut_up){ //printf("yes4");
						dh_ph_fl[k4_buf][0] = dh(d_ph_fl[i][j],q_ph,q_fl[j],eps,kappa,f_pi_eps0);
						dh_ph_fl[k4_buf][1] = i;
						dh_ph_fl[k4_buf][2] = j;
						dh_ph_fl[k4_buf][3] = d_ph_fl[i][j];
						k4_buf++;
						    }
				  }
			     }

	*k3 = k3_buf;
	*k4 = k4_buf;

}



void unit_matrix3(double **matrix)
{
int i,j;
for(i=0;i<3;i++){
for(j=0;j<3;j++){
		if(i==j){matrix[i][j] = 1.0;}
		else {matrix[i][j] = 0.0;}
		}
		}

}

void zero_matrix_n(double **matrix, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){
		matrix[i][j] = 0.0;
		}
		}

}

void rot_x(double **rot, double phi)
{
int i,j;
rot[0][1] = rot[1][0] = rot[0][2] = rot[2][0] = 0;
rot[0][0] = 1.0;
rot[1][1] = cos(phi);
rot[2][2] = cos(phi);
rot[2][1] = sin(phi);
rot[1][2] = -sin(phi);
}

void rot_y(double **rot, double phi)
{
int i,j;
rot[0][1] = rot[1][0] = rot[1][2] = rot[2][1] = 0;
rot[1][1] = 1.0;
rot[0][0] = cos(phi);
rot[2][2] = cos(phi);
rot[0][2] = sin(phi);
rot[2][0] = -sin(phi);
}


void rot_z(double **rot, double phi)
{
int i,j;
rot[0][2] = rot[2][0] = rot[1][2] = rot[2][1] = 0;
rot[2][2] = 1.0;
rot[0][0] = cos(phi);
rot[1][1] = cos(phi);
rot[1][0] = sin(phi);
rot[0][1] = -sin(phi);
}

void mat_mult_n(double **result, double **a, double **b, int n)
{
int i,j,k;

zero_matrix_n(result,n);

for(i=0;i<n;i++){
for(j=0;j<n;j++){
for(k=0;k<n;k++){
result[i][j] = result[i][j] + a[i][k]*b[k][j];
		}
		}
		}

}

void mat_vec_mult_n(double *result, double **matrix, double *vec, int n)
{
int i,j,k;

for(i=0;i<n;i++){result[i]=0;}

for(i=0;i<n;i++){
for(k=0;k<n;k++){
result[i] = result[i] + matrix[i][k]*vec[k];
		}
		}

}

void copy_2darray(double **a, double **b, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){
		a[i][j] = b[i][j];
		}
		}

}

void transp_matr_n(double **tr_matrix, double **matrix, int n)
{
int i,j;
for(i=0;i<n;i++){
for(j=0;j<n;j++){tr_matrix[i][j] = matrix[j][i];}}

}


void cart_rec_old(int itot, int no_rot_pars, double **xconf, double **store_r_mst, double ***store_orient_t_i, double ***rot_store, double *e_dist, double *e_dist2)
{
int i,j,k,m;
const float pi = 3.14159265358979323846;
double **rotz1, **roty, **rotz2;
double g,p,tw;	//gamma,phi,twist (for end-to-end distance)
double **res_buf1, **res_buf2, **result, **res_prev;
double *transl;
double **t_mst, *r_mst, *r_bppar, *end_r_mst, *end_r_mst2;
double ete_dist_buf;


rotz1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
rotz2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
roty = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf1 = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_buf2 = dbl_aloc2d(no_rot_pars,no_rot_pars);
result = dbl_aloc2d(no_rot_pars,no_rot_pars);
res_prev = dbl_aloc2d(no_rot_pars,no_rot_pars);
transl = malloc(no_rot_pars*sizeof(double));
t_mst = dbl_aloc2d(no_rot_pars,no_rot_pars);
r_mst = malloc(no_rot_pars*sizeof(double));
r_bppar = malloc(no_rot_pars*sizeof(double));
end_r_mst = malloc(no_rot_pars*sizeof(double));
end_r_mst2 = malloc(no_rot_pars*sizeof(double));


	unit_matrix3(res_prev);		// initialize matrix as unit matrix (= T_i)
	for(i=0;i<3;i++){end_r_mst[i]=0;}
	for(i=0;i<3;i++){end_r_mst2[i]=0;}

	m=0;
	for(i=0;i<itot;i++){	// result is multiplication of all

			
			p = atan(xconf[3][i]/xconf[4][i]); //in rad (mixture of tilt and roll from Hassan&Calladine 95)
			g = xconf[4][i]/cos(p)*pi/180;	//in rad (mixture of tilt and roll from Hassan&Calladine 95)
			tw = xconf[5][i]*pi/180;	//in rad
		

			//Get Mid-step triad
			rot_y(roty,g/2);
			rot_z(rotz1,p);
			rot_z(rotz2,tw/2-p);

			

			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);
			
			copy_2darray(rot_store[m],res_buf2,3);	// Store rot matrix for mid-step triad


			mat_mult_n(t_mst,res_prev,res_buf2,3);	//get mid-step matrix out of T_i (=T_i * R_z*R_y*R_z)


			//Get r-vector i+1 from mid-step triad
		
			r_bppar[0] = xconf[0][i];	//shift
			r_bppar[1] = xconf[1][i];	//slide	
			r_bppar[2] = xconf[2][i];	//rise

			mat_vec_mult_n(r_mst,t_mst,r_bppar,3);

			double buf;
			for(j=0;j<3;j++){buf = end_r_mst[j];
					end_r_mst[j] = buf + r_mst[j];	//printf("i %d j %d end_r_mst %f\n",i,j, end_r_mst[j]);
					store_r_mst[i][j] = end_r_mst[j]; //printf("2\n"); 
					}	// This gives the end r-vector
                       
                       if (2<i && i<13) {
                       double buf2;
                       for (j=0;j<3;j++) {buf2=end_r_mst2[j];
                                          end_r_mst2[j] = buf2 + r_mst[j];
                                        }
                       } // end 2 end distance for a part of the sequence just for checking and analysis purposes



				
			//Get T_i+1 out of T_i
			rot_y(roty,g);
			rot_z(rotz1,tw/2+p);
			rot_z(rotz2,tw/2-p);

			
			mat_mult_n(res_buf1,roty,rotz1,3);
			mat_mult_n(res_buf2,rotz2,res_buf1,3);

			
			copy_2darray(rot_store[m+1],res_buf2,3);	// Store rot matrix for T_i+1
			m=m+2;
			
			
			mat_mult_n(result,res_prev,res_buf2,3);	//get T_i+1 out of T_i

			copy_2darray(res_prev,result,3);	//Copy T_i+1 in res_prev

			for(j=0;j<3;j++){
			for(k=0;k<3;k++){
					store_orient_t_i[i][j][k] = res_prev[j][k]; // store orientation matrix T_i+1
					}
					}

				}

       ete_dist_buf = vec_len(end_r_mst2,3);
        *e_dist2=ete_dist_buf;


	ete_dist_buf = vec_len(end_r_mst,3);	// End-to-end distance
        *e_dist=ete_dist_buf;


free(rotz1);
free(rotz2);
free(roty);
free(res_buf1);
free(res_buf2);
free(res_prev);
free(result);
free(transl);
free(t_mst);
free(r_mst);
free(r_bppar);
free(end_r_mst);
free(end_r_mst2);

}




