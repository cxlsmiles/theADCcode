#include "id_self_energy.hpp"
#include "integral_table.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "blas_matrix.hpp"
#include <cmath>

#define OLD_ID_CODE
#undef OLD_ID_CODE


static const double SQRT_1_2 = 0.707106781186547;
static const double SQRT_1_6 = 0.408248290463863;
static const double SQRT_3_2 = 1.22474487139159;
static int chunk = 1;

#ifdef OLD_ID_CODE
static Blas_matrix f_store[8];
#endif

inline double idSelf_energy::V1212(unsigned a, unsigned b, unsigned c, unsigned d)
{
  return table_.integral(a,c,b,d);
}





// For testing of the integral deriven code
static inline void num2rowcol(unsigned num,  unsigned &a, unsigned &b)
{
  const static double eps = 0.001;
  
  a  = floor((-1. + sqrt(1. + 8. * num ) + eps) * .5);
  b = num - a*(a+1)/2;  
  
}

inline bool idSelf_energy::get_next_fourvir(unsigned& vir1, unsigned& vir2, unsigned& vir3, unsigned& vir4, double& vi)
{
  static unsigned call_count = 0;
  call_count++;
  
  unsigned pair_col, pair_row;
  num2rowcol(call_count - 1, pair_col, pair_row);
  num2rowcol(pair_row, vir2, vir1);
  num2rowcol(pair_col, vir4, vir3);
  

  vir1 += phis_.number_occupied();
  vir2 += phis_.number_occupied();
  vir3 += phis_.number_occupied();
  vir4 += phis_.number_occupied();


  if (vir4 >= phis_.number_orbitals()) return false;

  // phys. ordering
  vi = table_.integral(vir1,vir2,vir3,vir4);
  return true;
}





idSelf_energy::idSelf_energy(Integral_table& tb, SCF_data_reader& r)
  : table_(tb), phis_(r)
{

  //Group the occupied and virtual orbitals according to symmetry;
  for(unsigned i = 0; i < phis_.number_occupied();  i++) {
    occs_[phis_.irrep(i)].push_back(i); 
  }


  for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals();  i++)
    virs_[phis_.irrep(i)].push_back(i); 

  //Group all possible 2h1p states according to symmetry
  for(unsigned int sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned int k = 0; k < phis_.number_occupied(); k++)         
      for(unsigned int l = 0; l <= k; l++) 
	for(unsigned int a = phis_.number_occupied(); a < phis_.number_orbitals(); a++) {
	  
	  if(sym != 
	     phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)),phis_.irrep(l))) 
	    continue;
	  
	  if (k != l)
	    for (unsigned type=0; type < 2; type++) {
	      Conf c = {a,k,l,type};
	      sats_[sym].push_back(c);
	    } 
	  else {
	    Conf c = {a,k,l,0};
	    sats_[sym].push_back(c);
	  }
	}
}


void idSelf_energy::static_selfenergy_4plus(Blas_matrix& sigma)
{
  
  unsigned num_occ = phis_.number_occupied();
  unsigned num_vir = phis_.number_orbitals() - num_occ;
  
  //Get the density correction, 2nd order
  rho_hole_part_2(sigma);
  rho_holeparticle_2(sigma);
  rho_particle_part_2(sigma);
 
  // Compute the dynamic self-energy contributions:
  Blas_matrix m_ak(num_vir, num_occ);
  m_ak = 0.;
  dynamic_self_energy_3_res(m_ak);
  
  //Get the density correction, 3rd order
  rho_hole_part_3_res(sigma);
  rho_particle_part_3_res(sigma);

  unsigned a, b, c, d;
  double eri;

  while (get_next_fourvir(a,b,c,d,eri)) {
    if (eri) 
      contrib_4vir(sigma, m_ak, eri, a, b, c, d);
  }
  

  rho_holeparticle_3(sigma, m_ak);
  //Get the inhomogeneities
  rho2sigma(sigma);
  //Solve for sigma 4+
  linear_eq_selfenergy(sigma);
}



// This method implements the relation between 
// a lower order density to a higher order static self-energy.
// (Takes the density correction, i.e. without zeroth order, as input)
// eq.(A25), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
void idSelf_energy::rho2sigma(Blas_matrix& rho) {
  
  Blas_matrix sigma(phis_.number_orbitals(), phis_.number_orbitals());
  sigma = 0.;

  std::vector<unsigned> orbs[8];
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) {
    orbs[sym].insert(orbs[sym].begin(),virs_[sym].begin(),virs_[sym].end());
    orbs[sym].insert(orbs[sym].begin(),occs_[sym].begin(),occs_[sym].end());
  }
  

  //eq.(A25)
  for(unsigned sym_pq = 0; sym_pq < phis_.number_irreps(); sym_pq++)
    for(unsigned p_ = 0; p_ < orbs[sym_pq].size(); p_++) 
      for(unsigned q_ = 0; q_ < orbs[sym_pq].size(); q_++) {
	
	unsigned p = orbs[sym_pq][p_];
	unsigned q = orbs[sym_pq][q_];

	double sum = 0.;

	for(unsigned sym_rs = 0; sym_rs < phis_.number_irreps(); sym_rs++)
	  for(int r_ = 0; r_ < orbs[sym_rs].size();r_++)
	    for(int s_ = 0; s_ < orbs[sym_rs].size();s_++) {

	    unsigned r = orbs[sym_rs][r_];
	    unsigned s = orbs[sym_rs][s_];

	    sum += (2. * V1212(p,r,q,s) - V1212(p,r,s,q)) * rho(r,s);
	  }
	
	sigma(p,q) = sum;
      }
  
  rho = sigma;
  
}

// This method computes the second order density, hole part
// See eqs.(A15,A33), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been simplified and spin-adapted.
void  idSelf_energy::rho_hole_part_2(Blas_matrix& rho)
{
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned k1_ = 0; k1_ <=k_; k1_++) {
	unsigned k  = occs_[sym][k_];
	unsigned k1 = occs_[sym][k1_];
	
	//eq.(A15)
	double fkk = 0.;

	for(unsigned l = 0; l < phis_.number_occupied();l++)
	  for(unsigned a = phis_.number_occupied(); a < phis_.number_orbitals();a++) {
	    unsigned sym_b =
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(l)),phis_.irrep(k));
	    
	    for(unsigned b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	      
	      fkk += (2.*V1212(a,b,k,l)*V1212(a,b,k1,l)
		      -V1212(a,b,l,k)*V1212(a,b,k1,l))
		/((phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))
		  *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k1)-phis_.energy(l)));
	    }
	  }
	
	//eq.(A33)
	rho(k,k1) -= fkk; 
	if (k != k1) rho(k1,k) -= fkk; 
      }
}

// This method computes the second order density, hole/particle part
// See eqs.(A16-18,A34), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been simplified and spin-adapted.
void idSelf_energy::rho_holeparticle_2(Blas_matrix& rho)
{
  
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned a_ = 0; a_ < virs_[sym].size();a_++) {
	unsigned k = occs_[sym][k_];
	unsigned a = virs_[sym][a_];
	// eq.(A17)	
	double m_ak_p = 0.;

	for(int j = 0; j < phis_.number_occupied();j++)
	  for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(j)),phis_.irrep(b));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      
	      m_ak_p += ( 2.*V1212(a,j,b,c)*V1212(b,c,k,j) -
			  V1212(a,j,c,b)*V1212(b,c,k,j) )
		/(phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(j));
	      
	    }
	  }
	// eq.(A18)
	double m_ak_m = 0.;

	for(int i = 0; i < phis_.number_occupied();i++)
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_b = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(i),phis_.irrep(j)),phis_.irrep(k));
	    
	    for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	      
	      m_ak_m += ( 2.*V1212(i,j,k,b)*V1212(a,b,i,j) -
			  V1212(i,j,b,k)*V1212(a,b,i,j) )
		/(phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j));
	    }
	  }
	//eq.(A16)
	double fka = (m_ak_m - m_ak_p) / (phis_.energy(k)  - phis_.energy(a));
	//eq.(34)
	rho(k,a) += fka;
	rho(a,k) += fka;
      }
}

// This method computes the second order density, particle part
// See eqs.(A27,A36), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been spin-adapted, according to the table:
// <akl| 2h,1p states doublet spin eigenfunctions:
//   bab   bba   aaa   coef
//   -1     1     0   1/sqrt2 for k <> l
//    1     1     2   1/sqrt6 
//    1     0     0   1       for l = k
void idSelf_energy::rho_particle_part_2(Blas_matrix& rho)
{
  for(unsigned int sym = 0; sym < phis_.number_irreps(); sym++) {
    
    Blas_matrix f(sats_[sym].size(),virs_[sym].size());
    f = 0.;

    for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++)
      for(unsigned s_ = 0; s_ < sats_[sym].size(); s_++) {
	unsigned b = virs_[sym][b_];
	unsigned a = sats_[sym][s_].a;
	unsigned k = sats_[sym][s_].k;
	unsigned l = sats_[sym][s_].l;
	unsigned type = sats_[sym][s_].type;
	
	//eq.(A27)
	if (k == l)
	  f(s_,b_) = -V1212(a,b,l,k)
	    / (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
	else if (type == 0)
	  f(s_,b_) = -SQRT_1_2 * (V1212(a,b,l,k) + V1212(a,b,k,l))
	    /(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
	else
	  f(s_,b_) = SQRT_1_6 * (3.*V1212(a,b,l,k) - 3.*V1212(a,b,k,l))
	    /(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l)); 
      }
    
    Blas_matrix temp(virs_[sym].size(),virs_[sym].size());
    temp = 0.;
    //eq.(A36)
    temp.dgemm('T','N',f,f);
    //Reorder the matrix according to the original orbital numbering:
    for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++)
      for(unsigned a_ = 0; a_ < virs_[sym].size(); a_++) {
	unsigned b = virs_[sym][b_];
	unsigned a = virs_[sym][a_];
	rho(b,a) += temp(b_,a_);
      }
  }
}

// This method computes the third order density, hole part
// See eqs.(A19-A23,A33), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been simplified, spin-adapted and the dummy indexes swapped where necessary.
inline void idSelf_energy::rho_hole_part_3_res(Blas_matrix& rho) 
{
  
  for (unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned k1_ = 0; k1_ < occs_[sym].size() ; k1_++) {
	unsigned k = occs_[sym][k_];
	unsigned k1 = occs_[sym][k1_];
	
	//double fA = 0.;
	double fB = 0.; 
	double fC = 0.;   
	double fD = 0.;
	
	
	for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++) 
	  for(int m = 0; m < phis_.number_occupied();m++) {
	    
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(m)),phis_.irrep(k));

	    for(int c_ = 0; c_ < virs_[sym_c].size();c_++) {
	      int c = virs_[sym_c][c_];
	      
	      // Compute intermediate data for  fA-D:
	      double pole = 1. / ((phis_.energy(a)+phis_.energy(c)-phis_.energy(k1)-phis_.energy(m))
				  *(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(m)));
	      double v_ack1m = V1212(a,c,k1,m) * pole;
	      double v_acmk1 = V1212(a,c,m,k1) * pole;
	      double exp1 = v_acmk1 - 2.*v_ack1m;
	      double exp2 = v_ack1m - 2.*v_acmk1;
	      double exp3 = -2. * exp1;
	      double exp4 = -2. * exp2;
	      
// 	      // Eq.(A20)
// 	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
// 		unsigned sym_d = 
// 		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(k)),phis_.irrep(m));
// 		for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
// 		  int d = virs_[sym_d][d_];
		  
// 		  fA += V1212(d,b,k,m) * (V1212(a,c,d,b) * exp1 + V1212(a,c,b,d) * exp2) 
// 		    /(phis_.energy(k)+phis_.energy(m)-phis_.energy(d)-phis_.energy(b));
// 		}
// 	      }

	      // Eq.(A21)	      
	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
		unsigned sym_l = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(k));
		for(int l_ = 0; l_ < occs_[sym_l].size(); l_++) {
		  unsigned l = occs_[sym_l][l_];
		  
		  double v_lcbm = V1212(l,c,b,m);
		  double v_lcmb = V1212(l,c,m,b);
		    
		  fB += ( V1212(a,b,k,l) * (v_lcbm * exp3 + v_lcmb * exp1) + 
			  V1212(a,b,l,k) * (v_lcbm * exp1 + v_lcmb * exp2) )
		    /(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
		}
	      }
	      // Eq.(A22)	      	

	      for(int j = 0; j < phis_.number_occupied();j++) {
		unsigned sym_l = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(c)),phis_.irrep(j));
		for(int l_ = 0; l_ < occs_[sym_l].size(); l_++) {
		  unsigned l = occs_[sym_l][l_];
		  
		  fC += V1212(a,c,j,l) * ( V1212(j,l,m,k) * exp2 + V1212(j,l,k,m) * exp1)
		    /(phis_.energy(j)+phis_.energy(l)-phis_.energy(a)-phis_.energy(c));
		  
		}
	      }
	      // Eq.(A23)	      	

	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
		unsigned sym_l = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(a)),phis_.irrep(m));
		for(int l_ = 0; l_ < occs_[sym_l].size(); l_++) {
		  unsigned l = occs_[sym_l][l_];
		  
		  double v_lckb = V1212(l,c,k,b);
		  double v_lcbk = V1212(l,c,b,k);
		  
		  fD += ( V1212(b,a,l,m) * (v_lckb * exp2 + v_lcbk * exp4) + 
			  V1212(b,a,m,l) * (v_lckb * exp1 + v_lcbk * exp2) )
		    /(phis_.energy(b)+phis_.energy(a)-phis_.energy(l)-phis_.energy(m));
		}		 
	      }
	      
	      
	    }
	  }
	
	//Eq.(A19)
	double fkk = 0.5 * fC + (fB + fD);
	//Eq.(A33)
	rho(k,k1) += fkk; 
	rho(k1,k) += fkk;
      }
}




// This method computes the third order density, hole/particle part,
// but only with the dynamic self energy contributions
// See eqs.(A24,A34), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
void idSelf_energy::rho_holeparticle_3(Blas_matrix& rho, Blas_matrix& self_energy)
{ 

  unsigned num_occ = phis_.number_occupied();
  //eq.(A24)
  for (unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned a_ = 0; a_ < virs_[sym].size();a_++) {
	unsigned k = occs_[sym][k_];
	unsigned a = virs_[sym][a_];
      
	double fka = self_energy(a - num_occ,k) / (phis_.energy(k) - phis_.energy(a));
	rho(k,a) += fka;
	rho(a,k) += fka;
    }
}


// This method computes the third order density, particle part
// See eqs.(A27-A32,A36), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been spin-adapted (see rho_particle_part_2) and simplified
void idSelf_energy::rho_particle_part_3_res(Blas_matrix& rho)
{

  for(unsigned int sym = 0; sym < phis_.number_irreps(); sym++) {
    
    Blas_matrix f1(sats_[sym].size(),virs_[sym].size());
    f1 = 0.;
    
    for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++)
      for(unsigned s_ = 0; s_ < sats_[sym].size(); s_++) {
	unsigned b = virs_[sym][b_];
	unsigned a = sats_[sym][s_].a;
	unsigned k = sats_[sym][s_].k;
	unsigned l = sats_[sym][s_].l;
	unsigned type = sats_[sym][s_].type;
	
	
	// eq.(A27)
	if (k == l)
	  f1(s_,b_) = V1212(a,b,l,k)
	    /(phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(b));
	
	else if (type == 0)
	  f1(s_,b_) = SQRT_1_2 * (V1212(a,b,l,k) + V1212(a,b,k,l))
	    /(phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(b));
	
	else
	  f1(s_,b_) = SQRT_3_2 * (V1212(a,b,l,k) - V1212(a,b,k,l))
	    /(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
      }
#ifdef OLD_ID_CODE    
    f_store[sym].allocate(f1.rows(), f1.cols());
    f_store[sym] = f1;
#endif
    
    Blas_matrix f2(sats_[sym].size(),virs_[sym].size());
    f2 = 0.;

    
    for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++)
      for(unsigned s_ = 0; s_ < sats_[sym].size(); s_++) {
	unsigned b = virs_[sym][b_];
	unsigned a = sats_[sym][s_].a;
	unsigned k = sats_[sym][s_].k;
	unsigned l = sats_[sym][s_].l;
	unsigned type = sats_[sym][s_].type;
	
	double gA = 0.;
	double gB = 0.;
	double gC = 0.;
	double gD = 0.;
	double gC1 = 0.;
	double gD1 = 0.;

#pragma omp parallel
	{
	  if (k == l) {
	    
	    //eq. (A29)
#pragma omp for reduction(+:gA) schedule(static,chunk) nowait
	    for(int i = 0; i < phis_.number_occupied();i++) {
	      
	    unsigned sym_j = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(i));
	    
	    for(int j_ = 0; j_ < occs_[sym_j].size(); j_++) {
	      unsigned j = occs_[sym_j][j_];
	      
	      gA += V1212(a,b,i,j) * V1212(i,j,l,k)
		/(phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j));
	    }
	    
	    }
	    //eq. (A30)
// #pragma omp for reduction(+:gB)	 schedule(static,chunk) nowait
// 	    for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
// 	      unsigned sym_d = 
// 		phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	      
// 	      for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
// 		unsigned d = virs_[sym_d][d_];
		
// 		gB += V1212(c,d,l,k) * V1212(a,b,c,d)
// 		  / (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
// 	      }
// 	    }
	    // eq.(A31,A32)
#pragma omp for reduction(+:gC,gD) schedule(static,chunk)  nowait
	    for(int j = 0; j < phis_.number_occupied();j++) {
	      unsigned sym_c = 
		phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(j)),phis_.irrep(l));
	      
	      for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
		unsigned c = virs_[sym_c][c_];
		
		gC += V1212(b,c,j,l) * V1212(a,j,c,k) 
		  / (phis_.energy(j)+phis_.energy(l)-phis_.energy(b)-phis_.energy(c));
		
		gD += V1212(a,c,j,k) * V1212(b,j,c,l) 
		  / (phis_.energy(j)+phis_.energy(k)-phis_.energy(a)-phis_.energy(c));
	      }
	    } 
	    // eq.(A31,A32)	  
#pragma omp for reduction(+:gC1,gD1) schedule(static,chunk)  
	    for(int j = 0; j < phis_.number_occupied();j++) {
	      unsigned sym_c = 
		phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(j)),phis_.irrep(l));
	      
	      for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
		unsigned c = virs_[sym_c][c_];
		
		double v_bckj = V1212(b,c,k,j);
		
		gC1 += ( (2.*v_bckj - V1212(b,c,j,k)) * V1212(a,j,l,c)
			-v_bckj * V1212(a,j,c,l) ) 
		  / (phis_.energy(b)+phis_.energy(c)-phis_.energy(j)-phis_.energy(k));
		
		double v_bjkc = V1212(b,j,k,c);
		
		gD1 += (V1212(a,c,l,j)*(2.*v_bjkc-V1212(b,j,c,k))
		       -V1212(a,c,j,l)*v_bjkc)
		  / (phis_.energy(a)+phis_.energy(c)-phis_.energy(j)-phis_.energy(l));
	      }
	    }
	  } else if (type == 0) {
	    // eq.(A29)
#pragma omp for reduction(+:gA)  schedule(static,chunk) nowait
	  for(int i = 0; i < phis_.number_occupied();i++) {
	    
	    unsigned sym_j = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(i));
	    
	    for(int j_ = 0; j_ < occs_[sym_j].size(); j_++) {
	      unsigned j = occs_[sym_j][j_];
	      
	      gA += V1212(a,b,i,j) * (V1212(i,j,l,k) + V1212(i,j,k,l))
		/ (phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j));
	    }
	  }
	  // eq.(A30)
// #pragma omp for reduction(+:gB)  schedule(static,chunk) nowait
// 	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
//  	    unsigned sym_d = 
// 	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	    
// 	    for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
// 	      unsigned d = virs_[sym_d][d_];
	      
// 	      gB += V1212(a,b,c,d) * (V1212(c,d,l,k) + V1212(c,d,k,l))
// 		/ (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
// 	    }
// 	  }
	  // eq.(A31,A32)
#pragma omp for reduction(+:gC,gD)  schedule(static,chunk) nowait
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(j)),phis_.irrep(l));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      
	      double v_ajkc = V1212(a,j,k,c);
	      double v_ajck = V1212(a,j,c,k);
	      
	      gC += ( V1212(b,c,l,j) * (2.*v_ajkc - v_ajck)
		      -V1212(b,c,j,l) * (v_ajck + v_ajkc) ) 
		/ (phis_.energy(b)+phis_.energy(c)-phis_.energy(j)-phis_.energy(l));
	      
	      
	      double v_ackj = V1212(a,c,k,j);
	      double v_acjk = V1212(a,c,j,k);
	      
	      gD += ( (2.*v_ackj - v_acjk) * V1212(b,j,l,c)
		      -(v_acjk + v_ackj) * V1212(b,j,c,l))
		/ (phis_.energy(a)+phis_.energy(c)-phis_.energy(j)-phis_.energy(k));
	    }
	  }
	  // eq.(A31,A32)
#pragma omp for reduction(+:gC1,gD1)  schedule(static,chunk) 
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(j)),phis_.irrep(l));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      
	      double v_bckj = V1212(b,c,k,j);
	      double v_bcjk = V1212(b,c,j,k);
	      
	      gC1 += ( (2.*v_bckj - v_bcjk) * V1212(a,j,l,c)
		      -(v_bckj + v_bcjk) * V1212(a,j,c,l) )
		/ (phis_.energy(b)+phis_.energy(c)-phis_.energy(j)-phis_.energy(k));
	      
	      double v_aclj = V1212(a,c,l,j);
	      double v_acjl = V1212(a,c,j,l);
	      
	      gD1 += ( (2.*v_aclj - v_acjl) * V1212(b,j,k,c)
		      -(v_acjl + v_aclj) * V1212(b,j,c,k) )
		/ (phis_.energy(a)+phis_.energy(c)-phis_.energy(j)-phis_.energy(l));
	      
	    }
	    
	  }
#pragma omp single
	    {
	      gA *= SQRT_1_2;
	      gB *= SQRT_1_2;
	      gC *= SQRT_1_2;
	      gD *= SQRT_1_2;
	      gC1 *= SQRT_1_2;
	      gD1 *= SQRT_1_2;

	    }
	  } else {
	    // eq.(29)
#pragma omp for reduction(+:gA) schedule(static,chunk)  nowait
	    for(int i = 0; i < phis_.number_occupied();i++) {
	      
	      unsigned sym_j = 
		phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(i));
	      
	      for(int j_ = 0; j_ < occs_[sym_j].size(); j_++) {
		unsigned j = occs_[sym_j][j_];
		
		gA += V1212(a,b,i,j)*(V1212(i,j,k,l) - V1212(i,j,l,k))
		  / (phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j));
	      }
	    }
	    // eq.(A30)
// #pragma omp for reduction(+:gB) schedule(static,chunk)  nowait    
// 	    for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
// 	      unsigned sym_d = 
// 		phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	      
// 	      for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
// 		unsigned d = virs_[sym_d][d_];
		
		
// 		gB += (V1212(c,d,k,l) - V1212(c,d,l,k)) * V1212(a,b,c,d)
// 		  / (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
// 	      }
// 	    }
	    // eq.(A31,A32)
#pragma omp for reduction(+:gC,gD) schedule(static,chunk)  nowait
	    for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(j)),phis_.irrep(l));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      
	      
	      double v_ajck = V1212(a,j,c,k);
	      double v_ajkc = V1212(a,j,k,c);
	      
	      gC += ( V1212(b,c,j,l)*(v_ajck - v_ajkc)
		      +V1212(b,c,l,j)*(2.*v_ajkc - v_ajck) )
		/ (phis_.energy(b)+phis_.energy(c)-phis_.energy(j)-phis_.energy(l));
	      
	      double v_acjk = V1212(a,c,j,k);
	      double v_ackj = V1212(a,c,k,j);
	      
	      gD += ( (v_acjk - v_ackj) * V1212(b,j,c,l)
		      +(2.*v_ackj - v_acjk) * V1212(b,j,l,c))
		/ (phis_.energy(a)+phis_.energy(c)-phis_.energy(j)-phis_.energy(k));
	    }
	    } 
	    // eq.(A31,A32)
#pragma omp for reduction(+:gC1,gD1) schedule(static,chunk) 
	    for(int j = 0; j < phis_.number_occupied();j++) {
	      unsigned sym_c = 
		phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(j)),phis_.irrep(l));
	      
	      for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
		unsigned c = virs_[sym_c][c_];
		
		double v_ajcl = V1212(a,j,c,l);
		double v_ajlc = V1212(a,j,l,c);
		
		gC1 +=  ( V1212(b,c,j,k) * (v_ajcl - v_ajlc)
			 +(2.*v_ajlc - v_ajcl) * V1212(b,c,k,j))
		  / (phis_.energy(j)+phis_.energy(k)-phis_.energy(b)-phis_.energy(c));
		
		double v_acjl = V1212(a,c,j,l);
		double v_aclj = V1212(a,c,l,j);
		
		gD1 += ( (v_acjl - v_aclj) * V1212(b,j,c,k)
		      +(2.*v_aclj - v_acjl) * V1212(b,j,k,c) )
		  / (phis_.energy(j)+phis_.energy(l)-phis_.energy(a)-phis_.energy(c));
		
	      }
	    }
#pragma omp single
	    {
	      gA *= SQRT_3_2;
	      gB *= SQRT_3_2;
	      gC *= SQRT_3_2;
	      gD *= SQRT_3_2;
	      gC1 *= SQRT_3_2;
	      gD1 *= SQRT_3_2;

	    }
	  }
	  
	}
	f2(s_,b_) = (gA+gB+gC+gD+gC1+gD1)/(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
      }
    
    Blas_matrix temp(virs_[sym].size(),virs_[sym].size());
    temp = 0.;
    
    temp.dgemm('T','N',f1,f2);
    temp.dgemm('T','N',f2,f1);
    
    for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++)
      for(unsigned a_ = 0; a_ < virs_[sym].size(); a_++) {
	unsigned b = virs_[sym][b_];
	unsigned a = virs_[sym][a_];
	rho(b,a) += temp(b_,a_);
      }
  }
}







// This method computes the third order static self energy
// contributions eq.(C.6-C.11)in Appendix C, W. von Niessen et al.,
// Computer Physics Reports 1 (1984) 57-125.
void idSelf_energy::sigma_3diagrams(Blas_matrix& sigma)
{
  std::vector<unsigned> orbs[8];
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) {
    orbs[sym].insert(orbs[sym].begin(),virs_[sym].begin(),virs_[sym].end());
    orbs[sym].insert(orbs[sym].begin(),occs_[sym].begin(),occs_[sym].end());
  }
  
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned p_ = 0; p_ < orbs[sym].size(); p_++) 
      for(unsigned q_ = 0; q_ <= p_; q_++) {
	
	unsigned p = orbs[sym][p_];
	unsigned q = orbs[sym][q_];

	double A1pq = 0.; double A2pq = 0.; double A3pq = 0.;
	double A4pq = 0.; double A5pq = 0.; double A6pq = 0.;
	
	for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++)
	  for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++)
	    for(int i = 0; i < phis_.number_occupied();i++)
	      for(int j = 0; j < phis_.number_occupied();j++) {
		
		if (phis_.irrep(a) != phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(i)),phis_.irrep(j))) 
		  continue;
	      
		double term = (2. * V1212(j,i,a,b) - V1212(j,i,b,a)) /
		  (phis_.energy(j)+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
		
		for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) { 
		  
		  if ( phis_.irrep(i) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(j),phis_.irrep(c)),phis_.irrep(a))
		       && phis_.irrep(p) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(c),phis_.irrep(b)),phis_.irrep(q))) 
		    
		    A2pq += V1212(i,j,c,a) *
		      (2. * V1212(p,c,q,b) - V1212(p,c,b,q)) *
		      term /
		      (phis_.energy(j)+phis_.energy(i)-phis_.energy(a)-phis_.energy(c));
		
		  if ( phis_.irrep(a) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(c)),phis_.irrep(i))
		       && phis_.irrep(p) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(c),phis_.irrep(q)),phis_.irrep(j))) {
		    
		    double term1 = V1212(a,b,c,i) * term / (phis_.energy(j)-phis_.energy(c));
		  
		    A3pq += (2. * V1212(p,c,q,j) - V1212(p,c,j,q)) 
		    * term1;
		    
		    A4pq += (2. * V1212(p,j,q,c) - V1212(p,j,c,q)) 
		      * term1;
		  }
		}
		
		for(int k = 0; k < phis_.number_occupied();k++) {
		  
		  if ( phis_.irrep(a) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(k)),phis_.irrep(i))
		       && phis_.irrep(p) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(k),phis_.irrep(q)),phis_.irrep(j))) 
		    
		    A1pq += V1212(a,b,k,i) *
		      (2. * V1212(p,k,q,j) - V1212(p,k,j,q)) *
		      term / (phis_.energy(k)+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
		  
		  if ( phis_.irrep(i) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(j),phis_.irrep(k)),phis_.irrep(a))
		       && phis_.irrep(p) == phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(q)),phis_.irrep(k))) {
		    
		    double term1 =  V1212(i,j,k,a) * term / (phis_.energy(k)-phis_.energy(b));
		    
		    A5pq += (2. * V1212(p,b,q,k) - V1212(p,b,k,q)) 
		      * term1;
		    
		    A6pq += (2. * V1212(p,k,q,b) - V1212(p,k,b,q))
		      * term1;
		  }
		}
	      }
	
      
	sigma(p,q) = -A1pq + A2pq + A3pq + A4pq - A5pq - A6pq;
      
	sigma(q,p) = sigma(p,q);
      
    }
  
}


void idSelf_energy::dynamic_self_energy_3_res(Blas_matrix& M_ak)
{
  unsigned num_occ = phis_.number_occupied();

  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned p_ = 0; p_ < virs_[sym].size(); p_++) 
      for(unsigned q_ = 0; q_ < occs_[sym].size(); q_++) {
      unsigned p = virs_[sym][p_];
      unsigned q = occs_[sym][q_];
      
      //cout << p << " " << q << endl;
      double e_p = phis_.energy(p);
      double e_q = phis_.energy(q);

      double C1pq = 0.; double C2pq = 0.; double C3pq = 0.;
      double C4pq = 0.; double C5pq = 0.; double C6pq = 0.;

      double D1pq = 0.; double D2pq = 0.; double D3pq = 0.;
      double D4pq = 0.; double D5pq = 0.; double D6pq = 0.;

      // eqs. (C12,C13)
#pragma omp parallel 
      {
#pragma omp for reduction(+:C1pq, C2pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied(); i++) 
	  for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++) {
	    unsigned sym_b = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(a));
	    
	    for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	    
	      double v_piab = (2. * V1212(p,i,a,b) - V1212(p,i,b,a))
		/ (e_q+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
	      
// 	      for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
// 		unsigned sym_d = 
// 		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(c));
		
// 		for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
// 		  unsigned d = virs_[sym_d][d_];
		  
// 		  C1pq += v_piab * V1212(a,b,c,d) * V1212(q,i,c,d)
// 		    / (e_q+phis_.energy(i)-phis_.energy(c)-phis_.energy(d));
// 		}
// 	      }
	      
	      
	      for(int j = 0; j < phis_.number_occupied();j++) {
		unsigned sym_k = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(j));
		
		for(int k_ = 0; k_ < occs_[sym_k].size(); k_++) {
		  unsigned k = occs_[sym_k][k_];
		  
		  C2pq += v_piab * V1212(a,b,j,k) * V1212(q,i,j,k)
		    / (phis_.energy(j)+phis_.energy(k)-phis_.energy(a)-phis_.energy(b));
		}
	      }
	    }
	  }
	// eq.(C.14)
#pragma omp for reduction(+:C3pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied(); i++) 
	  for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++) {
	    unsigned sym_b = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(a));
	    
	    for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	      
	      double v_qiab = V1212(q,i,a,b)
		/ (e_q+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
	      
	      for(int k = 0; k < phis_.number_occupied();k++) {
		unsigned sym_j = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(k));
		
		for(int j_ = 0; j_ < occs_[sym_j].size(); j_++) {
		  unsigned j = occs_[sym_j][j_];
		  
		  C3pq += (2. * V1212(p,i,j,k) - V1212(p,i,k,j)) * V1212(a,b,j,k) * v_qiab
		    / (phis_.energy(j)+phis_.energy(k)-phis_.energy(a)-phis_.energy(b));
		}
	      }
	    }
	  } 
	// eqs. (C.15,C.17)
#pragma omp for reduction(+:C4pq, C6pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied();i++)
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(j));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      double v_paij = (2. * V1212(p,a,i,j) - V1212(p,a,j,i))
		/ (e_p+phis_.energy(a)-phis_.energy(i)-phis_.energy(j));
	      
	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
		unsigned sym_c = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),sym_a),phis_.irrep(b));
		
		for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
		  unsigned c = virs_[sym_c][c_];
		  
		  C4pq += v_paij * V1212(i,j,b,c) * V1212(q,a,b,c)
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
		}
	      }
	      
	      for(int k = 0; k < phis_.number_occupied();k++) {
		unsigned sym_l = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),sym_a),phis_.irrep(k));
		
		for(int l_ = 0; l_ < occs_[sym_l].size(); l_++) {
		  unsigned l = occs_[sym_l][l_];
		  
		  C6pq += v_paij * V1212(i,j,k,l) * V1212(q,a,k,l)
		    / (e_p+phis_.energy(a)-phis_.energy(k)-phis_.energy(l));
		  
		}
	      }		
	    }
	  }
// 	// eq.(C.16)
// #pragma omp for reduction(+:C5pq) schedule(static,chunk)  nowait
// 	for(int i = 0; i < phis_.number_occupied();i++)
// 	  for(int j = 0; j < phis_.number_occupied();j++) {
// 	    unsigned sym_a = 
// 	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(i)),phis_.irrep(j));
	    
// 	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
// 	      unsigned a = virs_[sym_a][a_];
	      
// 	      double v_qaij = V1212(q,a,i,j)
// 		/ (e_p+phis_.energy(a)-phis_.energy(i)-phis_.energy(j));
	      
// 	      for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
// 		unsigned sym_b = 
// 		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),sym_a),phis_.irrep(c));
		
// 		for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
// 		  unsigned b = virs_[sym_b][b_];
		  
// 		  C5pq += (2. * V1212(p,a,b,c) - V1212(p,a,c,b)) * V1212(i,j,b,c) * v_qaij
// 		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
// 		}
// 	      }
// 	    }
// 	  }
	
	
	//eqs.(C.18,C.20)
#pragma omp for reduction(+:D1pq, D3pq) schedule(static,chunk)  nowait
	for(int j = 0; j < phis_.number_occupied();j++) 
	  for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++){
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(j)),phis_.irrep(b));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      //intermediates for D1 and D3
	      double pole = 1. / (e_q + phis_.energy(j) - phis_.energy(b) - phis_.energy(c));
	      double v_qjcb = V1212(q,j,c,b) * pole;
	      double v_qjbc = V1212(q,j,b,c) * pole;
	      double exp1 = v_qjcb - 2. * v_qjbc;
	      double exp2 = v_qjbc - 2. * v_qjcb;
	      double exp3 = -2. * exp1;//( 4. * V1212(q,j,b,c) - 2. * V1212(q,j,c,b) ) * pole;
	      
	      for(int i = 0; i < phis_.number_occupied();i++) {
		unsigned sym_a = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(j),phis_.irrep(i)),sym_c);
		
		for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
		  unsigned a = virs_[sym_a][a_];
		  
		  double v_ajic = V1212(a,j,i,c);
		  double v_ajci = V1212(a,j,c,i);
		  
		  D1pq += (   V1212(p,i,a,b) * 
			      (v_ajic * exp1 + v_ajci * exp2) +
			      V1212(p,i,b,a) * 
			      (v_ajic * exp3 + v_ajci * exp1)  )  /
		    (e_q + phis_.energy(i) - phis_.energy(a) - phis_.energy(b));
		  
		  double v_ijac = V1212(i,j,a,c);
		  double v_ijca = V1212(i,j,c,a);
		  
		  D3pq += (   V1212(p,a,i,b) * 
			      (v_ijac * exp1 + v_ijca * exp2) +
			      V1212(p,a,b,i) * 
			      (v_ijac * exp3 + v_ijca * exp1)  )  /
		    (phis_.energy(i)+phis_.energy(j)-phis_.energy(c)-phis_.energy(a));
		}
	      }
	    }
	  }
	//eq.(C.19)
#pragma omp for reduction(+:D2pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied();i++) 
	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(c));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      double pole = 1. / (e_q + phis_.energy(i) - phis_.energy(a) - phis_.energy(c));
	      double v_pica = V1212(p,i,c,a) * pole;
	      double v_piac = V1212(p,i,a,c) * pole;
	      
	      for(int j = 0; j < phis_.number_occupied();j++) {
		unsigned sym_b = 
		  phis_.irrep_product(phis_.irrep_product(sym_a,phis_.irrep(i)),phis_.irrep(j));
		
		for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
		  unsigned b = virs_[sym_b][b_];
		  
		  double v_abij = V1212(a,b,i,j);
		  double v_abji = V1212(a,b,j,i);
		  double v_qbjc = V1212(q,b,j,c);
		  double v_qbcj = V1212(q,b,c,j);
		  double exp1 = v_qbjc - 2. * v_qbcj;
		  double exp2 = v_qbcj - 2. * v_qbjc;
		  double exp3 = -2. * exp1;
		  
		  D2pq += (  v_pica * ( v_abij * exp3 + 
					v_abji * exp1 ) +
			     v_piac * ( v_abij * exp1 +
					v_abji * exp2 )  ) /
		    (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(b));
		}
	      }
	    }
	  }
	//eqs(C.21,C.22)
#pragma omp for reduction(+:D6pq, D4pq, D5pq) schedule(static,chunk)  nowait
	for(int j = 0; j < phis_.number_occupied();j++) 
	  for(int k = 0; k < phis_.number_occupied();k++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(j)),phis_.irrep(k));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      // intermediates for D5 and D6
	      double pole = 1. / (e_p+phis_.energy(a)-phis_.energy(j)-phis_.energy(k));
	      double v_qajk = V1212(q,a,j,k) * pole;
	      double v_qakj = V1212(q,a,k,j) * pole;
	      double exp1 = v_qajk - 2. * v_qakj;
	      double exp2 = v_qakj - 2. * v_qajk;
	      double exp3 = -2. * exp1; //(4. * V1212(q,a,k,j) - 2. * V1212(q,a,j,k)) *pole;
	      //intermediates for D4
	      double v_pakj = V1212(p,a,k,j) * pole;
	      double v_pajk = V1212(p,a,j,k) * pole;
	      
	      for(int i = 0; i < phis_.number_occupied();i++) {
		unsigned sym_b = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(j),phis_.irrep(i)),sym_a);
		
		for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
		  unsigned b = virs_[sym_b][b_];
		  
		  double v_jiab = V1212(j,i,a,b);
		  double v_jiba = V1212(j,i,b,a);
		  
		  D5pq += (V1212(p,i,b,k) * (v_jiab * exp1 + v_jiba * exp2)
			   +V1212(p,i,k,b) * (v_jiab * exp3 + v_jiba * exp1))
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(b));
		  
		  double v_iabj = V1212(i,a,b,j);
		  double v_iajb = V1212(i,a,j,b);
		  
		  D6pq += (V1212(p,b,k,i) * (v_iabj * exp3 + v_iajb * exp1)
			   +V1212(p,b,i,k) * (v_iabj * exp1 + v_iajb * exp2) )
		    / (e_p+phis_.energy(b)-phis_.energy(i)-phis_.energy(k));
		  
		  double v_qibk = V1212(q,i,b,k);
		  double v_qikb = V1212(q,i,k,b);
		  double exp4 = v_qibk - 2. * v_qikb;
		  double exp5 = v_qikb - 2. * v_qibk;
		  double exp6 = -2. * exp4;
		  
		  D4pq += (v_pakj * (v_jiab * exp6
				     + v_jiba * exp4)
			   +v_pajk * (v_jiab * exp4
				      + v_jiba * exp5))
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(b));
		  
		}
	      }
	    }
	  }
      }
      M_ak(p-num_occ,q) = 
	(C1pq + C2pq) + (C3pq + C4pq) + (C5pq - C6pq) +
	(D1pq + D2pq) + (D3pq + D4pq) + (D5pq - D6pq);
      
      }  
}
 
 
 
// This method takes as input the inhomogeneity matrix B
// and produces the static self-interaction using the
// linear equation defined in Appendix B, W. von Niessen et al.,
// Computer Physics Reports 1 (1984) 57-125.
void idSelf_energy::linear_eq_selfenergy(Blas_matrix& b)  
{
  // Find the number of particle-hole pairs (1)
  // and hole-hole, particle-particle pairs (2)
  // skip all pairs leading to zero elements due to symmetry.
  unsigned int ph_count = 0, pphh_count = 0;
  for(int sym = 0; sym < phis_.number_irreps(); sym++) {
    pphh_count += 
      occs_[sym].size() * (occs_[sym].size() + 1) / 2
      + virs_[sym].size() * (virs_[sym].size() + 1) / 2;
    
    ph_count += occs_[sym].size() * virs_[sym].size() ;
  }
  // Compose the matrices A11 and A21 eq.(B.5) 
  // Also separate B and Sigma into B1 and B2 (and assign the latter to Sig2 directly)
  // and Sig1 and Sig2, respectively.
  Blas_matrix A11(ph_count, ph_count), A21(pphh_count, ph_count), 
    B1(ph_count, 1),
    Sig1(ph_count, 1), Sig2(pphh_count, 1);
  
  Sig1 = 0.;
  
  // Initialize the p-h pair elements
  unsigned row_count = 0;
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < virs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
	unsigned i = virs_[sym][i_];
	unsigned j = occs_[sym][j_];
	
	B1(row_count++) = b(i,j);
      }
  // and the h-h pair  and  p-p pair elements of B
  row_count = 0;
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ <= i_; j_++) {
	unsigned i = occs_[sym][i_];
	unsigned j = occs_[sym][j_];
	
	Sig2(row_count++) = b(i,j);
      }
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < virs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ <= i_; j_++) {
	unsigned i = virs_[sym][i_];
	unsigned j = virs_[sym][j_];
	
	Sig2(row_count++) = b(i,j);
      }
  
  // Initialize the p-h/p-h pair elements of A11
  row_count = 0;
  for(unsigned sym_pq = 0; sym_pq < phis_.number_irreps(); sym_pq++) 
    for(unsigned q_ = 0; q_ < virs_[sym_pq].size(); q_++)
      for(unsigned p_ = 0; p_ < occs_[sym_pq].size(); p_++) {
	unsigned q = virs_[sym_pq][q_];
	unsigned p = occs_[sym_pq][p_];
	
	unsigned col_count = 0;
	for(unsigned sym_lk = 0; sym_lk < phis_.number_irreps(); sym_lk++) 
	  for(unsigned l_ = 0; l_ < virs_[sym_lk].size(); l_++)
	    for(unsigned k_ = 0; k_ < occs_[sym_lk].size(); k_++) {
	      unsigned l = virs_[sym_lk][l_];
	      unsigned k = occs_[sym_lk][k_];
	      
	      A11(row_count, col_count) =
		(2.* (V1212(p,k,q,l)+V1212(p,l,q,k)) - (V1212(p,k,l,q)+V1212(p,l,k,q)))
		/(phis_.energy(k) - phis_.energy(l)); 
	      
	      col_count++;
	    }
	row_count++;
      }
  // and the h-h/p-h pair elements of A21
  row_count = 0;
  for(unsigned sym_pq = 0; sym_pq < phis_.number_irreps(); sym_pq++) 
    for(unsigned q_ = 0; q_ < occs_[sym_pq].size(); q_++)
      for(unsigned p_ = 0; p_ <= q_; p_++) {
	unsigned q = occs_[sym_pq][q_];
	unsigned p = occs_[sym_pq][p_];
	
	unsigned col_count = 0;
	for(unsigned sym_lk = 0; sym_lk < phis_.number_irreps(); sym_lk++) 
	  for(unsigned l_ = 0; l_ < virs_[sym_lk].size(); l_++)
	    for(unsigned k_ = 0; k_ < occs_[sym_lk].size(); k_++) {
	      unsigned l = virs_[sym_lk][l_];
	      unsigned k = occs_[sym_lk][k_];
	      
	      A21(row_count, col_count) =
		(2.* (V1212(p,k,q,l)+V1212(p,l,q,k)) - (V1212(p,k,l,q)+V1212(p,l,k,q)))
		/(phis_.energy(k) - phis_.energy(l));
	      col_count++;
	    }
	row_count++;
      }
  // and the p-p/p-h pair elements of A21
  for(unsigned sym_pq = 0; sym_pq < phis_.number_irreps(); sym_pq++) 
    for(unsigned q_ = 0; q_ < virs_[sym_pq].size(); q_++)
      for(unsigned p_ = 0; p_ <= q_; p_++) {
	unsigned q = virs_[sym_pq][q_];
	unsigned p = virs_[sym_pq][p_];
	
	unsigned col_count = 0;
	for(unsigned sym_lk = 0; sym_lk < phis_.number_irreps(); sym_lk++) 
	  for(unsigned l_ = 0; l_ < virs_[sym_lk].size(); l_++)
	    for(unsigned k_ = 0; k_ < occs_[sym_lk].size(); k_++) {
	      unsigned l = virs_[sym_lk][l_];
	      unsigned k = occs_[sym_lk][k_];
	      
	      A21(row_count, col_count) = 
		(2.* (V1212(p,k,q,l)+V1212(p,l,q,k)) - (V1212(p,k,l,q)+V1212(p,l,k,q)))
		/(phis_.energy(k) - phis_.energy(l));
	      col_count++;
	    }
	row_count++;
      }
  
  //Compute Sigma(inf) eq. (B.6)
  A11.add_diag(-1.);
  A11.inverse();
  
  Sig1.dgemv('N',A11, B1);
  Sig1 *= -1.;
  Sig2.dgemv('N', A21, Sig1);
  
  // Save the results
  // p-h pairs
  row_count = 0;
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < virs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
	unsigned i = virs_[sym][i_];
	unsigned j = occs_[sym][j_];
	
	b(j,i) = b(i,j) = Sig1(row_count++);
      }
  // h-h pairs
  row_count = 0;
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ <= i_; j_++) {
	unsigned i = occs_[sym][i_];
	unsigned j = occs_[sym][j_];

	b(j,i) = b(i,j) = Sig2(row_count++);
      }
  // p-p pairs
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned i_ = 0; i_ < virs_[sym].size(); i_++)
      for(unsigned j_ = 0; j_ <= i_; j_++) {
	unsigned i = virs_[sym][i_];
	unsigned j = virs_[sym][j_];
	
	b(j,i) = b(i,j) = Sig2(row_count++);
      }
}
 



/////////////////////////// integral driven code addtions




inline static void swap(unsigned& a, unsigned& b) {unsigned c=a; a=b;b=c;}


inline double idSelf_energy::eq20term_1(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{
  double sum1=0.;
  
  unsigned sym = phis_.irrep(k);
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum1 += V1212(a,a,k,l) * V1212(a,a,k1,l)
      /(
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	);
  }  
  return sum1;

}


inline double idSelf_energy::eq20term_2(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.,sum2=0.;
  unsigned sym = phis_.irrep(k);
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum1 += (V1212(a,d,k,l)+V1212(a,d,l,k))* V1212(a,a,k1,l)
      /(
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	);

    sum2 += V1212(a,a,k,l)*(V1212(a,d,k1,l)+ V1212(a,d,l,k1))
      /(
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);
  }  
  return sum1+sum2;

}

inline double idSelf_energy::eq20term_3(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.;
  unsigned sym = 
    phis_.irrep_product(phis_.irrep(k),
			phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum1 += (V1212(a,d,l,k)*(2.*V1212(a,d,l,k1)-V1212(a,d,k1,l))+V1212(a,d,k,l)*(2.*V1212(a,d,k1,l)-V1212(a,d,l,k1)))
      /(
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);
  }  
  return sum1;

}



inline double idSelf_energy::eq20term_4(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.,sum2=0.,sum3=0.;
    unsigned sym = 
      phis_.irrep(k);
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];


    sum1 += V1212(d,d,k,l) * V1212(a,a,k1,l)
      /(
	(phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	);


    sum3 += V1212(a,a,k,l) * V1212(d,d,k1,l)
      /(
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(d)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);

  }
  
  sym =   phis_.irrep_product(phis_.irrep(k),
			      phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];
    
    sum2 += (V1212(a,d,l,k)*(2.*V1212(a,d,k1,l)-V1212(a,d,l,k1))+V1212(a,d,k,l)*(2.*V1212(a,d,l,k1)-V1212(a,d,k1,l)))
      /(
	(phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);
    
  }  
  
  return sum1+sum2+sum3;
  
}




inline double idSelf_energy::eq20term_6(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.,sum2=0.;
  
  unsigned sym =   phis_.irrep_product(phis_.irrep(k),
			      phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];
    

    sum1 += (V1212(a,d,k,l)*(2.*V1212(a,c,k1,l)-V1212(a,c,l,k1))+V1212(a,d,l,k)*(2.*V1212(a,c,l,k1)-V1212(a,c,k1,l)))
      /(
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	);


    sum2 += (V1212(a,c,k,l)*(2.*V1212(a,d,k1,l)-V1212(a,d,l,k1))+V1212(a,c,l,k)*(2.*V1212(a,d,l,k1)-V1212(a,d,k1,l)))
      /(
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);


  }  
  return sum1+sum2;

}

inline double idSelf_energy::eq20term_7(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.,sum2=0.,sum3=0.,sum4=0.;
  unsigned sym = 
    phis_.irrep(k);
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum1 += (V1212(b,d,k,l)+V1212(b,d,l,k))* V1212(a,a,k1,l)
      /(
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	);

    sum4 +=  V1212(a,a,k,l)*(V1212(b,d,k1,l)+V1212(b,d,l,k1))
      /(
	(phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);
  }

  sym =   phis_.irrep_product(phis_.irrep(k),
			      phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum2 += (V1212(a,b,k,l)*(2.*V1212(a,d,l,k1)-V1212(a,d,k1,l))+V1212(a,b,l,k)*(2.*V1212(a,d,k1,l)-V1212(a,d,l,k1)))
      /(
	(phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);

    sum3 += (V1212(a,d,k,l)*(2.*V1212(a,b,l,k1)-V1212(a,b,k1,l))+V1212(a,d,l,k)*(2.*V1212(a,b,k1,l)-V1212(a,b,l,k1)))
      /(
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(a)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	);



  }  
  return sum1+sum2+sum3+sum4;

}




inline double idSelf_energy::eq20term_12(unsigned  a, unsigned b, unsigned c, unsigned d, unsigned k, unsigned k1)
{

  double sum1=0.,sum2=0.,sum3=0.,sum4=0.;

  unsigned sym =   phis_.irrep_product(phis_.irrep(k),
			      phis_.irrep_product(phis_.irrep(a),phis_.irrep(c)));
  
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];

    sum1 += (V1212(b,d,k,l)*(2.*V1212(a,c,k1,l)-V1212(a,c,l,k1))+V1212(d,b,k,l)*(2.*V1212(a,c,l,k1)-V1212(a,c,k1,l)))
      /(
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	);


    sum4 += (V1212(a,c,k,l)*(2.*V1212(b,d,k1,l)-V1212(b,d,l,k1))+V1212(a,c,l,k)*(2.*V1212(b,d,l,k1)-V1212(b,d,k1,l)))
      /(
	(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);
  }

  sym =   phis_.irrep_product(phis_.irrep(k),
			      phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  
  for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
    unsigned l = occs_[sym][l_];


    sum2 += (V1212(b,c,k,l)*(2.*V1212(a,d,k1,l)-V1212(a,d,l,k1))+V1212(b,c,l,k)*(2.*V1212(a,d,l,k1)-V1212(a,d,k1,l)))
      /(
	(phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	);


    sum3 += (V1212(a,d,k,l)*(2.*V1212(b,c,k1,l)-V1212(b,c,l,k1))+V1212(a,d,l,k)*(2.*V1212(b,c,l,k1)-V1212(b,c,k1,l)))
      /(
	(phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(c)-phis_.energy(k1)-phis_.energy(l))*
	(phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	);





  }  
  return sum1+sum2+sum3+sum4;

}



inline void idSelf_energy::eq20(Blas_matrix& rho,  Term func, double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{

  for (unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned k1_ = 0; k1_ < occs_[sym].size() ; k1_++) {
	unsigned k = occs_[sym][k_];
	unsigned k1 = occs_[sym][k1_];
	double fA = 0.;
	
	fA = eri*(this->*func)(a, b, c, d, k, k1);
	
	
	double fkk = fA;
	
	rho(k,k1) += fkk; 
	rho(k1,k) += fkk;
      }
}


double idSelf_energy::eq_c12term_1(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  
  double sum1 = 0.; 
  
  unsigned sym =   phis_.irrep(p);
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];
    
    sum1 += V1212(p,i,a,a) * V1212(q,i,a,a)
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(a))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(a))
	 );
  }
  
  return sum1;

}



double idSelf_energy::eq_c12term_2(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0.; 
  
  unsigned sym =   phis_.irrep(p);
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];
    
    sum1 += (V1212(p,i,a,a) * (V1212(q,i,a,d)+V1212(q,i,d,a))
	     +(V1212(p,i,a,d)+V1212(p,i,d,a))*V1212(q,i,a,a))
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(d)-phis_.energy(a))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(a))
	 );
  }
  
  return sum1;

}



double idSelf_energy::eq_c12term_3(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0.; 
  
  
  
  unsigned sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];
    
    sum1 += ((2.*V1212(p,i,a,d)-V1212(p,i,d,a))*V1212(q,i,a,d)
	     +(2.*V1212(p,i,d,a)-V1212(p,i,a,d))*V1212(q,i,d,a))
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(d))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(d))
	 );
  }
  
  return sum1;

}


double idSelf_energy::eq_c12term_4(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0., sum2=0.; 
  
  
  unsigned sym =   phis_.irrep(p);
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];
  

    
    sum1 += (V1212(p,i,a,a)*V1212(q,i,d,d)+V1212(p,i,d,d)*V1212(q,i,a,a))
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(a))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(d)-phis_.energy(d))
	 );
  }

  sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];

    sum2 += ((2.*V1212(p,i,a,d)-V1212(p,i,d,a))*V1212(q,i,d,a)
	     +(2.*V1212(p,i,d,a)-V1212(p,i,a,d))*V1212(q,i,a,d))
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(d))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(d)-phis_.energy(a))
	 );
  }
  
  return sum1+sum2;

}



double idSelf_energy::eq_c12term_6(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0.; 
  
  unsigned sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];

    sum1 += (
	     (2.*V1212(p,i,a,c)-V1212(p,i,c,a))*V1212(q,i,a,d)+
	     (2.*V1212(p,i,a,d)-V1212(p,i,d,a))*V1212(q,i,a,c)+
	     (2.*V1212(p,i,c,a)-V1212(p,i,a,c))*V1212(q,i,d,a)+
	     (2.*V1212(p,i,d,a)-V1212(p,i,a,d))*V1212(q,i,c,a)
	     )
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(c))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(d))
	 );
  }
  
  return sum1;

}


double idSelf_energy::eq_c12term_7(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0., sum2 =0.; 
  
  unsigned sym =   phis_.irrep(p);
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];

    

    sum1 += (
	     V1212(p,i,a,a)*(V1212(q,i,b,d)+V1212(q,i,d,b))+
	     (V1212(p,i,b,d)+V1212(p,i,d,b))*V1212(q,i,a,a)
	     )
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(a))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(b)-phis_.energy(d))
	 );
  }

  sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];

    sum2 += (
	     (2.*V1212(p,i,a,d)-V1212(p,i,d,a))*V1212(q,i,b,a)+
	     (2.*V1212(p,i,b,a)-V1212(p,i,a,b))*V1212(q,i,a,d)+
	     (2.*V1212(p,i,a,b)-V1212(p,i,b,a))*V1212(q,i,d,a)+
	     (2.*V1212(p,i,d,a)-V1212(p,i,a,d))*V1212(q,i,a,b)
	     )
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(d)-phis_.energy(a))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(b))
	 );

  }
  
  return sum1+sum2;

}



double idSelf_energy::eq_c12term_12(unsigned a, unsigned b, unsigned c, unsigned d, unsigned p, unsigned q)
{

  double e_q = phis_.energy(q);
  
  double sum1 = 0., sum2 =0.; 
  
  

  unsigned sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(c)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];
  

    sum1 += (
	     (2.*V1212(p,i,a,c)-V1212(p,i,c,a))*V1212(q,i,b,d)+
	     (2.*V1212(p,i,d,b)-V1212(p,i,b,d))*V1212(q,i,c,a)+
	     (2.*V1212(p,i,c,a)-V1212(p,i,a,c))*V1212(q,i,d,b)+
	     (2.*V1212(p,i,b,d)-V1212(p,i,d,b))*V1212(q,i,a,c)
	     )
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(c))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(b)-phis_.energy(d))
	 );

  }

  sym =   phis_.irrep_product(phis_.irrep(p),phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
  for(unsigned i_ = 0; i_ < occs_[sym].size(); i_++) {
    unsigned i = occs_[sym][i_];

    sum2 += (
	     (2.*V1212(p,i,a,d)-V1212(p,i,d,a))*V1212(q,i,b,c)+
	     (2.*V1212(p,i,b,c)-V1212(p,i,c,b))*V1212(q,i,a,d)+
	     (2.*V1212(p,i,c,b)-V1212(p,i,b,c))*V1212(q,i,d,a)+
	     (2.*V1212(p,i,d,a)-V1212(p,i,a,d))*V1212(q,i,c,b)
	     )
      / (
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(a)-phis_.energy(d))*
	 (phis_.energy(q)+phis_.energy(i)-phis_.energy(b)-phis_.energy(c))
	 );

    


  }
  
  return sum1+sum2;

}




void idSelf_energy::eq_c12(Blas_matrix& M_ak, Term func, double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{
  unsigned num_occ = phis_.number_occupied();
  
  // eq.C12
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned p_ = 0; p_ < virs_[sym].size(); p_++) 
      for(unsigned q_ = 0; q_ < occs_[sym].size(); q_++) {
	unsigned p = virs_[sym][p_];
	unsigned q = occs_[sym][q_];
	
	double C1pq = 0.; 

	C1pq = eri * (this->*func)(a,b,c,d,p,q);
	
	M_ak(p-num_occ,q) += C1pq;

      }
}


void idSelf_energy::eq_c16term_1(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{

  

  
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);
  if (!sym_aq) return;

  double sum1 = 0.;

  
  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym = phis_.irrep(i);
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];
      
      
      sum1 += V1212(i,j,a,a)*V1212(q,a,i,j)
	
	/ (
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(a))
	   );
    }
  }
  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * sum1;
}




void idSelf_energy::eq_c16term_2(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{


  
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);
  if (!sym_aq) return;

  double sum1 = 0., sum2=0.;


  
  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym = phis_.irrep(i);
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];
      
      sum1 += ((V1212(i,j,a,d)+V1212(i,j,d,a))*V1212(q,a,i,j)+V1212(i,j,a,a)*V1212(q,d,i,j))
	/ (
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	   );

      sum2 += V1212(i,j,a,a)*V1212(q,a,i,j)
	/ (
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(a))
	   );
    }
  }
  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * sum1;
  mat(d-num_occ,q) += vi * sum2;
}


void idSelf_energy::eq_c16term_3(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{

  

  
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);
  bool sym_dq = phis_.irrep(d) == phis_.irrep(q);

  if (!sym_aq && !sym_dq) return;

  double sum1 = 0., sum2 = 0.;
    
    
  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym =  phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];
      
      if (sym_aq)
	sum1 += (2.*V1212(i,j,a,d)-V1212(i,j,d,a))*V1212(q,d,i,j)
	/
	(
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	 );
      
      if (sym_dq)
	sum2 += (2.*V1212(i,j,d,a)-V1212(i,j,a,d))*V1212(q,a,i,j)
	  /
	  (
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	   );
      
      
    }
  }  
  mat(a-phis_.number_occupied(),q) += vi * sum1;
  mat(d-phis_.number_occupied(),q) += vi * sum2;

}



void idSelf_energy::eq_c16term_4(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{

  

  bool sym_dq = phis_.irrep(d) == phis_.irrep(q);
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);

  if (!sym_aq && !sym_dq) return;

  double sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0.;

  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym =  phis_.irrep(i);
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];


      if (sym_aq) 
	sum1 += V1212(i,j,d,d)*V1212(q,a,i,j) 
	  /
	  (
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(d)-phis_.energy(d))
	   );


      if (sym_dq) 
	sum4 += V1212(i,j,a,a)*V1212(q,d,i,j) 
	  /
	  (
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(d)-phis_.energy(d))
	   );



    }
      
    sym =  phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];

      if (sym_aq)
	sum2 += (2.*V1212(i,j,d,a)-V1212(i,j,a,d))*V1212(q,d,i,j) 
	  /
	(
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	 );

      if (sym_dq)
	sum3 += (2.*V1212(i,j,a,d)-V1212(i,j,d,a))*V1212(q,a,i,j)
	  /
	  (
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	   (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	   );
      
    }
  }

  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * (sum1+sum2);
  mat(d-num_occ,q) += vi * (sum3+sum4);
}



void idSelf_energy::eq_c16term_6(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{
  
  

  bool sym_dq = phis_.irrep(d) == phis_.irrep(q);
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);

  if (!sym_aq  && !sym_dq ) return;

  double sum1 = 0., sum2 = 0., sum3 = 0.;



  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym =   phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];


      
      double denum1 =  1 /
	(
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(d))
	 );

      if (sym_aq) {
	sum1 += ((2.*V1212(i,j,a,d)-V1212(i,j,d,a))*V1212(q,c,i,j)
		 +(2.*V1212(i,j,a,c)-V1212(i,j,c,a))*V1212(q,d,i,j))
	  * denum1;
      }
      
      if (sym_dq){
	sum2 += (2.*V1212(i,j,c,a)-V1212(i,j,a,c))*V1212(q,a,i,j) * denum1;
	sum3 += (2.*V1212(i,j,d,a)-V1212(i,j,a,d))*V1212(q,a,i,j) * denum1;
      }

    }
  }
  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * sum1;
  mat(d-num_occ,q) += vi * sum2;
  mat(c-num_occ,q) += vi * sum3;

}



void idSelf_energy::eq_c16term_7(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{


  bool sym_dq = phis_.irrep(d) == phis_.irrep(q);
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);

  if (!sym_aq  && !sym_dq ) return;

  double sum1 = 0., sum2 = 0., sum3 = 0.;
  double sum4 = 0., sum5 = 0., sum6 = 0.;
  

      

  for(int i = 0; i < phis_.number_occupied();i++) {
    unsigned sym =   phis_.irrep(i);
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];



      double denum2 =  1 /
	(
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(d)-phis_.energy(b))
	 );



      if (sym_aq) {
	sum1 += (V1212(i,j,b,d)+V1212(i,j,d,b))*V1212(q,a,i,j)
	  * denum2;
      }


      if (sym_dq){

	sum4 += V1212(i,j,a,a)*V1212(q,b,i,j) * denum2;

	sum6 += V1212(i,j,a,a)*V1212(q,d,i,j) * denum2;
      }
    }

    sym =   phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];


      double denum1 =  1 /
	(
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(a))
	 );

      if (sym_aq) {

	sum2 += ((2.*V1212(i,j,b,a)-V1212(i,j,a,b))*V1212(q,d,i,j)
		 +(2.*V1212(i,j,d,a)-V1212(i,j,a,d))*V1212(q,b,i,j))
	  * denum1;
      }
      
      if (sym_dq){
	sum3 += (2.*V1212(i,j,a,b)-V1212(i,j,b,a))*V1212(q,a,i,j) * denum1;

	sum5 += (2.*V1212(i,j,a,d)-V1212(i,j,d,a))*V1212(q,a,i,j) * denum1;
      }

    }
  }
  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * (sum1+sum2);
  mat(d-num_occ,q) += vi * (sum3+sum4);
  mat(b-num_occ,q) += vi * (sum5+sum6);

}


void idSelf_energy::eq_c16term_12(unsigned a, unsigned b, unsigned c, unsigned d, unsigned q, Blas_matrix& mat, double vi)
{


  bool sym_dq = phis_.irrep(d) == phis_.irrep(q);
  bool sym_aq = phis_.irrep(a) == phis_.irrep(q);
  bool sym_bq = phis_.irrep(b) == phis_.irrep(q);
  bool sym_cq = phis_.irrep(c) == phis_.irrep(q);


  if (!sym_aq  && !sym_dq && !sym_bq && !sym_cq ) return;

  double sum1 = 0., sum2 = 0., sum3 = 0.;
  double sum4 = 0., sum5 = 0., sum6 = 0.;
  double sum7 = 0., sum8 = 0.;
  

  for(int i = 0; i < phis_.number_occupied();i++) {

    unsigned sym =   phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(c)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];
      

      
      double denum1 =  1 /
	(
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(d))
	 );

      if (sym_aq) {
	sum1 += (2.*V1212(i,j,b,d)-V1212(i,j,d,b))*V1212(q,c,i,j)
	  * denum1;
      }
      if (sym_bq){
	sum4 += (2.*V1212(i,j,a,c)-V1212(i,j,c,a))*V1212(q,d,i,j) * denum1;
      }
      if (sym_cq){
	sum5 += (2.*V1212(i,j,d,b)-V1212(i,j,b,d))*V1212(q,a,i,j) * denum1;
      }
      if (sym_dq){
	sum8 += (2.*V1212(i,j,c,a)-V1212(i,j,a,c))*V1212(q,b,i,j) * denum1;
      }

    }

    sym =   phis_.irrep_product(phis_.irrep(i), phis_.irrep_product(phis_.irrep(a),phis_.irrep(d)));
    
    for(unsigned j_ = 0; j_ < occs_[sym].size(); j_++) {
      unsigned j = occs_[sym][j_];
      


      double denum2 =  1 /
	(
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(i)-phis_.energy(j))*
	 (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c))
	 );




      if (sym_aq) {
	sum2 += (2.*V1212(i,j,b,c)-V1212(i,j,c,b))*V1212(q,d,i,j)
	  * denum2;
      }
      
      if (sym_bq){
	sum3 += (2.*V1212(i,j,a,d)-V1212(i,j,d,a))*V1212(q,c,i,j) * denum2;
      }
      if (sym_cq){
	sum6 += (2.*V1212(i,j,d,a)-V1212(i,j,a,d))*V1212(q,b,i,j) * denum2;
      }

      if (sym_dq){
	sum7 += (2.*V1212(i,j,c,b)-V1212(i,j,b,c))*V1212(q,a,i,j) * denum2;
      }


    }
  }

  unsigned num_occ = phis_.number_occupied();

  mat(a-num_occ,q) += vi * (sum1+sum2);
  mat(b-num_occ,q) += vi * (sum3+sum4);
  mat(c-num_occ,q) += vi * (sum5+sum6);
  mat(d-num_occ,q) += vi * (sum7+sum8);

}


void idSelf_energy::eq_c16(Blas_matrix& M_ak, Term1 func, double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{
  unsigned sym = phis_.irrep(a);
  
  // eq.C16
  for(int q = 0; q < phis_.number_occupied();q++)
    (this->*func)(a,b,c,d,q, M_ak, eri);

  

}



void idSelf_energy::eq_30term_1(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{

  
  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  if (!sym_av) return;

  double sum1 = 0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {

    sum1 += V1212(a,v,k,k) * V1212(a,a,k,k)
      /
      (
       (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
       (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
       (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
       );



    unsigned sym =   phis_.irrep_product(phis_.irrep(v), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;
    
      sum1 += (V1212(a,v,k,l)+V1212(a,v,l,k)) * V1212(a,a,k,l)
	/
	(
	 (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	 );
      
      
    }  
  }
  

  mat(a,v) += vi * sum1;
  mat(v,a) += vi * sum1;
}



void idSelf_energy::eq_30term_2(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{

  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  if (!sym_av) return;

  
  double sum1 = 0., sum2 = 0., sum3 = 0.;

  for(int k = 0; k < phis_.number_occupied();k++) {

      sum1 += V1212(d,v,k,k) * V1212(a,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );
      
      sum2 += 2.*V1212(a,v,k,k) * V1212(a,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	 );

      sum3 += V1212(a,v,k,k) * V1212(a,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	 );

      unsigned sym =   phis_.irrep_product(phis_.irrep(v), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
      
      for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
	unsigned l = occs_[sym][l_];
	if (l >= k ) break;
	
	sum1 += (V1212(a,v,k,l)+V1212(a,v,l,k)) * (V1212(a,d,k,l)+V1212(a,d,l,k))
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
	
	sum2 += (V1212(d,v,k,l)+V1212(d,v,l,k)) * V1212(a,a,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );

	sum3 += (V1212(a,v,k,l)+V1212(a,v,l,k)) * V1212(a,a,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );

    }  
  }
  

  mat(a,v) += vi * (sum1+sum2);
  mat(v,a) += vi * (sum1+sum2);

  mat(d,v) += vi * sum3;
  mat(v,d) += vi * sum3;
}




void idSelf_energy::eq_30term_3(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{


  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  bool sym_dv = phis_.irrep(d) == phis_.irrep(v);
  bool sym_da = phis_.irrep(d) == phis_.irrep(a);

  if (!sym_av && !sym_dv) return;

  double sum1 = 0., sum2 = 0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {

    if (sym_da) {

      sum1 += V1212(d,v,k,k) * V1212(d,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );
      
      sum2 += V1212(a,v,k,k) * V1212(a,d,k,k)
	/
	(
       (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
       (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
       (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
       );
    }
    
    
    unsigned sym =   phis_.irrep_product(phis_.irrep(d), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;
      
      if (sym_av) {
      sum1 +=  (
		V1212(d,v,l,k)*(2.*V1212(d,a,l,k)-V1212(d,a,k,l)) +
		V1212(d,v,k,l)*(2.*V1212(d,a,k,l)-V1212(d,a,l,k))
		)
	/
	(
	 (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	 );
      }
      if (sym_dv) {

	sum2 += (
		 V1212(a,v,l,k)*(2.*V1212(a,d,l,k)-V1212(a,d,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(a,d,k,l)-V1212(a,d,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }
    
    }  
  }
  

    mat(a,v) += vi * sum1;
    mat(v,a) += vi * sum1;


    mat(d,v) += vi * sum2;
    mat(v,d) += vi * sum2;

}



void idSelf_energy::eq_30term_4(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{


  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  bool sym_dv = phis_.irrep(d) == phis_.irrep(v);
  bool sym_da = phis_.irrep(d) == phis_.irrep(a);

  if (!sym_av && !sym_dv) return;

  double sum1 = 0., sum2 = 0.,sum3=0.,sum4=0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {
    if (sym_av) {
      sum1 += V1212(a,v,k,k) * V1212(d,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	 );
    if (sym_da) 
      sum2 += V1212(d,v,k,k) * V1212(a,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );

    }

    if (sym_dv) {

      if (sym_da) 
	sum3 += V1212(a,v,k,k) * V1212(d,a,k,k)
	  /
	  (
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
 	 );
      
      sum4 += V1212(d,v,k,k) * V1212(a,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
 	 );

      
    }


    
    unsigned sym =   phis_.irrep_product(phis_.irrep(d), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;

      if (sym_av) {

	sum2 +=  (
		  V1212(d,v,l,k)*(2.*V1212(a,d,l,k)-V1212(a,d,k,l)) +
		  V1212(d,v,k,l)*(2.*V1212(a,d,k,l)-V1212(a,d,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );

	
      }
      
      if (sym_dv) {
	sum3 += (
		 V1212(a,v,l,k)*(2.*V1212(d,a,l,k)-V1212(d,a,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(d,a,k,l)-V1212(d,a,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }
    }

    sym =   phis_.irrep(k);
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;

      if (sym_dv) {

	sum4 +=  (V1212(d,v,k,l)+V1212(d,v,l,k))*V1212(a,a,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );


      }
      if (sym_av) {
	sum1 +=  (V1212(a,v,k,l)+V1212(a,v,l,k))*V1212(d,d,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }

    }  
  }
  
  if (sym_av) {
    mat(a,v) += vi * (sum1+sum2);
    mat(v,a) += vi * (sum1+sum2);
  }


  if (sym_dv) {
    mat(d,v) += vi * (sum3+sum4);
    mat(v,d) += vi * (sum3+sum4);
  }
}



void idSelf_energy::eq_30term_6(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{


  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  bool sym_dv = phis_.irrep(d) == phis_.irrep(v);
  bool sym_da = phis_.irrep(d) == phis_.irrep(a);

  if (!sym_av && !sym_dv) return;

  double sum1 = 0., sum2 = 0.,sum3=0.,sum4=0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {
    if (sym_da) {
      sum1 += V1212(c,v,k,k) * V1212(d,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(c)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))
	 );

      sum2 += V1212(d,v,k,k) * V1212(c,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );


      sum3 += V1212(a,v,k,k) * V1212(a,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
 	 );

      sum4 += V1212(a,v,k,k) * V1212(a,c,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
 	 );
    }


    
    unsigned sym =   phis_.irrep_product(phis_.irrep(d), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;

      if (sym_av) {
	sum1 +=  (
		  V1212(c,v,l,k)*(2.*V1212(d,a,l,k)-V1212(d,a,k,l)) +
		  V1212(c,v,k,l)*(2.*V1212(d,a,k,l)-V1212(d,a,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(c)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	   );


	sum2 +=  (
		  V1212(d,v,l,k)*(2.*V1212(c,a,l,k)-V1212(c,a,k,l)) +
		  V1212(d,v,k,l)*(2.*V1212(c,a,k,l)-V1212(c,a,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );

	
      }
      
      if (sym_dv) {
	sum3 += (
		 V1212(a,v,l,k)*(2.*V1212(a,d,l,k)-V1212(a,d,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(a,d,k,l)-V1212(a,d,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );

	sum4 += (
		 V1212(a,v,l,k)*(2.*V1212(a,c,l,k)-V1212(a,c,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(a,c,k,l)-V1212(a,c,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }


    }  
  }
  
  if (sym_av) {
    mat(a,v) += vi * (sum1+sum2);
    mat(v,a) += vi * (sum1+sum2);
  }

  if (sym_dv) {
    mat(c,v) += vi * sum3;
    mat(v,c) += vi * sum3;

    mat(d,v) += vi * sum4;
    mat(v,d) += vi * sum4;
  }
}



void idSelf_energy::eq_30term_7(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{


  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  bool sym_dv = phis_.irrep(d) == phis_.irrep(v);
  bool sym_da = phis_.irrep(d) == phis_.irrep(a);

  if (!sym_av && !sym_dv) return;


  double sum1 = 0., sum2 = 0.,sum3=0.,sum4=0.,sum5=0.,sum6=0.,sum7=0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {
    if (sym_av) {
      sum1 += V1212(a,v,k,k) * (V1212(b,d,k,k)+V1212(b,d,k,k))
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	 );
      if (sym_da) {
	sum2 += V1212(b,v,k,k) * V1212(a,d,k,k)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(k)-phis_.energy(b)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	   (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))
	   );
      
	sum3 += V1212(d,v,k,k) * V1212(a,b,k,k)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	   );
      }
    }

    if (sym_dv) {
      if (sym_da) {
	
	sum4 += V1212(a,v,k,k) * V1212(d,a,k,k)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	   (phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	   );

	sum6 += V1212(a,v,k,k) * V1212(b,a,k,k)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
	   );
      
      }

      sum5 += V1212(d,v,k,k) * V1212(a,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );


      sum7 += V1212(b,v,k,k) * V1212(a,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(b)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))
 	 );
      
    }
    

    
    unsigned sym =   phis_.irrep(k);
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;

      if (sym_av) {
	sum1 +=  ( V1212(a,v,k,l)+V1212(a,v,l,k))*(V1212(d,b,k,l)+V1212(d,b,l,k))
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }


      if (sym_dv) {

	sum5 += (V1212(d,v,k,l)+V1212(d,v,l,k))*V1212(a,a,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );


	sum7 += (V1212(b,v,k,l)+V1212(b,v,l,k))*V1212(a,a,k,l)
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(b)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))
	   );
	

      }

    }
    

    sym =   phis_.irrep_product(phis_.irrep(d), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;

      
      if (sym_av) {
	
	sum2 +=  (
		  V1212(b,v,l,k)*(2.*V1212(a,d,l,k)-V1212(a,d,k,l)) +
		  V1212(b,v,k,l)*(2.*V1212(a,d,k,l)-V1212(a,d,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(b)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))
	   );

	sum3 +=  (
		  V1212(d,v,l,k)*(2.*V1212(a,b,l,k)-V1212(a,b,k,l)) +
		  V1212(d,v,k,l)*(2.*V1212(a,b,k,l)-V1212(a,b,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );
	
      }
      
      if (sym_dv) {
	sum4 += (
		 V1212(a,v,l,k)*(2.*V1212(d,a,l,k)-V1212(d,a,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(d,a,k,l)-V1212(d,a,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );


	sum6 += (
		 V1212(a,v,l,k)*(2.*V1212(b,a,l,k)-V1212(b,a,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(b,a,k,l)-V1212(b,a,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(b)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      
      }


    }  
  }
  
  if (sym_av) {
    mat(a,v) += vi * (sum1+sum2+sum3);
    mat(v,a) += vi * (sum1+sum2+sum3);
  }

  if (sym_dv) {
    mat(b,v) += vi * (sum4+sum5);
    mat(v,b) += vi * (sum4+sum5);

    mat(d,v) += vi * (sum6+sum7);
    mat(v,d) += vi * (sum6+sum7);
  }
}


void idSelf_energy::eq_30term_12(unsigned a, unsigned b, unsigned c, unsigned d, unsigned v, Blas_matrix& mat, double vi)
{


  bool sym_av = phis_.irrep(a) == phis_.irrep(v);
  bool sym_dv = phis_.irrep(d) == phis_.irrep(v);
  bool sym_bv = phis_.irrep(b) == phis_.irrep(v);
  bool sym_cv = phis_.irrep(c) == phis_.irrep(v);

  
  if (!sym_av && !sym_dv  && !sym_cv  && !sym_bv) return;

  bool sym_ba = phis_.irrep(a) == phis_.irrep(b);
  bool sym_da = phis_.irrep(a) == phis_.irrep(d);
  bool sym_ca = phis_.irrep(a) == phis_.irrep(c);



  double sum1 = 0., sum2 = 0.,sum3=0.,sum4=0.,sum5=0.,sum6=0.,sum7=0.,sum8=0.;
  
  for(int k = 0; k < phis_.number_occupied();k++) {
    if (sym_av) {
      if(sym_ca)
      sum1 += V1212(c,v,k,k) * V1212(d,b,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(c)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))
	 );
      if(sym_da)
      sum2 += V1212(d,v,k,k) * V1212(c,b,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(c)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
	 );
    }

    if (sym_bv) {
      if(sym_da)
      sum3 += V1212(c,v,k,k) * V1212(d,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(c)-phis_.energy(v))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))
	 );
      if(sym_ca)
      sum4 += V1212(d,v,k,k) * V1212(c,a,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(d)-phis_.energy(v))*
	 (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))
 	 );



    }

    if (sym_cv) {
      if(sym_ca)
      sum6 += V1212(a,v,k,k) * V1212(b,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
 	 );
      if(sym_da)
      sum5 += V1212(b,v,k,k) * V1212(a,d,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(b)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(c)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))
 	 );
    }
    
    if (sym_dv) {

      if(sym_da)
      sum7 += V1212(a,v,k,k) * V1212(b,c,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(a)-phis_.energy(v))*
	 (phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(k))
 	 );

      if(sym_ca)
	sum8 += V1212(b,v,k,k) * V1212(a,c,k,k)
	/
	(
	 (phis_.energy(k)+phis_.energy(k)-phis_.energy(b)-phis_.energy(v))*
	 (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(k))*
	 (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k))
 	 );
    }

    unsigned sym =   phis_.irrep_product(phis_.irrep(c), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;


      if (sym_cv) {
	

	sum5 += (
		 V1212(a,v,l,k)*(2.*V1212(b,d,l,k)-V1212(b,d,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(b,d,k,l)-V1212(b,d,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }


      if (sym_av) {


	sum1 +=  (
		  V1212(c,v,l,k)*(2.*V1212(d,b,l,k)-V1212(d,b,k,l)) +
		  V1212(c,v,k,l)*(2.*V1212(d,b,k,l)-V1212(d,b,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(c)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	   );	
      }
      if (sym_bv) {

	sum4 += (
		 V1212(d,v,l,k)*(2.*V1212(c,a,l,k)-V1212(c,a,k,l))+
		 V1212(d,v,k,l)*(2.*V1212(c,a,k,l)-V1212(c,a,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(c)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(b)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );
      }

      if (sym_dv) {

	sum8 += (
		 V1212(b,v,l,k)*(2.*V1212(a,c,l,k)-V1212(a,c,k,l))+
		 V1212(b,v,k,l)*(2.*V1212(a,c,k,l)-V1212(a,c,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(b)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))
	   );
      }
    }
    sym =   phis_.irrep_product(phis_.irrep(d), phis_.irrep_product(phis_.irrep(a),phis_.irrep(k)));
    for(unsigned l_ = 0; l_ < occs_[sym].size(); l_++) {
      unsigned l = occs_[sym][l_];
      if (l >= k ) break;


      if (sym_av) {

	sum2 +=  (
		  V1212(d,v,l,k)*(2.*V1212(c,b,l,k)-V1212(c,b,k,l)) +
		  V1212(d,v,k,l)*(2.*V1212(c,b,k,l)-V1212(c,b,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(d)-phis_.energy(v))*
	   (phis_.energy(c)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))
	   );	
      }
      
      if (sym_bv) {

	sum3 +=  (
		  V1212(c,v,l,k)*(2.*V1212(d,a,l,k)-V1212(d,a,k,l)) +
		  V1212(c,v,k,l)*(2.*V1212(d,a,k,l)-V1212(d,a,l,k))
		  )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(c)-phis_.energy(v))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))
	   );
      }
      if (sym_cv) {

	sum6 += (
		 V1212(b,v,l,k)*(2.*V1212(a,d,l,k)-V1212(a,d,k,l))+
		 V1212(b,v,k,l)*(2.*V1212(a,d,k,l)-V1212(a,d,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(b)-phis_.energy(v))*
	   (phis_.energy(a)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(c)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l))
	   );
      }

      if (sym_dv) {
	
	sum7 += (
		 V1212(a,v,l,k)*(2.*V1212(b,c,l,k)-V1212(b,c,k,l))+
		 V1212(a,v,k,l)*(2.*V1212(b,c,k,l)-V1212(b,c,l,k)) 
		 )
	  /
	  (
	   (phis_.energy(k)+phis_.energy(l)-phis_.energy(a)-phis_.energy(v))*
	   (phis_.energy(b)+phis_.energy(c)-phis_.energy(k)-phis_.energy(l))*
	   (phis_.energy(d)+phis_.energy(a)-phis_.energy(k)-phis_.energy(l))
	   );
      }
      
    }


  }  


  if (sym_av) {
    mat(a,v) += vi * (sum1+sum2);
    mat(v,a) += vi * (sum1+sum2);
  }

  if (sym_bv) {
    mat(b,v) += vi * (sum3+sum4);
    mat(v,b) += vi * (sum3+sum4);
  }

  if (sym_cv) {
    mat(c,v) += vi * (sum6+sum5);
    mat(v,c) += vi * (sum6+sum5);
  }

  if (sym_dv) {
    mat(d,v) += vi * (sum8+sum7);
    mat(v,d) += vi * (sum8+sum7);
  }
}




void idSelf_energy::eq_30(Blas_matrix& M_ak, Term1 func, double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{

  

  for(unsigned v = phis_.number_occupied(); v < phis_.number_orbitals(); v++) {
    (this->*func)(a,b,c,d,v, M_ak, eri);
  }
  

}





void idSelf_energy::contrib_4vir(Blas_matrix& rho, Blas_matrix& m_ak,double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{
  
  if ( a < b) swap(a,b);
  if ( c < d) swap(c,d);

  if ((a < c) || ((a == c) && (b < d))) {
    swap(a,c);
    swap(b,d);
  }


 
  if ((a == b) && (c == d) && (a == c)) { //1 
    //cout << "case1" << endl;
#ifndef OLD_ID_CODE    
    //cout << "eq20" << endl;
    eq20(rho, &idSelf_energy::eq20term_1, eri, a, a, a, a);
    //cout << "eq12" << endl;
    eq_c12(m_ak, &idSelf_energy::eq_c12term_1, eri, a, a, a, a);
    //cout << "eq16" << endl;
    eq_c16(m_ak, &idSelf_energy::eq_c16term_1, eri, a, a, a, a);
    //cout << "eq30" << endl;
    eq_30(rho, &idSelf_energy::eq_30term_1, eri, a, a, a, a);
#endif
#ifdef OLD_ID_CODE
    eq20term(rho, eri, eri, a, a, a, a); 
    c12c16term(m_ak, eri, eri, a, a, a, a); 
    eq30term(rho, eri, a, a, a, a); 
#endif

  } else if ((a == b) && (b == c) && (c > d)) { //2 
    //cout << "case2" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_2, eri, a, a, a, d); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_2, eri, a, a, a, d);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_2, eri, a, a, a, d);
   eq_30(rho, &idSelf_energy::eq_30term_2, eri, a, a, a, d);
#endif
#ifdef OLD_ID_CODE
    eq20term(rho, eri, eri, a, a, a, d); 
    eq20term(rho, eri, eri, a, a, d, a); 
    eq20term(rho, eri, eri, a, d, a, a); 
    eq20term(rho, eri, eri, d, a, a, a); 


    c12c16term(m_ak, eri, eri, a, a, a, d); 
    c12c16term(m_ak, eri, eri, a, a, d, a); 
    c12c16term(m_ak, eri, eri, a, d, a, a); 
    c12c16term(m_ak, eri, eri, d, a, a, a); 

    eq30term(rho, eri, a, a, a, d); 
    eq30term(rho, eri, a, a, d, a); 
    eq30term(rho, eri, a, d, a, a); 
    eq30term(rho, eri, d, a, a, a); 
#endif    

    
  } else if ((a == b) && (c == d) && (b > c)) { //3
    //cout << "case3" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_3, eri, a, d, a, d); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_3, eri, a, d, a, d);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_3, eri, a, d, a, d);
   eq_30(rho, &idSelf_energy::eq_30term_3, eri, a, d, a, d);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, 0., a, c, a, c); 
    eq20term(rho, 0., eri, c, a, a, c); 
    eq20term(rho, eri, 0., c, a, c, a); 
    eq20term(rho, 0., eri, a, c, c, a); 

   c12c16term(m_ak, eri, 0., a, c, a, c); 
   c12c16term(m_ak, 0., eri, c, a, a, c); 
   c12c16term(m_ak, eri, 0., c, a, c, a); 
   c12c16term(m_ak, 0., eri, a, c, c, a); 

    eq30term(rho, eri, a, c, a, c); 
    eq30term(rho, eri, c, a, c, a); 
#endif

  } else if ((a == c) && (b == d) && (c > b)) { //4
    //cout << "case4" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_4, eri, a, b, a, b); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_4, eri, a, b, a, b);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_4, eri, a, b, a, b);
   eq_30(rho, &idSelf_energy::eq_30term_4, eri, a, b, a, b);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, eri, b, b, a, a); 
    eq20term(rho, eri, 0.,  b, a, a, b); 
    eq20term(rho, 0., eri,  a, b, a, b); 
    eq20term(rho, eri, 0.,  a, b, b, a); 
    eq20term(rho, 0., eri,  b, a, b, a); 
    eq20term(rho, eri, eri, a, a, b, b); 

    c12c16term(m_ak, eri, eri, b, b, a, a); 
    c12c16term(m_ak, eri, 0.,  b, a, a, b); 
    c12c16term(m_ak, 0., eri,  a, b, a, b); 
    c12c16term(m_ak, eri, 0.,  a, b, b, a); 
    c12c16term(m_ak, 0., eri,  b, a, b, a); 
    c12c16term(m_ak, eri, eri, a, a, b, b); 

    eq30term(rho, eri, b, b, a, a); 
    eq30term(rho, eri, b, a, a, b); 
    eq30term(rho, eri, a, b, b, a); 
    eq30term(rho, eri, a, a, b, b); 
#endif
    
  } else if ((a > b) && (b == c) && (c == d)) { //5
    //cout << "case5" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_2, eri, d, d, d, a); 
    eq_c12(m_ak, &idSelf_energy::eq_c12term_2, eri, d, d, d, a);
    eq_c16(m_ak, &idSelf_energy::eq_c16term_2, eri, d, d, d, a);
    eq_30(rho, &idSelf_energy::eq_30term_2, eri, d, d, d, a);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, eri,  b, b, a, b); 
    eq20term(rho, eri, eri,  a, b, b, b); 
    eq20term(rho, eri, eri,  b, a, b, b); 
    eq20term(rho, eri, eri,  b, b, b, a); 


    c12c16term(m_ak, eri, eri,  b, b, a, b); 
    c12c16term(m_ak, eri, eri,  a, b, b, b); 
    c12c16term(m_ak, eri, eri,  b, a, b, b); 
    c12c16term(m_ak, eri, eri,  b, b, b, a); 

     
    eq30term(rho, eri,  b, b, a, b); 
    eq30term(rho, eri,  a, b, b, b); 
    eq30term(rho, eri,  b, a, b, b); 
    eq30term(rho, eri,  b, b, b, a); 
#endif

  } else if ((a == b) && (b > c) && (c > d)) { //6
    //cout << "case6" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_6, eri, a, a, c, d); 
    eq_c12(m_ak, &idSelf_energy::eq_c12term_6, eri, a, a, c, d);
    eq_c16(m_ak, &idSelf_energy::eq_c16term_6, eri, a, a, c, d);
    eq_30(rho, &idSelf_energy::eq_30term_6, eri, a, a, c, d);
#endif  
#ifdef OLD_ID_CODE

    eq20term(rho, eri, 0., a, d, a, c); 
    eq20term(rho, 0., eri, d, a, a, c); 
    eq20term(rho, eri, 0., a, c, a, d); 
    eq20term(rho, 0., eri, c, a, a, d); 
    eq20term(rho, eri, 0., d, a, c, a); 
    eq20term(rho, 0., eri, a, d, c, a); 
    eq20term(rho, eri, 0., c, a, d, a); 
    eq20term(rho, 0., eri, a, c, d, a); 

    c12c16term(m_ak, eri, 0., a, d, a, c); 
    c12c16term(m_ak, 0., eri, d, a, a, c); 
    c12c16term(m_ak, eri, 0., a, c, a, d); 
    c12c16term(m_ak, 0., eri, c, a, a, d); 
    c12c16term(m_ak, eri, 0., d, a, c, a); 
    c12c16term(m_ak, 0., eri, a, d, c, a); 
    c12c16term(m_ak, eri, 0., c, a, d, a); 
    c12c16term(m_ak, 0., eri, a, c, d, a); 

    
    eq30term(rho, eri, a, d, a, c); 
    eq30term(rho, eri, a, c, a, d); 
    eq30term(rho, eri, d, a, c, a); 
    eq30term(rho, eri, c, a, d, a); 
#endif

  } else if ((a == c) && (c > b) && (b > d)) { //7
    //cout << "case7" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_7, eri, a, b, a, d); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_7, eri, a, b, a, d);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_7, eri, a, b, a, d);
   eq_30(rho, &idSelf_energy::eq_30term_7, eri, a, b, a, d);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri,eri, b, d, a, a); 
    eq20term(rho, eri,eri, d, b, a, a); 
    eq20term(rho, eri,eri, a, a, b, d); 
    eq20term(rho, eri,eri, a, a, d, b); 
    eq20term(rho, eri, 0., b, a, a, d); 
    eq20term(rho, 0., eri, b, a, d, a); 
    eq20term(rho, 0., eri, a, b, a, d); 
    eq20term(rho, eri, 0., a, b, d, a); 
    eq20term(rho, eri, 0., a, d, b, a); 
    eq20term(rho, 0., eri, a, d, a, b); 
    eq20term(rho, 0., eri, d, a, b, a); 
    eq20term(rho, eri, 0., d, a, a, b); 



    c12c16term(m_ak, eri,eri, b, d, a, a); 
    c12c16term(m_ak, eri,eri, d, b, a, a); 
    c12c16term(m_ak, eri,eri, a, a, b, d); 
    c12c16term(m_ak, eri,eri, a, a, d, b); 
    c12c16term(m_ak, eri, 0., b, a, a, d); 
    c12c16term(m_ak, 0., eri, b, a, d, a); 
    c12c16term(m_ak, 0., eri, a, b, a, d); 
    c12c16term(m_ak, eri, 0., a, b, d, a); 
    c12c16term(m_ak, eri, 0., a, d, b, a); 
    c12c16term(m_ak, 0., eri, a, d, a, b); 
    c12c16term(m_ak, 0., eri, d, a, b, a); 
    c12c16term(m_ak, eri, 0., d, a, a, b); 


    eq30term(rho, eri, b, d, a, a); 
    eq30term(rho, eri, d, b, a, a); 
    eq30term(rho, eri, a, a, b, d); 
    eq30term(rho, eri, a, a, d, b); 
    eq30term(rho, eri, b, a, a, d); 
    eq30term(rho, eri, a, b, d, a); 
    eq30term(rho, eri, a, d, b, a); 
    eq30term(rho, eri, d, a, a, b); 

#endif

    

  } else if ((a > c) && (c == b) && (c > d)) { //8
    //cout << "case8" << endl;
#ifndef OLD_ID_CODE 
    eq20(rho, &idSelf_energy::eq20term_7, eri, b, a, b, d); 
    eq_c12(m_ak, &idSelf_energy::eq_c12term_7, eri, b, a, b, d);
    eq_c16(m_ak, &idSelf_energy::eq_c16term_7, eri, b, a, b, d);
    eq_30(rho, &idSelf_energy::eq_30term_7, eri, b, a, b, d);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri,eri, a, d, b, b); 
    eq20term(rho, eri,eri, d, a, b, b); 
    eq20term(rho, eri,eri, b, b, a, d); 
    eq20term(rho, eri,eri, b, b, d, a); 
    eq20term(rho, eri, 0., a, b, b, d); 
    eq20term(rho, 0., eri, a, b, d, b); 
    eq20term(rho, 0., eri, b, a, b, d); 
    eq20term(rho, eri, 0., b, a, d, b); 
    eq20term(rho, eri, 0., b, d, a, b); 
    eq20term(rho, 0., eri, b, d, b, a); 
    eq20term(rho, 0., eri, d, b, a, b); 
    eq20term(rho, eri, 0., d, b, b, a); 


    c12c16term(m_ak, eri,eri, a, d, b, b); 
    c12c16term(m_ak, eri,eri, d, a, b, b); 
    c12c16term(m_ak, eri,eri, b, b, a, d); 
    c12c16term(m_ak, eri,eri, b, b, d, a); 
    c12c16term(m_ak, eri, 0., a, b, b, d); 
    c12c16term(m_ak, 0., eri, a, b, d, b); 
    c12c16term(m_ak, 0., eri, b, a, b, d); 
    c12c16term(m_ak, eri, 0., b, a, d, b); 
    c12c16term(m_ak, eri, 0., b, d, a, b); 
    c12c16term(m_ak, 0., eri, b, d, b, a); 
    c12c16term(m_ak, 0., eri, d, b, a, b); 
    c12c16term(m_ak, eri, 0., d, b, b, a); 


    eq30term(rho, eri, a, d, b, b); 
    eq30term(rho, eri, d, a, b, b); 
    eq30term(rho, eri, b, b, a, d); 
    eq30term(rho, eri, b, b, d, a); 
    eq30term(rho, eri, a, b, b, d); 
    eq30term(rho, eri, b, a, d, b); 
    eq30term(rho, eri, b, d, a, b); 
    eq30term(rho, eri, d, b, b, a); 
#endif
    
  } else if ((a > c) && (c > b) && (b == d)) { //9
    //cout << "case9" << endl;
#ifndef OLD_ID_CODE
    eq20(rho, &idSelf_energy::eq20term_7, eri, b, a, b, c); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_7, eri, b, a, b, c);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_7, eri, b, a, b, c);
   eq_30(rho, &idSelf_energy::eq_30term_7, eri, b, a, b, c);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri,eri, a, c, b, b); 
    eq20term(rho, eri,eri, c, a, b, b); 
    eq20term(rho, eri,eri, b, b, a, c); 
    eq20term(rho, eri,eri, b, b, c, a); 
    eq20term(rho, eri, 0., a, b, b, c); 
    eq20term(rho, 0., eri, a, b, c, b); 
    eq20term(rho, 0., eri, b, a, b, c); 
    eq20term(rho, eri, 0., b, a, c, b); 
    eq20term(rho, eri, 0., b, c, a, b); 
    eq20term(rho, 0., eri, b, c, b, a); 
    eq20term(rho, 0., eri, c, b, a, b); 
    eq20term(rho, eri, 0., c, b, b, a); 

    c12c16term(m_ak, eri,eri, a, c, b, b); 
    c12c16term(m_ak, eri,eri, c, a, b, b); 
    c12c16term(m_ak, eri,eri, b, b, a, c); 
    c12c16term(m_ak, eri,eri, b, b, c, a); 
    c12c16term(m_ak, eri, 0., a, b, b, c); 
    c12c16term(m_ak, 0., eri, a, b, c, b); 
    c12c16term(m_ak, 0., eri, b, a, b, c); 
    c12c16term(m_ak, eri, 0., b, a, c, b); 
    c12c16term(m_ak, eri, 0., b, c, a, b); 
    c12c16term(m_ak, 0., eri, b, c, b, a); 
    c12c16term(m_ak, 0., eri, c, b, a, b); 
    c12c16term(m_ak, eri, 0., c, b, b, a); 

    eq30term(rho, eri, a, c, b, b); 
    eq30term(rho, eri, c, a, b, b); 
    eq30term(rho, eri, b, b, a, c); 
    eq30term(rho, eri, b, b, c, a); 
    eq30term(rho, eri, a, b, b, c); 
    eq30term(rho, eri, b, a, c, b); 
    eq30term(rho, eri, b, c, a, b); 
    eq30term(rho, eri, c, b, b, a); 
#endif


    

  } else if ((a > c) && (c == d) && (d > b)) { //10
        //cout << "case10" << endl;
#ifndef OLD_ID_CODE
    eq20(rho, &idSelf_energy::eq20term_6, eri, c, c, a, b); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_6, eri, c, c, a, b);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_6, eri, c, c, a, b);
   eq_30(rho, &idSelf_energy::eq_30term_6, eri, c, c, a, b);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, 0., c, b, c, a); 
    eq20term(rho, 0., eri, b, c, c, a); 
    eq20term(rho, eri, 0., c, a, c, b); 
    eq20term(rho, 0., eri, a, c, c, b); 
    eq20term(rho, eri, 0., b, c, a, c); 
    eq20term(rho, 0., eri, c, b, a, c); 
    eq20term(rho, eri, 0., a, c, b, c); 
    eq20term(rho, 0., eri, c, a, b, c); 


    c12c16term(m_ak, eri, 0., c, b, c, a); 
    c12c16term(m_ak, 0., eri, b, c, c, a); 
    c12c16term(m_ak, eri, 0., c, a, c, b); 
    c12c16term(m_ak, 0., eri, a, c, c, b); 
    c12c16term(m_ak, eri, 0., b, c, a, c); 
    c12c16term(m_ak, 0., eri, c, b, a, c); 
    c12c16term(m_ak, eri, 0., a, c, b, c); 
    c12c16term(m_ak, 0., eri, c, a, b, c); 

    eq30term(rho, eri, c, b, c, a); 
    eq30term(rho, eri, c, a, c, b); 
    eq30term(rho, eri, b, c, a, c); 
    eq30term(rho, eri, a, c, b, c); 
#endif
    
  } else if ((a > b) && (b > c) && (c == d)) { //11
    //cout << "case11" << endl;
#ifndef OLD_ID_CODE   
    eq20(rho, &idSelf_energy::eq20term_6, eri, c, c, a, b); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_6, eri, c, c, a, b);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_6, eri, c, c, a, b);
       eq_30(rho, &idSelf_energy::eq_30term_6, eri, c, c, a, b);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, 0., c, b, c, a); 
    eq20term(rho, 0., eri, b, c, c, a); 
    eq20term(rho, eri, 0., c, a, c, b); 
    eq20term(rho, 0., eri, a, c, c, b); 
    eq20term(rho, eri, 0., b, c, a, c); 
    eq20term(rho, 0., eri, c, b, a, c); 
    eq20term(rho, eri, 0., a, c, b, c); 
    eq20term(rho, 0., eri, c, a, b, c); 


    c12c16term(m_ak, eri, 0., c, b, c, a); 
    c12c16term(m_ak, 0., eri, b, c, c, a); 
    c12c16term(m_ak, eri, 0., c, a, c, b); 
    c12c16term(m_ak, 0., eri, a, c, c, b); 
    c12c16term(m_ak, eri, 0., b, c, a, c); 
    c12c16term(m_ak, 0., eri, c, b, a, c); 
    c12c16term(m_ak, eri, 0., a, c, b, c); 
    c12c16term(m_ak, 0., eri, c, a, b, c); 


    eq30term(rho, eri, c, b, c, a); 
    eq30term(rho, eri, c, a, c, b); 
    eq30term(rho, eri, b, c, a, c); 
    eq30term(rho, eri, a, c, b, c); 
#endif  

  } else if ((a > b) && (b > c) && (c > d)) {  //12
    //cout << "case12" << endl;
#ifndef OLD_ID_CODE  
    eq20(rho, &idSelf_energy::eq20term_12, eri, a, b, c, d); 
    eq_c12(m_ak, &idSelf_energy::eq_c12term_12, eri, a, b, c, d);
   eq_c16(m_ak, &idSelf_energy::eq_c16term_12, eri, a, b, c, d);
   eq_30(rho, &idSelf_energy::eq_30term_12, eri, a, b, c, d);
#endif
#ifdef OLD_ID_CODE

    eq20term(rho, eri, 0., b, d, a, c); 
    eq20term(rho, 0., eri, d, b, a, c); 
    eq20term(rho, eri, 0., b, c, a, d); 
    eq20term(rho, 0., eri, c, b, a, d); 
    eq20term(rho, eri, 0., a, d, b, c); 
    eq20term(rho, 0., eri, d, a, b, c); 
    eq20term(rho, eri, 0., a, c, b, d); 
    eq20term(rho, 0., eri, c, a, b, d); 
    eq20term(rho, eri, 0., d, b, c, a); 
    eq20term(rho, 0., eri, b, d, c, a); 
    eq20term(rho, eri, 0., d, a, c, b); 
    eq20term(rho, 0., eri, a, d, c, b); 
    eq20term(rho, eri, 0., c, b, d, a); 
    eq20term(rho, 0., eri, b, c, d, a); 
    eq20term(rho, eri, 0., c, a, d, b); 
    eq20term(rho, 0., eri, a, c, d, b); 


    c12c16term(m_ak, eri, 0., b, d, a, c); 
    c12c16term(m_ak, 0., eri, d, b, a, c); 
    c12c16term(m_ak, eri, 0., b, c, a, d); 
    c12c16term(m_ak, 0., eri, c, b, a, d); 
    c12c16term(m_ak, eri, 0., a, d, b, c); 
    c12c16term(m_ak, 0., eri, d, a, b, c); 
    c12c16term(m_ak, eri, 0., a, c, b, d); 
    c12c16term(m_ak, 0., eri, c, a, b, d); 
    c12c16term(m_ak, eri, 0., d, b, c, a); 
    c12c16term(m_ak, 0., eri, b, d, c, a); 
    c12c16term(m_ak, eri, 0., d, a, c, b); 
    c12c16term(m_ak, 0., eri, a, d, c, b); 
    c12c16term(m_ak, eri, 0., c, b, d, a); 
    c12c16term(m_ak, 0., eri, b, c, d, a); 
    c12c16term(m_ak, eri, 0., c, a, d, b); 
    c12c16term(m_ak, 0., eri, a, c, d, b); 


    eq30term(rho, eri, b, d, a, c); 
    eq30term(rho, eri, b, c, a, d); 
    eq30term(rho, eri, a, d, b, c); 
    eq30term(rho, eri, a, c, b, d); 
    eq30term(rho, eri, d, b, c, a); 
    eq30term(rho, eri, d, a, c, b); 
    eq30term(rho, eri, c, b, d, a); 
    eq30term(rho, eri, c, a, d, b); 

#endif

    
  } else if ((a > c) && (c > d) && (d > b)) {  //13
    //cout << "case13" << endl;
#ifndef OLD_ID_CODE  
    eq20(rho, &idSelf_energy::eq20term_12, eri, a, b, c, d); 
    eq_c12(m_ak, &idSelf_energy::eq_c12term_12, eri, a, b, c, d);
    eq_c16(m_ak, &idSelf_energy::eq_c16term_12, eri, a, b, c, d);
    eq_30(rho, &idSelf_energy::eq_30term_12, eri, a, b, c, d);
#endif
#ifdef OLD_ID_CODE
    eq20term(rho, eri, 0., b, d, a, c); 
    eq20term(rho, 0., eri, d, b, a, c); 
    eq20term(rho, eri, 0., b, c, a, d); 
    eq20term(rho, 0., eri, c, b, a, d); 
    eq20term(rho, eri, 0., a, d, b, c); 
    eq20term(rho, 0., eri, d, a, b, c); 
    eq20term(rho, eri, 0., a, c, b, d); 
    eq20term(rho, 0., eri, c, a, b, d); 
    eq20term(rho, eri, 0., d, b, c, a); 
    eq20term(rho, 0., eri, b, d, c, a); 
    eq20term(rho, eri, 0., d, a, c, b); 
    eq20term(rho, 0., eri, a, d, c, b); 
    eq20term(rho, eri, 0., c, b, d, a); 
    eq20term(rho, 0., eri, b, c, d, a); 
    eq20term(rho, eri, 0., c, a, d, b); 
    eq20term(rho, 0., eri, a, c, d, b); 


    c12c16term(m_ak, eri, 0., b, d, a, c); 
    c12c16term(m_ak, 0., eri, d, b, a, c); 
    c12c16term(m_ak, eri, 0., b, c, a, d); 
    c12c16term(m_ak, 0., eri, c, b, a, d); 
    c12c16term(m_ak, eri, 0., a, d, b, c); 
    c12c16term(m_ak, 0., eri, d, a, b, c); 
    c12c16term(m_ak, eri, 0., a, c, b, d); 
    c12c16term(m_ak, 0., eri, c, a, b, d); 
    c12c16term(m_ak, eri, 0., d, b, c, a); 
    c12c16term(m_ak, 0., eri, b, d, c, a); 
    c12c16term(m_ak, eri, 0., d, a, c, b); 
    c12c16term(m_ak, 0., eri, a, d, c, b); 
    c12c16term(m_ak, eri, 0., c, b, d, a); 
    c12c16term(m_ak, 0., eri, b, c, d, a); 
    c12c16term(m_ak, eri, 0., c, a, d, b); 
    c12c16term(m_ak, 0., eri, a, c, d, b); 

    eq30term(rho, eri, b, d, a, c); 
    eq30term(rho, eri, b, c, a, d); 
    eq30term(rho, eri, a, d, b, c); 
    eq30term(rho, eri, a, c, b, d); 
    eq30term(rho, eri, d, b, c, a); 
    eq30term(rho, eri, d, a, c, b); 
    eq30term(rho, eri, c, b, d, a); 
    eq30term(rho, eri, c, a, d, b); 
#endif


    
    
  } else if ((a > c) && (c > b) && (b > d)) {  //14

    //cout << "case14" << endl;
#ifndef OLD_ID_CODE  
    eq20(rho, &idSelf_energy::eq20term_12, eri, a, b, c, d); 
   eq_c12(m_ak, &idSelf_energy::eq_c12term_12, eri, a, b, c, d);
    eq_c16(m_ak, &idSelf_energy::eq_c16term_12, eri, a, b, c, d);
   eq_30(rho, &idSelf_energy::eq_30term_12, eri, a, b, c, d);
#endif
#ifdef OLD_ID_CODE
    
    eq20term(rho, eri, 0., b, d, a, c); 
    eq20term(rho, 0., eri, d, b, a, c); 
    eq20term(rho, eri, 0., b, c, a, d); 
    eq20term(rho, 0., eri, c, b, a, d); 
    eq20term(rho, eri, 0., a, d, b, c); 
    eq20term(rho, 0., eri, d, a, b, c); 
    eq20term(rho, eri, 0., a, c, b, d); 
    eq20term(rho, 0., eri, c, a, b, d); 
    eq20term(rho, eri, 0., d, b, c, a); 
    eq20term(rho, 0., eri, b, d, c, a); 
    eq20term(rho, eri, 0., d, a, c, b); 
    eq20term(rho, 0., eri, a, d, c, b); 
    eq20term(rho, eri, 0., c, b, d, a); 
    eq20term(rho, 0., eri, b, c, d, a); 
    eq20term(rho, eri, 0., c, a, d, b); 
    eq20term(rho, 0., eri, a, c, d, b); 

    c12c16term(m_ak, eri, 0., b, d, a, c); 
    c12c16term(m_ak, 0., eri, d, b, a, c); 
    c12c16term(m_ak, eri, 0., b, c, a, d); 
    c12c16term(m_ak, 0., eri, c, b, a, d); 
    c12c16term(m_ak, eri, 0., a, d, b, c); 
    c12c16term(m_ak, 0., eri, d, a, b, c); 
    c12c16term(m_ak, eri, 0., a, c, b, d); 
    c12c16term(m_ak, 0., eri, c, a, b, d); 
    c12c16term(m_ak, eri, 0., d, b, c, a); 
    c12c16term(m_ak, 0., eri, b, d, c, a); 
    c12c16term(m_ak, eri, 0., d, a, c, b); 
    c12c16term(m_ak, 0., eri, a, d, c, b); 
    c12c16term(m_ak, eri, 0., c, b, d, a); 
    c12c16term(m_ak, 0., eri, b, c, d, a); 
    c12c16term(m_ak, eri, 0., c, a, d, b); 
    c12c16term(m_ak, 0., eri, a, c, d, b); 
    
    eq30term(rho, eri, b, d, a, c); 
    eq30term(rho, eri, b, c, a, d); 
    eq30term(rho, eri, a, d, b, c); 
    eq30term(rho, eri, a, c, b, d); 
    eq30term(rho, eri, d, b, c, a); 
    eq30term(rho, eri, d, a, c, b); 
    eq30term(rho, eri, c, b, d, a); 
    eq30term(rho, eri, c, a, d, b); 
#endif
    
  } 
}



#ifdef OLD_ID_CODE

void idSelf_energy::eq20term(Blas_matrix& rho, double eri1, double eri2, unsigned a, unsigned c, unsigned d, unsigned b) 
{


  for (unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned k1_ = 0; k1_ < occs_[sym].size() ; k1_++) {
	unsigned k = occs_[sym][k_];
	unsigned k1 = occs_[sym][k1_];
	double fA = 0.;
	
	
	for(int m = 0; m < phis_.number_occupied();m++) {
	  double pole = 1. / ((phis_.energy(a)+phis_.energy(c)-phis_.energy(k1)-phis_.energy(m))
			      *(phis_.energy(a)+phis_.energy(c)-phis_.energy(k)-phis_.energy(m)));
	  double v_ack1m = V1212(a,c,k1,m) * pole;
	  double v_acmk1 = V1212(a,c,m,k1) * pole;
	  double exp1 = v_acmk1 - 2.*v_ack1m;
	  double exp2 = v_ack1m - 2.*v_acmk1;
	  
	  
	  fA += V1212(d,b,k,m) * (eri1 * exp1 )
 	    /(phis_.energy(k)+phis_.energy(m)-phis_.energy(d)-phis_.energy(b));
	  
	  
	  fA += V1212(d,b,k,m) *( eri2 * exp2) 
 	    /(phis_.energy(k)+phis_.energy(m)-phis_.energy(d)-phis_.energy(b));
	  
	}
	
	double fkk = 0.5 * fA;
	
	rho(k,k1) += fkk; 
	rho(k1,k) += fkk;
      }
}






void idSelf_energy::c12c16term(Blas_matrix& M_ak, double eri1, double eri2, unsigned a, unsigned b, unsigned c, unsigned d) 
{
  unsigned num_occ = phis_.number_occupied();
  
  // eq.C12
  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++) 
    for(unsigned p_ = 0; p_ < virs_[sym].size(); p_++) 
      for(unsigned q_ = 0; q_ < occs_[sym].size(); q_++) {
	unsigned p = virs_[sym][p_];
	unsigned q = occs_[sym][q_];
	
	double e_q = phis_.energy(q);
	
	double C1pq = 0.; 

	
	// eqs. (C12,C13)
	for(int i = 0; i < phis_.number_occupied(); i++) {
	  
	  double v_piab = (2. * V1212(p,i,a,b) - V1212(p,i,b,a))
	    / (e_q+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
	  
	  C1pq += v_piab * eri1 * V1212(q,i,c,d)
	    / (e_q+phis_.energy(i)-phis_.energy(c)-phis_.energy(d));
	}

	
	M_ak(p-num_occ,q) += C1pq;

      }

  //eq.C16
  unsigned p = a;
  a = b; b = c; c = d;
  unsigned sym = phis_.irrep(p);
  double e_p = phis_.energy(p);
  
  for(unsigned q_ = 0; q_ < occs_[sym].size(); q_++) {
    unsigned q = occs_[sym][q_];
    double C5pq = 0.;      
    for(int i = 0; i < phis_.number_occupied();i++)
      for(int j = 0; j < phis_.number_occupied();j++) {
	
	double v_qaij = V1212(q,a,i,j)
	  / (e_p+phis_.energy(a)-phis_.energy(i)-phis_.energy(j));
	
	C5pq += (2. * eri1) * V1212(i,j,b,c) * v_qaij
	  / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
	
	C5pq += (- eri2) * V1212(i,j,b,c) * v_qaij
	  / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
	
      }
    
    M_ak(p-num_occ,q) += C5pq;
  }
  
  
}




void idSelf_energy::eq30term(Blas_matrix& rho, double eri, unsigned a, unsigned b, unsigned c, unsigned d) 
{
  
  unsigned sym = phis_.irrep(b);  
  // f_store must have been initialized by particle_3_res
  Blas_matrix& f1 = f_store[sym];
  
  
  for(unsigned b_ = 0; b_ < virs_[sym].size(); b_++) {
    
    unsigned b_test = virs_[sym][b_];
    if (b_test != b) continue;
    
    
    for(unsigned s_ = 0; s_ < sats_[sym].size(); s_++) {
      
      unsigned a_test = sats_[sym][s_].a;
      
      if (a_test != a) continue;
      
      unsigned k = sats_[sym][s_].k;
      unsigned l = sats_[sym][s_].l;
      unsigned type = sats_[sym][s_].type;
       
      double gB = 0.;
      
      if (k == l) {
	gB += V1212(c,d,l,k) * eri
	  / (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
      } else if (type == 0) {
	gB += eri * (V1212(c,d,l,k) + V1212(c,d,k,l))
	  / (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
	gB *= SQRT_1_2;
      } else {
	gB += (V1212(c,d,k,l) - V1212(c,d,l,k)) * eri
	  / (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
	gB *= SQRT_3_2;
      }
      gB /= (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
      
      for(unsigned v_ = 0; v_ < virs_[sym].size(); v_++) {
	unsigned v = virs_[sym][v_];
	double mult = gB*f1(s_,v_);
	rho(v,b) += mult;
	rho(b,v) += mult;
      }
    }
  }
}


#endif



