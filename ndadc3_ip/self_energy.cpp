#include "self_energy.hpp"
#include "integral_table.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "blas_matrix.hpp"


static const double SQRT_1_2 = 0.707106781186547;
static const double SQRT_1_6 = 0.408248290463863;
static const double SQRT_3_2 = 1.22474487139159;
static int chunk = 1;

//Result from the pole test: appears to speed up the computation around 10%
//vector<vector<vector<vector<double> > > > poles;

inline double Self_energy::V1212(unsigned a, unsigned b, unsigned c, unsigned d)
{
  return table_.integral(a,c,b,d);
}


Self_energy::Self_energy(Integral_table& tb, SCF_data_reader& r)
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

 
//   vector<double> v1(phis_.number_orbitals(), 0.);
//   vector<vector<double> > v2(phis_.number_orbitals(), v1);
//   vector<vector<vector<double> > > v3(phis_.number_orbitals(), v2);
//   poles.assign(phis_.number_orbitals(), v3);
  
//   for(unsigned a = phis_.number_occupied(); a < phis_.number_orbitals(); a++) 
//     for(unsigned b = phis_.number_occupied(); b < phis_.number_orbitals(); b++) 
//       for(unsigned i = 0; i < phis_.number_occupied(); i++)
//  	for(unsigned j = 0; j < phis_.number_occupied(); j++) 
//  	  poles[i][j][a][b] = 1. / (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(b));





}

//Computes the density up to 2nd order
void  Self_energy::density_2(Blas_matrix& rho) 
{

  //Second order corrections
  rho_hole_part_2(rho);
  rho_holeparticle_2(rho);
  rho_particle_part_2(rho);  

  //Zeroth order density
  for(unsigned k = 0; k < phis_.number_occupied(); k++)
    rho(k,k) += 1.;
}
//Computes the density up to 3rd order
void  Self_energy::density_3(Blas_matrix& rho) 
{
  unsigned num_occ = phis_.number_occupied();
  unsigned num_vir = phis_.number_orbitals() - num_occ;

  //Second order corrections
  rho_hole_part_2(rho);
  rho_holeparticle_2(rho);
  rho_particle_part_2(rho);  

  // Compute the static self-energy contributions:
  Blas_matrix sigma3(rho);
  rho2sigma(sigma3);
  // Compute the dynamic self-energy contributions:
  Blas_matrix m_ak(num_vir, num_occ);
  m_ak = 0.;
  dynamic_self_energy_3(m_ak);
  m_ak(0,num_vir,0,num_occ).daxpy(1., sigma3(num_occ,num_vir, 0, num_occ));

  //Third order corrections
  rho_hole_part_3(rho);
  rho_holeparticle_3(rho, m_ak);
  rho_particle_part_3(rho);  

  //Zeroth order density
  for(unsigned k = 0; k < phis_.number_occupied(); k++)
    rho(k,k) += 1.;  
}
//Computes the density up to 3nd order, but uses sigma4+ 
// for its static self-energy contribution
void Self_energy::density_3plus(Blas_matrix& rho)
{
  unsigned num_occ = phis_.number_occupied();
  unsigned num_vir = phis_.number_orbitals() - num_occ;
  
  //Get the density correction, 2nd order
  rho_hole_part_2(rho);
  rho_holeparticle_2(rho);
  rho_particle_part_2(rho);
  
  // Compute the dynamic self-energy contributions:
  Blas_matrix m_ak(num_vir, num_occ);
  m_ak = 0.;
  dynamic_self_energy_3(m_ak);

  //Get the density correction, 3rd order
  rho_hole_part_3(rho);  
  rho_holeparticle_3(rho, m_ak); // no static contribution here
  rho_particle_part_3(rho);  

  //Get the inhomogeneities
  Blas_matrix sigma4p(rho);
  rho2sigma(sigma4p);

  //Solve for sigma 4+
  linear_eq_selfenergy(sigma4p);

  //Add the static self-energy contribution of sigma4+
  m_ak = 0.;
  m_ak(0,num_vir,0,num_occ).daxpy(1., sigma4p(num_occ,num_vir, 0, num_occ));
  rho_holeparticle_3(rho, m_ak); // just the static contribution here

  //Zeroth order density
  for(unsigned k = 0; k < phis_.number_occupied(); k++)
    rho(k,k) += 1.;
}
// Computes the 3rd order static self energy
void Self_energy::static_selfenergy_3(Blas_matrix& sigma)
{
  //Get the density correction, 2nd order
  rho_hole_part_2(sigma);
  rho_holeparticle_2(sigma);
  rho_particle_part_2(sigma);
  //Get the third order static self-energy
  rho2sigma(sigma);
}

// Computes the 4th order static self energy
void Self_energy::static_selfenergy_4(Blas_matrix& sigma)
{
  unsigned num_occ = phis_.number_occupied();
  unsigned num_vir = phis_.number_orbitals() - num_occ;

  //Get the density correction, 2nd order
  rho_hole_part_2(sigma);
  rho_holeparticle_2(sigma);
  rho_particle_part_2(sigma);
  
  // Compute the static self-energy contributions:
  Blas_matrix sigma3(sigma);
  rho2sigma(sigma3);
  // Compute the dynamic self-energy contributions:
  Blas_matrix m_ak(num_vir, num_occ);
  m_ak = 0.;
  dynamic_self_energy_3(m_ak);
  m_ak(0,num_vir,0,num_occ).daxpy(1., sigma3(num_occ,num_vir, 0, num_occ));
    
  //Get the density correction, 3rd order
  rho_hole_part_3(sigma);
  rho_holeparticle_3(sigma, m_ak); //both static and dynamic contributions here
  rho_particle_part_3(sigma);  
  
  //Get the fourth order static self-energy
  rho2sigma(sigma);
}
// Computes the 4+ static self energy
void Self_energy::static_selfenergy_4plus(Blas_matrix& sigma)
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
  dynamic_self_energy_3(m_ak);


  //Get the density correction, 3rd order
  rho_hole_part_3(sigma);
  rho_holeparticle_3(sigma, m_ak);
  rho_particle_part_3(sigma);
  
  //Get the inhomogeneities
  rho2sigma(sigma);


  //Solve for sigma 4+
  linear_eq_selfenergy(sigma);


}

// This method implements the relation between 
// a lower order density to a higher order static self-energy.
// (Takes the density correction, i.e. without zeroth order, as input)
// eq.(A25), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
void Self_energy::rho2sigma(Blas_matrix& rho) {
  
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
void  Self_energy::rho_hole_part_2(Blas_matrix& rho)
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
void Self_energy::rho_holeparticle_2(Blas_matrix& rho)
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
void Self_energy::rho_particle_part_2(Blas_matrix& rho)
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
inline void Self_energy::rho_hole_part_3(Blas_matrix& rho) 
{
  
  for (unsigned sym = 0; sym < phis_.number_irreps(); sym++)
    for(unsigned k_ = 0; k_ < occs_[sym].size();k_++)
      for(unsigned k1_ = 0; k1_ < occs_[sym].size() ; k1_++) {
	unsigned k = occs_[sym][k_];
	unsigned k1 = occs_[sym][k1_];
	
	double fA = 0.;
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
	      
	      // Eq.(A20)

	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
		unsigned sym_d = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b),phis_.irrep(k)),phis_.irrep(m));
		for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
		  int d = virs_[sym_d][d_];

		  fA += V1212(d,b,k,m) * (V1212(a,c,d,b) * exp1 + V1212(a,c,b,d) * exp2) 
		    /(phis_.energy(k)+phis_.energy(m)-phis_.energy(d)-phis_.energy(b));
		  
		}
	      }
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
	double fkk = 0.5 * (fA + fC) + (fB + fD);
	//Eq.(A33)

	rho(k,k1) += fkk; 
	rho(k1,k) += fkk;


      }
}


// This method computes the third order density, hole/particle part,
// but only with the dynamic self energy contributions
// See eqs.(A24,A34), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
void Self_energy::rho_holeparticle_3(Blas_matrix& rho, Blas_matrix& self_energy)
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

//TODO: remove the f1, f2 quantites 
// This method computes the third order density, particle part
// See eqs.(A27-A32,A36), Schirmer et al, J. Chem. Phys., 109 p.4734 (1998).
// The expressions have been spin-adapted (see rho_particle_part_2) and simplified
void Self_energy::rho_particle_part_3(Blas_matrix& rho)
{

  for(unsigned int sym = 0; sym < phis_.number_irreps(); sym++) {
    
    Blas_matrix f1(sats_[sym].size(),virs_[sym].size());
    f1 = 0.;
#pragma omp parallel for schedule(static,chunk)
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
    

    Blas_matrix f2(sats_[sym].size(),virs_[sym].size());
    f2 = 0.;
#pragma omp parallel for schedule(static,chunk)
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

	if (k == l) {
	  
	  //eq. (A29)
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
	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
	    unsigned sym_d = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	    
	    for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
	      unsigned d = virs_[sym_d][d_];
	      
	      gB += V1212(c,d,l,k) * V1212(a,b,c,d)
		/ (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
	    }
	  }
	  // eq.(A31,A32)
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
	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
 	    unsigned sym_d = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	    
	    for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
	      unsigned d = virs_[sym_d][d_];
	      
	      gB += V1212(a,b,c,d) * (V1212(c,d,l,k) + V1212(c,d,k,l))
		/ (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
	    }
	  }
	  // eq.(A31,A32)
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
	      gA *= SQRT_1_2;
	      gB *= SQRT_1_2;
	      gC *= SQRT_1_2;
	      gD *= SQRT_1_2;
	      gC1 *= SQRT_1_2;
	      gD1 *= SQRT_1_2;
	      
	} else {
	  // eq.(29)
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
	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
	    unsigned sym_d = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),phis_.irrep(b)),phis_.irrep(c));
	    
	    for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
	      unsigned d = virs_[sym_d][d_];
	      
	      gB += (V1212(c,d,k,l) - V1212(c,d,l,k)) * V1212(a,b,c,d)
		/ (phis_.energy(c)+phis_.energy(d)-phis_.energy(k)-phis_.energy(l));
	    }
	  }
	  // eq.(A31,A32)
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
	  gA *= SQRT_3_2;
	  gB *= SQRT_3_2;
	  gC *= SQRT_3_2;
	  gD *= SQRT_3_2;
	  gC1 *= SQRT_3_2;
	  gD1 *= SQRT_3_2;
	}
	  

	f2(s_,b_) = (gA+gC+gB+gD+gC1+gD1)/(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(l));
	
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
void Self_energy::sigma_3diagrams(Blas_matrix& sigma)
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


// This method computes the third order dynamic self energy
// contributions eq.(C.12-C.23)in Appendix C, W. von Niessen et al.,
// Computer Physics Reports 1 (1984) 57-125.
// M^{(3)+}(\omega) includes eqs.(C.12-C.14,C.18-C.20)
// M^{(3)-}(\omega) includes eqs.(C.15-C.17,C.21-C.23)
void Self_energy::dynamic_self_energy_3(Blas_matrix& M_ak)
{
  unsigned num_occ = phis_.number_occupied();

  for(unsigned sym = 0; sym < phis_.number_irreps(); sym++)  
#pragma omp parallel for schedule(static,chunk) 
    for(unsigned p_ = 0; p_ < virs_[sym].size(); p_++)  {
      for(unsigned q_ = 0; q_ < occs_[sym].size(); q_++) {
      unsigned p = virs_[sym][p_];
      unsigned q = occs_[sym][q_];
      

      double e_p = phis_.energy(p);
      double e_q = phis_.energy(q);

      double C1pq = 0.; double C2pq = 0.; double C3pq = 0.;
      double C4pq = 0.; double C5pq = 0.; double C6pq = 0.;

      double D1pq = 0.; double D2pq = 0.; double D3pq = 0.;
      double D4pq = 0.; double D5pq = 0.; double D6pq = 0.;

      // eqs. (C12,C13)
      //#pragma omp parallel 
      //{
      //#pragma omp for reduction(+:C1pq, C2pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied(); i++) 
	  for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++) {
	    unsigned sym_b = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(a));
	    
	    for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	    
	      double v_piab = (2. * V1212(p,i,a,b) - V1212(p,i,b,a))
		/ (e_q+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
		//*poles[q][i][a][b];
	      
	      for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
		unsigned sym_d = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(c));
		
		for(int d_ = 0; d_ < virs_[sym_d].size(); d_++) {
		  unsigned d = virs_[sym_d][d_];
		  
		  C1pq += v_piab * V1212(a,b,c,d) * V1212(q,i,c,d)
		    / (e_q+phis_.energy(i)-phis_.energy(c)-phis_.energy(d));
		    //*poles[q][i][c][d];

		}
	      }
	      

	      
	      for(int j = 0; j < phis_.number_occupied();j++) {
		unsigned sym_k = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(j));
		
		for(int k_ = 0; k_ < occs_[sym_k].size(); k_++) {
		  unsigned k = occs_[sym_k][k_];
		  
		  C2pq += v_piab * V1212(a,b,j,k) * V1212(q,i,j,k)
		    / (phis_.energy(j)+phis_.energy(k)-phis_.energy(a)-phis_.energy(b));
		    //*poles[j][k][a][b];
		  
		}
	      }
	    }
	  }

	// eq.(C.14)
	//#pragma omp for reduction(+:C3pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied(); i++) 
	  for(int a = phis_.number_occupied(); a < phis_.number_orbitals();a++) {
	    unsigned sym_b = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(a));
	    
	    for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
	      unsigned b = virs_[sym_b][b_];
	      
	      double v_qiab = V1212(q,i,a,b)
		/ (e_q+phis_.energy(i)-phis_.energy(a)-phis_.energy(b));
		//*poles[q][i][a][b];

	      for(int k = 0; k < phis_.number_occupied();k++) {
		unsigned sym_j = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a),sym_b),phis_.irrep(k));
		
		for(int j_ = 0; j_ < occs_[sym_j].size(); j_++) {
		  unsigned j = occs_[sym_j][j_];
		  
		  C3pq += (2. * V1212(p,i,j,k) - V1212(p,i,k,j)) * V1212(a,b,j,k) * v_qiab
		    / (phis_.energy(j)+phis_.energy(k)-phis_.energy(a)-phis_.energy(b));
		    //*poles[j][k][a][b];

		}
	      }
	    }
	  } 
	// eqs. (C.15,C.17)
	//#pragma omp for reduction(+:C4pq, C6pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied();i++)
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(j));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      double v_paij = 
		(2. * V1212(p,a,i,j) - V1212(p,a,j,i))
		/ (e_p+phis_.energy(a)-phis_.energy(i)-phis_.energy(j));
		//(-2. * V1212(p,a,i,j) + V1212(p,a,j,i))
		//*(poles[i][j][p][a]);

	      
	      for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++) {
		unsigned sym_c = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),sym_a),phis_.irrep(b));
		
		for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
		  unsigned c = virs_[sym_c][c_];
		  
		  C4pq += v_paij * V1212(i,j,b,c) * V1212(q,a,b,c)
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
		    //*poles[i][j][b][c];
		  
		}
	      }
	      
	      for(int k = 0; k < phis_.number_occupied();k++) {
		unsigned sym_l = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),sym_a),phis_.irrep(k));
		
		for(int l_ = 0; l_ < occs_[sym_l].size(); l_++) {
		  unsigned l = occs_[sym_l][l_];
		  
		  C6pq += v_paij * V1212(i,j,k,l) * V1212(q,a,k,l)
		    / (e_p+phis_.energy(a)-phis_.energy(k)-phis_.energy(l));
		    //*(-poles[k][l][p][a]);

		  
		}
	      }		
	    }
	  }
	// eq.(C.16)
	//#pragma omp for reduction(+:C5pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied();i++)
	  for(int j = 0; j < phis_.number_occupied();j++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(i)),phis_.irrep(j));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      double v_qaij = V1212(q,a,i,j)
		/ (e_p+phis_.energy(a)-phis_.energy(i)-phis_.energy(j));
		//*(-poles[i][j][p][a]);

	      for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
		unsigned sym_b = 
		  phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),sym_a),phis_.irrep(c));
		
		for(int b_ = 0; b_ < virs_[sym_b].size(); b_++) {
		  unsigned b = virs_[sym_b][b_];
		  
		  C5pq += (2. * V1212(p,a,b,c) - V1212(p,a,c,b)) * V1212(i,j,b,c) * v_qaij
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(b)-phis_.energy(c));
		    //*poles[i][j][b][c];

		}
	      }
	    }
	  }
	
	
	//eqs.(C.18,C.20)
	//#pragma omp for reduction(+:D1pq, D3pq) schedule(static,chunk)  nowait
	for(int j = 0; j < phis_.number_occupied();j++) 
	  for(int b = phis_.number_occupied(); b < phis_.number_orbitals();b++){
	    unsigned sym_c = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(j)),phis_.irrep(b));
	    
	    for(int c_ = 0; c_ < virs_[sym_c].size(); c_++) {
	      unsigned c = virs_[sym_c][c_];
	      //intermediates for D1 and D3
	      double pole = 
		1. / (e_q + phis_.energy(j) - phis_.energy(b) - phis_.energy(c));
		//poles[q][j][b][c];

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
			      (v_ajic * exp3 + v_ajci * exp1)  )  
		    /(e_q + phis_.energy(i) - phis_.energy(a) - phis_.energy(b));
		    //*poles[q][i][a][b];

		  double v_ijac = V1212(i,j,a,c);
		  double v_ijca = V1212(i,j,c,a);
		  
		  D3pq += (   V1212(p,a,i,b) * 
			      (v_ijac * exp1 + v_ijca * exp2) +
			      V1212(p,a,b,i) * 
			      (v_ijac * exp3 + v_ijca * exp1)  )  
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(c)-phis_.energy(a));
		    //*poles[i][j][c][a];

		}
	      }
	    }
	  }
	//eq.(C.19)
	//#pragma omp for reduction(+:D2pq) schedule(static,chunk)  nowait
	for(int i = 0; i < phis_.number_occupied();i++) 
	  for(int c = phis_.number_occupied(); c < phis_.number_orbitals();c++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(p),phis_.irrep(i)),phis_.irrep(c));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      double pole = 
		1. / (e_q + phis_.energy(i) - phis_.energy(a) - phis_.energy(c));
		//poles[q][i][a][c];

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
					v_abji * exp2 )  ) 
		    / (phis_.energy(i)+phis_.energy(j)-phis_.energy(a)-phis_.energy(b));
		    //*poles[i][j][a][b];

		}
	      }
	    }
	  }
	//eqs(C.21,C.22)
	//#pragma omp for reduction(+:D6pq, D4pq, D5pq) schedule(static,chunk)  nowait
	for(int j = 0; j < phis_.number_occupied();j++) 
	  for(int k = 0; k < phis_.number_occupied();k++) {
	    unsigned sym_a = 
	      phis_.irrep_product(phis_.irrep_product(phis_.irrep(q),phis_.irrep(j)),phis_.irrep(k));
	    
	    for(int a_ = 0; a_ < virs_[sym_a].size(); a_++) {
	      unsigned a = virs_[sym_a][a_];
	      
	      // intermediates for D5 and D6
	      double pole = 
		1. / (e_p+phis_.energy(a)-phis_.energy(j)-phis_.energy(k));
		//-poles[j][k][p][a];

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
		    //*poles[i][j][a][b];

		  double v_iabj = V1212(i,a,b,j);
		  double v_iajb = V1212(i,a,j,b);
		  
		  D6pq += (V1212(p,b,k,i) * (v_iabj * exp3 + v_iajb * exp1)
			   +V1212(p,b,i,k) * (v_iabj * exp1 + v_iajb * exp2) )
		    / (e_p+phis_.energy(b)-phis_.energy(i)-phis_.energy(k));
		    //*(-poles[i][k][p][b]);

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
		    //*poles[i][j][a][b];

		}
	      }
	    }
	  }
	//}
      M_ak(p-num_occ,q) = 
	(C1pq + C2pq) + (C3pq + C4pq) + (C5pq - C6pq) +
	(D1pq + D2pq) + (D3pq + D4pq) + (D5pq - D6pq);
      
      }
   }  
}
 
 
 
// This method takes as input the inhomogeneity matrix B
// and produces the static self-interaction using the
// linear equation defined in Appendix B, W. von Niessen et al.,
// Computer Physics Reports 1 (1984) 57-125.
void Self_energy::linear_eq_selfenergy(Blas_matrix& b)  
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
 
