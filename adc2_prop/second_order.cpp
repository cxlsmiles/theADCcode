#include "second_order.hpp"
#include "adc2_dip/config.hpp"
#include "matrices.hpp"
#include "ndadc3_ip/self_energy.hpp"
#include "blas_matrix.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

//static const double SQRT_2 =  1.41421356237310;
static const double SQRT_1_2 = 0.707106781186548;





inline double Second_order_cap::V1212(unsigned a, unsigned b, unsigned c, unsigned d)
{
  return V1122(a,c,b,d);
}


Second_order_cap::Second_order_cap(SCF_data_reader& phis, Integral_table& tab, Triangular_matrix<double>* mat) 
  : phis_(phis), gs_expt_value2_(0.), d(mat->rows()),  ADC2_DIP_blocks(phis, tab)
{
  
  for(unsigned row = 0; row < mat->rows(); row++)
    for(unsigned col = 0; col <= row; col++) {
      d(row, col) = (*mat)(row,col);
    }
  
  //Group the occupied and virtual orbitals according to symmetry;
  for(unsigned i = 0; i < phis_.number_occupied();  i++) 
    occs_[phis_.irrep(i)].push_back(i); 

  for(unsigned i = phis_.number_occupied(); i < phis_.number_orbitals();  i++)
    virs_[phis_.irrep(i)].push_back(i); 

  rho = new Blas_matrix(phis_.number_orbitals(), phis_.number_orbitals());
  *rho = 0.;

  
  Self_energy(*int_block_, phis_).density_3plus(*rho);
  //Self_energy(*int_block_, phis_).density_2(*rho);

 

  for(unsigned i = 0; i < phis_.number_occupied();  i++) 
    rho->operator()(i,i) -= 1.; // get the density correction


 

  for(unsigned p = 0; p < phis_.number_orbitals();  p++) 
    for(unsigned q = 0; q < phis_.number_orbitals();  q++) {
      gs_expt_value2_ += d(p,q)*rho->operator()(q,p);
    }
  gs_expt_value2_ *= 2.;

}


inline double Second_order_cap::d2_1(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != i1) return 0.;
  
  unsigned j_sym = phis_.irrep(j);
  unsigned j1_sym = phis_.irrep(j1);
  
  if (j_sym != j1_sym) return 0.;
  
  double d21 = 0.;
  for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) 
    for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++) {
      
      unsigned a1_sym = phis_.irrep(a);
      for(unsigned a1_ = 0; a1_ < virs_[a1_sym].size() ; a1_++) {
	unsigned a1 = virs_[a1_sym][a1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, phis_.irrep(b)), j1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  

 	  d21 += (  2.*V1212(a1,b,k,j1)*V1212(a,b,k,j)
 		    + 2.*V1212(a1,b,j1,k)*V1212(a,b,j,k)
 		    - V1212(a1,b,j1,k)*V1212(a,b,k,j)
 		    - V1212(a1,b,k,j1)*V1212(a,b,j,k)) * d(a1,a) /
 	    (  (phis_.energy(a1)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j1)) 
 	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));


	  
	}
      }
    }
  
  d21 *= -.5;
  
  if ((i == j) && (i1 == j1)) return d21;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d21;
  
  return d21;
  
}


inline double Second_order_cap::d2_2(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != i1) return 0.;
  
  unsigned j_sym = phis_.irrep(j);
  unsigned j1_sym = phis_.irrep(j1);
  
  if (j_sym != j1_sym) return 0.;

  
  double d22 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = j1_sym;
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), k1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d22 += (  2.*V1212(a,b,k1,k)*V1212(a,b,j,k)
		    + 2.*V1212(a,b,k,k1)*V1212(a,b,k,j)
		    - V1212(a,b,k1,k)*V1212(a,b,k,j)
		    - V1212(a,b,k,k1)*V1212(a,b,j,k)) * d(j1,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));
	  
	}
      }
    }
  
  d22 *= .25;
  
  if ((i == j) && (i1 == j1)) return d22;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d22;
  
  return d22;
  
}


inline double Second_order_cap::d2_3(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != i1) return 0.;
  
  unsigned j_sym = phis_.irrep(j);
  unsigned j1_sym = phis_.irrep(j1);
  
  if (j_sym != j1_sym) return 0.;

  
  double d23 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), j1_sym);
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = k1_sym;
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d23 += (  2.*V1212(a,b,j1,k1)*V1212(a,b,j,k)
		    + 2.*V1212(a,b,k1,j1)*V1212(a,b,k,j)
		    - V1212(a,b,j1,k1)*V1212(a,b,k,j)
		    - V1212(a,b,k1,j1)*V1212(a,b,j,k)) * d(k,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k1)-phis_.energy(j1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));

	}
      }
    }
  
  d23 *= .25;
  
  if ((i == j) && (i1 == j1)) return d23;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d23;
  
  return d23;
  
}



inline double Second_order_cap::d2_4(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (j != j1) return 0.;
  
  unsigned i_sym = phis_.irrep(i);
  unsigned i1_sym = phis_.irrep(i1);
  
  if (i_sym != i1_sym) return 0.;

  //cout << "iji1j1 " << i << ' ' << j << ' ' << i1 << ' ' << j1 << endl;
  
  double d24 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned a1_sym = phis_.irrep(a);
      for(unsigned a1_ = 0; a1_ < virs_[a1_sym].size() ; a1_++) {
	unsigned a1 = virs_[a1_sym][a1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, phis_.irrep(b)), i1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d24 += (  2.*V1212(a1,b,k,i1)*V1212(a,b,k,i)
		    + 2.*V1212(a1,b,i1,k)*V1212(a,b,i,k)
		    - V1212(a1,b,i1,k)*V1212(a,b,k,i)
		    - V1212(a1,b,k,i1)*V1212(a,b,i,k)) * d(a1,a) /
	    (  (phis_.energy(a1)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  //cout << a << ' ' << a1 << ' '<< d(a1,a) << ' '<< d24 << endl;
	}
      }
    }

  
  d24 *= -.5;
  
  if ((i == j) && (i1 == j1)) return d24;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d24;
  
  return d24;
  
}




inline double Second_order_cap::d2_5(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (j != j1) return 0.;
  
  unsigned i_sym = phis_.irrep(i);
  unsigned i1_sym = phis_.irrep(i1);

  if (i_sym != i1_sym) return 0.;
  
  double d25 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = i1_sym;
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), k1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];

	  d25 += (  2.*V1212(a,b,k,k1)*V1212(a,b,k,i)
		    + 2.*V1212(a,b,k1,k)*V1212(a,b,i,k)
		    - V1212(a,b,k,k1)*V1212(a,b,i,k)
		    - V1212(a,b,k1,k)*V1212(a,b,k,i)) * d(i1,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  
	}
      }
    }
  
  d25 *= .25;
  
  if ((i == j) && (i1 == j1)) return d25;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d25;
  
  return d25;
  
}



inline double Second_order_cap::d2_6(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (j != j1) return 0.;

  unsigned i_sym = phis_.irrep(i);
  unsigned i1_sym = phis_.irrep(i1);

  if (i_sym != i1_sym) return 0.;
  
  double d26 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), i1_sym);
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = k1_sym;
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d26+= (  2.*V1212(a,b,k1,i1)*V1212(a,b,k,i)
		   + 2.*V1212(a,b,i1,k1)*V1212(a,b,i,k)
		   - V1212(a,b,k1,i1)*V1212(a,b,i,k)
		    - V1212(a,b,i1,k1)*V1212(a,b,k,i)) * d(k,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k1)-phis_.energy(i1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  
	}
      }
    }
  
  d26 *= .25;
  
  if ((i == j) && (i1 == j1)) return d26;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d26;
  
  return d26;
  
}


inline double Second_order_cap::d2_7(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != j1) return 0.;

  if ((i == j) &&(i1 == j1)) return 0.;

  
  unsigned j_sym = phis_.irrep(j);
  unsigned i1_sym = phis_.irrep(i1);

 
  if (j_sym != i1_sym) return 0.;
  
  double d27 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned a1_sym = phis_.irrep(a);
      for(unsigned a1_ = 0; a1_ < virs_[a1_sym].size() ; a1_++) {
	unsigned a1 = virs_[a1_sym][a1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, phis_.irrep(b)), i1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d27+= (  -2.*V1212(a1,b,k,i1)*V1212(a,b,k,j)
		   - 2.*V1212(a1,b,i1,k)*V1212(a,b,j,k)
		   + V1212(a1,b,k,i1)*V1212(a,b,j,k)
		   + V1212(a1,b,i1,k)*V1212(a,b,k,j)) * d(a1,a) /
	    (  (phis_.energy(a1)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));
	  
	}
      }
    }
  
  
  d27 *= .5;
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d27;
  
  return d27;
  
}




inline double Second_order_cap::d2_8(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != j1) return 0.;
  if ((i == j) &&(i1 == j1)) return 0.;
  
  unsigned i1_sym = phis_.irrep(i1);
  unsigned j_sym = phis_.irrep(j);
 
  if (i1_sym != j_sym) return 0.;
  
  double d28 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = i1_sym;
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), k1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d28+= (  -2.*V1212(a,b,k,k1)*V1212(a,b,k,j)
		   - 2.*V1212(a,b,k1,k)*V1212(a,b,j,k)
		   + V1212(a,b,k,k1)*V1212(a,b,j,k)
		   + V1212(a,b,k1,k)*V1212(a,b,k,j)) * d(i1,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));
	  
	}
      }
    }
  
  d28 *= -.25;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d28;
  
  return d28;
  
}



inline double Second_order_cap::d2_9(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i != j1) return 0.;
  if ((i == j) &&(i1 == j1)) return 0.;
  
  unsigned i1_sym = phis_.irrep(i1);
  unsigned j_sym = phis_.irrep(j);
 
  if (i1_sym != j_sym) return 0.;
  
  double d29 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym =  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), i1_sym);
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = k1_sym;
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d29+= (  -2.*V1212(a,b,k1,i1)*V1212(a,b,k,j)
		   - 2.*V1212(a,b,i1,k1)*V1212(a,b,j,k)
		   + V1212(a,b,k1,i1)*V1212(a,b,j,k)
		   + V1212(a,b,i1,k1)*V1212(a,b,k,j)) * d(k,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k1)-phis_.energy(i1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j)));
	  
	}
      }
    }
  
  d29 *= -.25;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d29;
  
  return d29;
  
}



inline double Second_order_cap::d2_10(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (j != i1) return 0.;
  if ((i == j) &&(i1 == j1)) return 0.;
  
  unsigned i_sym = phis_.irrep(i);
  unsigned j1_sym = phis_.irrep(j1);
  
  if (j1_sym != i_sym) return 0.;
  
  double d210 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned a1_sym = phis_.irrep(a);
      for(unsigned a1_ = 0; a1_ < virs_[a1_sym].size() ; a1_++) {
	unsigned a1 = virs_[a1_sym][a1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(a1_sym, phis_.irrep(b)), j1_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d210 += (  -2.*V1212(a1,b,k,j1)*V1212(a,b,k,i)
		     - 2.*V1212(a1,b,j1,k)*V1212(a,b,i,k)
		     + V1212(a1,b,k,j1)*V1212(a,b,i,k)
		     + V1212(a1,b,j1,k)*V1212(a,b,k,i)) * d(a1,a) /
	    (  (phis_.energy(a1)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  
	}
      }
    }
  
  d210 *= .5;

  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d210;
  
  return d210;
  
}



inline double Second_order_cap::d2_11(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (i1 != j) return 0.;
  if ((i == j) &&(i1 == j1)) return 0.;
  
  unsigned j1_sym = phis_.irrep(j1);
  unsigned i_sym = phis_.irrep(i);
  
  
  if (i_sym != j1_sym) return 0.;
  
  double d211 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym = j1_sym;
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), i_sym);
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d211+= (  -2.*V1212(a,b,k1,k)*V1212(a,b,i,k)
		    - 2.*V1212(a,b,k,k1)*V1212(a,b,k,i)
		    + V1212(a,b,k1,k)*V1212(a,b,k,i)
		    + V1212(a,b,k,k1)*V1212(a,b,i,k)) * d(j1,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(k1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  
	}
      }
    }
  
  
  d211 *= -.25;

  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d211;
  
  return d211;
  
}



inline double Second_order_cap::d2_12(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  if (j != i1) return 0.;
  if ((i == j) &&(i1 == j1)) return 0.;
  
  unsigned j1_sym = phis_.irrep(j1);
  unsigned i_sym = phis_.irrep(i);

  
  if (j1_sym != i_sym) return 0.;
  
  double d212 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++)   
    for(unsigned a = phis_.number_occupied(); a < number_orbitals(); a++) {
      unsigned k1_sym =  phis_.irrep_product(phis_.irrep_product(phis_.irrep(a), phis_.irrep(b)), j1_sym);
      for(unsigned k1_ = 0; k1_ < occs_[k1_sym].size(); k1_++) {
	unsigned k1 = occs_[k1_sym][k1_];
	unsigned k_sym = k1_sym;
	for(unsigned k_ = 0; k_ < occs_[k_sym].size(); k_++) {
	  unsigned k = occs_[k_sym][k_];
	  
	  d212+= (  -2.*V1212(a,b,k1,j1)*V1212(a,b,k,i)
		    - 2.*V1212(a,b,j1,k1)*V1212(a,b,i,k)
		    + V1212(a,b,k1,j1)*V1212(a,b,i,k)
		    + V1212(a,b,j1,k1)*V1212(a,b,k,i)) * d(k,k1) /
	    (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k1)-phis_.energy(j1)) 
	       *(phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(i)));
	  
	}
      }
    }
  
  
  d212 *= -.25;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d212;
  
  return d212;
  
}




inline double Second_order_cap::d2_13(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  unsigned i_sym = phis_.irrep(i);
  unsigned j_sym = phis_.irrep(j);
  unsigned i1_sym = phis_.irrep(i1);
  unsigned j1_sym = phis_.irrep(j1);
  

  if (phis_.irrep_product(i_sym,j_sym) != 
      phis_.irrep_product(i1_sym,j1_sym)) return 0.;
  
  
  double d213 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++) {
    unsigned a_sym =  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b), i_sym), j_sym);
    for(unsigned a_ = 0; a_ < virs_[a_sym].size(); a_++) {
      unsigned a = virs_[a_sym][a_];
      unsigned a1_sym = a_sym;
      for(unsigned a1_ = 0; a1_ < virs_[a1_sym].size() ; a1_++) {
	unsigned a1 = virs_[a1_sym][a1_];

	  
	d213 += (  V1212(a1,b,i1,j1)*V1212(a,b,i,j)
		   + V1212(a1,b,j1,i1)*V1212(a,b,j,i)
		   + V1212(a1,b,i1,j1)*V1212(a,b,j,i)
		   + V1212(a1,b,j1,i1)*V1212(a,b,i,j)) * d(a1,a) /
	  (  (phis_.energy(a1)+phis_.energy(b)-phis_.energy(i1)-phis_.energy(j1)) 
	     *(phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j)));
	
	
      }
    }
  
  }

  d213 *= .5;

  if ((i == j) && (i1 == j1)) return .5* d213;
  
  if ((i == j) || (i1 == j1)) return  SQRT_1_2 * d213;
  
  return d213;
  
}



inline double Second_order_cap::d2_14(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  
  
  unsigned i_sym = phis_.irrep(i);
  unsigned j_sym = phis_.irrep(j);
  unsigned i1_sym = phis_.irrep(i1);
  unsigned j1_sym = phis_.irrep(j1);
  

  if (phis_.irrep_product(i_sym,j_sym) != 
      phis_.irrep_product(i1_sym,j1_sym)) return 0.;

  

  
  double d214 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++) {
    unsigned a_sym =  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b), i_sym), j_sym);
    for(unsigned a_ = 0; a_ < virs_[a_sym].size(); a_++) {
      unsigned a = virs_[a_sym][a_];
      unsigned k_sym = j1_sym;
      for(unsigned k_ = 0; k_ < occs_[k_sym].size() ; k_++) {
	unsigned k = occs_[k_sym][k_];
	
	  
 	d214 += (  V1212(a,b,k,i1)*V1212(a,b,j,i)
 		   + V1212(a,b,i1,k)*V1212(a,b,i,j)
 		   + V1212(a,b,k,i1)*V1212(a,b,i,j)
 		   + V1212(a,b,i1,k)*V1212(a,b,j,i)) * d(j1,k) /
 	  (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(i1)-phis_.energy(k)) 
 	     *(phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j)));
       
	
      }
    }
  }
  d214 *= -.25;

  if ((i == j) && (i1 == j1)) return .5 * d214;
  
  if ((i == j) || (i1 == j1)) return SQRT_1_2 * d214;
  
  return d214;
  
}



inline double Second_order_cap::d2_15(unsigned i, unsigned j, unsigned i1, unsigned j1)
{
  


  unsigned i_sym = phis_.irrep(i);
  unsigned j_sym = phis_.irrep(j);
  unsigned i1_sym = phis_.irrep(i1);
  unsigned j1_sym = phis_.irrep(j1);
  

  if (phis_.irrep_product(i_sym,j_sym) != 
      phis_.irrep_product(i1_sym,j1_sym)) return 0.;

  
  
  double d215 = 0.;
  for(unsigned b = phis_.number_occupied(); b < number_orbitals(); b++) {
    unsigned a_sym =  phis_.irrep_product(phis_.irrep_product(phis_.irrep(b), i_sym), j_sym);
    for(unsigned a_ = 0; a_ < virs_[a_sym].size(); a_++) {
      unsigned a = virs_[a_sym][a_];
      unsigned k_sym = i1_sym;
      for(unsigned k_ = 0; k_ < occs_[k_sym].size() ; k_++) {
	unsigned k = occs_[k_sym][k_];
	  
	d215 += (  V1212(a,b,j1,k)*V1212(a,b,j,i)
		   + V1212(a,b,k,j1)*V1212(a,b,i,j)
		   + V1212(a,b,j1,k)*V1212(a,b,i,j)
		   + V1212(a,b,k,j1)*V1212(a,b,j,i)) * d(i1,k) /
	  (  (phis_.energy(a)+phis_.energy(b)-phis_.energy(k)-phis_.energy(j1)) 
	     *(phis_.energy(a)+phis_.energy(b)-phis_.energy(i)-phis_.energy(j)));
	
	
      }
    }
  }
  
  d215 *= -.25;
  
   if ((i == j) && (i1 == j1)) return .5 * d215;
  
   if ((i == j) || (i1 == j1)) return SQRT_1_2 * d215;
  
  return d215;
  
}





inline double Second_order_cap::d_ij_i1j1(unsigned i, unsigned j, unsigned i1, unsigned j1)
{

  double term1 = 0., term2 = 0., term3 = 0., term4 = 0.;
  term1 += d2_1(i,j,i1,j1);
  term2 += d2_2(i,j,i1,j1);
  term3 += d2_3(i,j,i1,j1);
  term4 += d2_4(i,j,i1,j1);
  term1 += d2_5(i,j,i1,j1);
  term2 += d2_6(i,j,i1,j1);
  term3 += d2_7(i,j,i1,j1);
  term4 += d2_8(i,j,i1,j1);
  term1 += d2_9(i,j,i1,j1);
  term2 += d2_10(i,j,i1,j1);
  term3 += d2_11(i,j,i1,j1);
  term4 += d2_12(i,j,i1,j1);
  term1 += d2_13(i,j,i1,j1);
  term2 += d2_14(i,j,i1,j1);
  term3 += d2_15(i,j,i1,j1);


  return (term1 + term2) + (term3 + term4);
  
}


inline double Second_order_cap::vir_sum_rho_d(unsigned i, unsigned i1)
{
  
  unsigned i_sym = phis_.irrep(i);
  if (i_sym != phis_.irrep(i1)) return 0.;
  
  double sum = 0.;
  for(unsigned a_ = 0; a_ < virs_[i_sym].size(); a_++) {
    unsigned a = virs_[i_sym][a_];

    sum += rho->operator()(a,i) * d(i1,a);
  }


  return sum;
}

bool Second_order_cap::block_ii_jj(const Config &row, const Config &col, double &element)
{
    
  unsigned i = row.occ[0];
  unsigned j = col.occ[0];

  bool delta_ij = i == j;

  if (delta_ij) {
    element +=   gs_expt_value2_;
    element -= 2. * (vir_sum_rho_d(i,j)+vir_sum_rho_d(j,i));
  }

  element += d_ij_i1j1(i,i,j,j);
  element += d_ij_i1j1(j,j,i,i);


  return true;
  
}


bool Second_order_cap::block_ij_kk(const Config &row, const Config &col, double &element)
{

  unsigned i = row.occ[0];
  unsigned j = row.occ[1];  
  unsigned k = col.occ[0];

  bool delta_ik = i == k;
  bool delta_jk = j == k;

  if (delta_ik) 
    element -= 2.* SQRT_1_2 * (vir_sum_rho_d(j,k)+vir_sum_rho_d(k,j));

  if (delta_jk) 
    element -= 2.* SQRT_1_2 * (vir_sum_rho_d(i,k)+vir_sum_rho_d(k,i));
 
 
  element += d_ij_i1j1(i,j,k,k);
  element += d_ij_i1j1(k,k,i,j);

  return true;
  
}

bool Second_order_cap::block_ij_kl(const Config &row, const Config &col, double &element)
{

  unsigned i = row.occ[0];
  unsigned j = row.occ[1];
  
  unsigned k = col.occ[0];
  unsigned l = col.occ[1];


  bool delta_ik = i == k;
  bool delta_jk = j == k;
  bool delta_il = i == l;
  bool delta_jl = j == l;


  if (delta_ik && delta_jl)
    element +=   gs_expt_value2_;

  if (delta_ik) 
    element -= vir_sum_rho_d(j,l) + vir_sum_rho_d(l,j);

  if (delta_jk) 
    element -= vir_sum_rho_d(i,l) + vir_sum_rho_d(l,i);

  if (delta_il) 
    element -= vir_sum_rho_d(j,k) + vir_sum_rho_d(k,j);

  if (delta_jl) 
    element -= vir_sum_rho_d(i,k) + vir_sum_rho_d(k,i);


  element += d_ij_i1j1(i,j,k,l);
  element += d_ij_i1j1(k,l,i,j);

  return true;
  
}


