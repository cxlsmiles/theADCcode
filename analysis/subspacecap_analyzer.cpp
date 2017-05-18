#include "subspacecap_analyzer.hpp"
#include "ndadc3_ip/nd_adc3_matrix.hpp"
#include "adc2_dip/adc2_matrix.hpp"
#include "adc2_dip/singlet.hpp"
#include "adc2_dip/triplet.hpp"
#include "blas_matrix.hpp"
#include "ndadc3_ip/self_energy.hpp"
#include "scf_data/scf_data_reader.hpp"
#include "input_struct.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace std;

#if defined INTEL
#define DOCAP     main_cap_mp_docap_
#elif defined PGI
#define DOCAP     main_cap_docap_
#elif defined GCC
#define DOCAP     __main_cap_MOD_docap
#endif

extern Integral_table* integral_table;

extern "C" {
  // Input parameters used by DOCAP
  static struct _Info { 
    int nmo,nocc, nicap, ioffset, incr, nholes, mult; 
    double srcap, sicap, eadcmax, eadcmin;
    char output[64];
  } info_;
  void DOCAP(struct _Info*, double*, double*, double*);
}


RSCAP_analyzer::RSCAP_analyzer(SCF_data_reader& phis, struct Diag_info& diag, 
			       struct Eigen_info& eig, struct Cap_info& cap)
  : phis_(phis), ADC_analyzer(diag, eig), cap_(cap) 
{

  if (!cap_.sicap) throw string("Subspace_CAP_analyzer:: step size of CAP not selected.");
  if (!cap_.nicap) throw string("Subspace_CAP_analyzer:: iteration number of CAP not selected.");
  if (cap_.ioffset < 0) throw string("Subspace_CAP_analyzer:: iteration offset not a non-negative integer.");
  if (!cap_.boxx || !cap_.boxy || !cap_.boxz) throw string("Subspace_CAP_analyzer:: box size of CAP not selected.");

}



// This method prepares the input parameters for Nicolas'
// FORTRAN routine that builds the subspace CAP trajectories
int RSCAP_analyzer::analyze( ADC_matrix& adc_mat)
{


     cout << "IN RSCAP" << endl;
  // Initialize the structure that is to be passed to CAP
  int nmo = info_.nmo = phis_.number_orbitals();
  info_.nocc = phis_.number_occupied();
  info_.nicap = cap_.nicap;
  info_.ioffset = cap_.ioffset;
  info_.incr = cap_.increase;  
  info_.srcap = cap_.srcap;  info_.sicap = cap_.sicap;  
  info_.eadcmax = cap_.eadcmax;  info_.eadcmin = cap_.eadcmin;
  for (int i = 0; i < 64; i++) info_.output[i] = cap_.output[i];

  // Get the density matrix
  Blas_matrix dens(nmo,nmo);
  dens = 0.;
// JUST FOR TESTING
  Self_energy(*integral_table,phis_).density_3plus(dens);
  
  // Get the orbital energies
  double* engs = new double[phis_.number_orbitals()];
  for (unsigned i = 0; i < phis_.number_orbitals(); i++) engs[i] = phis_.energy(i);
  
  Blas_matrix cap(nmo, nmo);
  phis_.get_mo_cap(&cap_.boxx, &cap_.boxy, &cap_.boxz, &nmo, &cap(0,0));
  cap.print();
  if (nmo < 0)
	    throw string (" RSCAP_analyzer:: error reading CAP matrix.");
  
  // Determine what type of vectors  we need IP or DIP:
  Adc2_matrix* is_adc2 = dynamic_cast<Adc2_matrix*>(&adc_mat);
  ND_ADC3_matrix* is_ndadc3 = dynamic_cast<ND_ADC3_matrix*>(&adc_mat);
  if (is_ndadc3){
    print_nicolas_input_ip(adc_mat);
    info_.nholes = 1; // Ionization
    info_.mult = 2;
  } else if (is_adc2){

    print_nicolas_input_dip(adc_mat);
    info_.nholes = 2; // Double ionization
    
//    Singlet* is_singlet = dynamic_cast<Singlet*>(is_adc2);
//    Triplet* is_triplet = dynamic_cast<Triplet*>(is_adc2);


//    if (is_singlet) 
//      info_.mult = 1;
//    else if (is_triplet)
//      info_.mult = 3;

      info_.mult = adc_mat.spin()+1;
    
  } else
    throw string("Unknown propagator\n.");

  
  DOCAP(&info_, &dens(0,0),engs,&cap(0,0));  
  
   // Nico 24.02.2011
   // send dipMO to docap
   //  int n =  0;
  //   n= phis_.number_orbitals();
  //   Blas_matrix x(n,n), y(n,n), z(n,n);
  //   phis_.get_dip(&x(0,0), &y(0,0), &z(0,0), &n);
  //   z.print();
  //   DOCAP(&info_, &dens(0,0),engs,&z(0,0));  
  // // end Nico part
  
  delete [] engs;

  return 0;
    
}


void RSCAP_analyzer::print_nicolas_input_ip( ADC_matrix& adc_mat)
{



  ofstream capdump("capdumpfile.dat", ofstream::out | ofstream::trunc);

  Blas_matrix eig_vecs, eig_vals, resnoms;

  diagonalize(adc_mat, eig_vecs, eig_vals, resnoms);



  eig_vals *= 27.211396;
  capdump << eig_vecs.cols() << " !nb of ADC eigenvecs" << endl;

  

  for(int i = 0; i < eig_vecs.cols(); i++) {

    capdump << endl << eig_vals(i) << " !eigenval" << endl;

    
    int comp_count = 0;
    ostringstream vec_comps;

    for(int j = 0; j < eig_vecs.rows(); j++) {
      
      // // printing threshold
      //if (fabs(eig_vecs(j,i)) < 10e-10) continue;
      comp_count++;
    
      if ( j < adc_mat.main_block_size())  
	vec_comps << 1 << ' ' << 0 << ' ';
      else 
	vec_comps << 2 << ' ' << 1 << ' ';
      
      
      //capdump << adc_mat.get_conf(j) << ' ';
      // change the format of the configurations to make Nicolas' life easier
      {
	ostringstream new_conf;
	istringstream conf(adc_mat.get_conf(j));
	char symb = 0;
	int a;
	//       <         num        
	conf >> symb; conf >> a; 
	// , or |
	conf >> symb;
	if (!conf.eof() && symb == '|') new_conf << "<" << a << "| 0 ";
	else {
	  int b,c;
	  //     num         ,        num  
	  conf >> b;conf >> symb; conf >> c;
	  conf >> symb;
	  if (!conf.eof() && symb == '|') new_conf << "<" << a << "," << b << "," << c << "| 1 ";

	  // ,
	  conf >> symb;

	  conf>> symb;
	  if (!conf.eof() && symb == '|') new_conf << "<" << a << "," << b << "," << c << "| 21 ";

	  conf >> symb;
	  if (!conf.eof() && symb == '|') new_conf << "<" << a << "," << b << "," << c << "| 22 ";


	}
	
	vec_comps << new_conf.str();
      }// end
      
      
      
      vec_comps << eig_vecs(j,i) << endl;
      
      
      
    }
    
    capdump << comp_count << endl;
    capdump << vec_comps.str();
  }

  capdump.close();

}



void RSCAP_analyzer::print_nicolas_input_dip( ADC_matrix& adc_mat)
{


  ofstream capdump("capdumpfile.dat", ofstream::out | ofstream::trunc);

  Blas_matrix eig_vecs, eig_vals, renoms;

  diagonalize(adc_mat, eig_vecs, eig_vals, renoms);



  eig_vals *= 27.211396;
  capdump << eig_vecs.cols() << " !nb of ADC eigenvecs" << endl;



  for(int i = 0; i < eig_vecs.cols(); i++) {

    capdump << endl << eig_vals(i) << " !eigenval" << endl;
    //capdump << eig_vecs.rows() << endl;
    
    int comp_count = 0;
    ostringstream vec_comps;

    for(int j = 0; j < eig_vecs.rows(); j++) {

      //      if (fabs(eig_vecs(j,i)) < 10e-10) continue;
      comp_count++;
    
      if ( j < adc_mat.main_block_size()) 
	vec_comps << 2 << ' ' << 0 << ' ';
      else 
	vec_comps << 3 << ' ' << 1 << ' ';
      
      
      //capdump << adc_mat.get_conf(j) << ' ';
      // change the format of the configurations to make Nicolas' life easier
      {
	ostringstream new_conf;
	istringstream conf(adc_mat.get_conf(j));
	char symb;
	int a,b;
	//       <         num        ,      num
	conf >> symb; conf >> a; conf >> symb; conf >> b;
	// , or |
	conf >> symb;
	if (symb == '|') new_conf << "<" << a << "," << b << "| 0 ";
	else {
	  int c,d;
	  //     num    ,             num
	  conf >> c; conf >> symb; conf >> d;
	  conf >> symb;
	  if (symb == '|' && !conf.eof()) new_conf << "<" << a << "," << b << "," << c << "," << d << "| 1 ";

	  // ,
	  conf >> symb;

	  conf>> symb;
	  if (symb == '|' && !conf.eof()) new_conf << "<" << a << "," << b << "," << c << "," << d << "| 21 ";

	  conf >> symb;
	  if (symb == '|' && !conf.eof()) new_conf << "<" << a << "," << b << "," << c << "," << d << "| 22 ";

	  conf >> symb;
	  if (symb == '|' && !conf.eof()) new_conf << "<" << a << "," << b << "," << c << "," << d << "| 23 ";
	  

	}
	
	vec_comps << new_conf.str();
      }// end
      

      
      vec_comps << eig_vecs(j,i) << endl;
      
      
      
    }

    capdump << comp_count << endl;
    capdump << vec_comps.str();
    
  }


  capdump.close();

}
 
