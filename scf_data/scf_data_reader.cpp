#include "scf_data_reader.hpp"
//#include "phis_guk.hpp"
//#include "phis_molcas.hpp"
#include "gus_reader.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;


#define SCF_THRESH 1e-6


const unsigned int SCF_data_reader::mult_table_[][8] = 
  {
    {0, 1, 2, 3, 4, 5, 6, 7},
    {1, 0, 3, 2, 5, 4, 7, 6},
    {2, 3, 0, 1, 6, 7, 4, 5},
    {3, 2, 1, 0, 7, 6, 5, 4},
    {4, 5, 6, 7, 0, 1, 2, 3},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {6, 7, 4, 5, 2, 3, 0, 1},
    {7, 6, 5, 4, 3, 2, 1, 0}
  };

// Similar to vpqrs_, the global pointer and the Multab_ function here
// are to be used by the Fortran and C parts of the code


extern "C" int Multab_(int a, int b)
{
  return SCF_data_reader::irrep_product(a-1, b-1) + 1;
}


#include <iostream>
// Prerequisite: GAMESS-UK SCF calculation
SCF_data_reader::SCF_data_reader(const char *inp_backend)
{
  
  string backend_name(inp_backend);
  bool backend_found  = false;
  string backend;
//#ifdef GUK
//  backend.append("guk ");
//  if(!backend_name.compare("guk")) {
//    scf = new Phis_guk(inp_backend);
//    backend_name = "Gamess UK";
//    backend_found = true;
//  } 
//#endif
//#ifdef MOLC
//  backend.append("molcas ");
//  if(!backend_name.compare("molcas") ) {
//    scf = new Phis_molcas(inp_backend);
//    backend_name = "Molcas";
//    backend_found = true;
//  } 
//#endif
#ifdef GUS
  backend.append("gus ");
  if(!backend_name.compare("gus")) {
    scf = new GUS_reader;
    backend_name = "Gamess US";
    backend_found = true;
  } 
#endif

  if ( !backend_found )
    throw string("SCF_data_reader::SCF_data_reader(): Uknown Backend, Available: "+ backend + "\n");

  get_info(inp_backend);
  get_sym();
  get_epsi();
  get_occ();

  cout << "  SCF data obtained from " << backend_name << '.' << endl;
  cout << "  m.o.\tsym\tenergy(a.u.)\n";
  cout << "  -----------------------------\n";
  for(int i = 0; i < number_orbitals(); i++)
        cout <<"  "<< i+1 << "\t" << irrep(i)+1 << "\t" << energy(i) << endl;
  cout << endl;

}


SCF_data_reader::~SCF_data_reader()
{
  
  delete scf;

}


void SCF_data_reader::get_info (int *nSym, int *nBas, int *nCenters)
{
  scf->get_info(nSym, nBas, nCenters);
}

void SCF_data_reader::get_epsi (double *E_hf, double *e, int *n)
{
  scf->get_epsi(E_hf, e, n);
}

void SCF_data_reader::get_sym (int *s, int *n)
{
  scf->get_sym(s, n);
}

void SCF_data_reader::get_occ (double *o, int *n)
{
  scf->get_occ(o, n);
}

void SCF_data_reader::get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint)
{
  scf->get_next_Vpqrs(p, q, r, s, vint);
}

void SCF_data_reader::get_scfvec (double *C, int *n, int *len)
{
  scf->get_scfvec(C, n, len);
}

void SCF_data_reader::get_overlap (double *S, int *n)
{
  scf->get_overlap(S, n);
}

void SCF_data_reader::get_dip (double *x, double *y, double *z, int *n)
{
  scf->get_dip (x, y, z, n);
}
void SCF_data_reader::get_vel (double *x, double *y, double *z, int *n)
{
  scf->get_vel (x, y, z, n);
}
void SCF_data_reader::get_quad(double *xx, double *xy, double *xz,
	      double *yy, double *yz, double *zz, int *n)
{
  scf->get_quad(xx, xy, xz, yy, yz, zz, n);
}
void SCF_data_reader::get_oneel(double *h, double *t, int *n)
{
  scf->get_oneel(h, t, n);
}

void SCF_data_reader::get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc)
{
  scf->get_geometry(n, geometry, Z_nuc, E_nuc);
}


void SCF_data_reader::get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo)
{
  scf->get_mo_cap (boxx, boxy, boxz, nmo, capmo);
}

// Initializes the phis backend and obtains the number of symmetries, 
// and the number of orbitals
inline void SCF_data_reader::get_info(const char *inp_backend) {
  
  int atoms;
  get_info((int*) &number_irreps_, (int*) &number_orbitals_, &atoms);
  
}

// Obtains the orbital symmetries
inline void SCF_data_reader::get_sym()
{ 
  
  int *syms = new int[number_orbitals_];
  int test = number_orbitals_;
  
  get_sym(syms, &test);
  if(test != number_orbitals_) 
    throw string("SCF_data_reader::get_sym(): Not enough space\
                         to save  the symmetry labels.");
  
  orbital_symmetry_.assign(syms, syms + number_orbitals_);
  delete [] syms;


  // The numbering of the irreps in "phis" starts from 1. Shifts to 0
  // and checks if there is an orbital without symmetry.
  for(vector<int>::iterator it = orbital_symmetry_.begin(); 
      it != orbital_symmetry_.end(); ++it) 
    if (0 == (*it)--) 
      throw string("get_sym(): An orbital has no symmetry!\n(Don't use \"bypass hf\" in integral transformation)\n");

}



// Obtains the orbital energies
inline void SCF_data_reader::get_epsi()
{
  
  int test = number_orbitals_;
  // Obtain orbital energy information
  
  double *epsi = new double[number_orbitals_];
  double e_hf;
  
  get_epsi(&e_hf, epsi, &test);
  if(test != number_orbitals_)
    throw string("get_epsi(): Not enough space to save the orbital energies\n");
  
  orbital_energy_.assign(epsi, epsi + number_orbitals_);
  delete [] epsi;

}


// Determins the number of occupied orbitals
// from the fetched occupation numbers 
inline void SCF_data_reader::get_occ()
{
  
  int test = number_orbitals_;
  // Obtain orbital occupancy 
  double *occ = new double [number_orbitals_];  
  get_occ(occ,&test);
  if(test != number_orbitals_) 
    throw string("read_scf_data: Not enough space to save the occupation numbers\n");
  
  // phis/molcas appears to return the orbitals ordered by symmetry
  // rather than energy (unlike phis/guk). If, however, we sort the arrays here
  // the code will work for both cases.
  //sort(occ);

  unsigned int i = 0;
  while (fabs(occ[i] - 2.0) < SCF_THRESH) i++;
  number_occupied_ = i;
  for(; i < number_orbitals_; i++) 
    if (occ[i] > SCF_THRESH)
      throw string("SCF_data_reader::read_scf_data: Occupied orbital after HOMO; Not a ground-state closed-shell calculation.");
  
  delete [] occ;
  
}

// Extracts the next integral from the guk integral file
double SCF_data_reader::get_integral(int& p, int& q, int& r, int& s)
{
  double integral;
  get_next_Vpqrs(&p, &q, &r, &s, &integral);
  // Numbering shifted by -1
  p--;q--;r--;s--;

  //  cout << p << ' ' << q << ' ' << r << ' ' << s << ' ' << integral<< endl;

  return integral;
}





