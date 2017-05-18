#ifdef MOLC
#include "phis_molcas.hpp"
#include "blas_matrix.hpp"
#include <string>
#include <algorithm>
#include <limits>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/regex.hpp> 

extern "C" {
#include "/cvos/shared/TC/phis/molcas/74.ifc/include/phis.h"
}

#include <iostream>
using namespace std;

extern "C" void compute_cap_mo_(double*, double*, double*, int*, int*);

void Phis_molcas::molcas2capinput()
{
  ifstream infile ("SCFDATA");
  
  
  if (!infile.is_open())
    throw string("Molcas_scf : \"SCFDATA\" does not exist.");


  // find the lines containing [N_ATOMS], [ATOMS] (AU), [GTO] and [MO]

  stringstream buffer;
  buffer << infile.rdbuf();
  infile.close();   
  string scfdata(buffer.str());

  string::const_iterator start, end;
  start = scfdata.begin(); end = scfdata.end();
  boost::match_results<std::string::const_iterator> what;
  //Regex for the lines after [ATOMS]
  boost::regex expr("[A-Za-z]{1,2}\\d*\\s+\\d+\\s+(\\d+)((?:\\s+[+-]?\\d+\\.\\d+){3})");   
  
  vector<string> elements;
  vector<string> xyz_lines;
  while (regex_search(start, end, what, expr, boost::format_perl)) {
    
    istringstream elem( string (what[1].first,what[1].second) );
    int charge;
    elem >> charge;
    
    switch (charge) {
    case 1: elements.push_back("HYDROGEN"); break;
    case 2: elements.push_back("HELIUM"); break;
    case 3: elements.push_back("LITHIUM"); break;
    case 4: elements.push_back("BERYLLIUM"); break;
    case 5: elements.push_back("BORON"); break;
    case 6: elements.push_back("CARBON"); break;
    case 7: elements.push_back("NITROGEN"); break;
    case 8: elements.push_back("OXYGEN"); break;
    case 9: elements.push_back("FLOURINE"); break;
    case 10: elements.push_back("NEON"); break;
    case 11: elements.push_back("SODIUM"); break;
    case 12: elements.push_back("MAGNESIUM"); break;
    case 13: elements.push_back("ALUMINIUM"); break;
    case 14: elements.push_back("SILICON"); break;
    case 15: elements.push_back("PHOSPHORUS"); break;
    case 16: elements.push_back("SULFUR"); break;
    case 17: elements.push_back("CHLORINE"); break;
    case 18: elements.push_back("ARGON"); break;
    case 19: elements.push_back("POTASSIUM"); break;
    case 20: elements.push_back("CALCIUM"); break;
    case 21: elements.push_back("SCANDIUM"); break;
    case 22: elements.push_back("TITANIUM"); break;
    case 23: elements.push_back("VANADIUM"); break;
    case 24: elements.push_back("CHROMIUM"); break;
    case 25: elements.push_back("MANGANESE"); break;
    case 26: elements.push_back("IRON"); break;
    case 27: elements.push_back("COBALT"); break;
    case 28: elements.push_back("NICKEL"); break;
    case 29: elements.push_back("COPPER"); break;
    case 30: elements.push_back("ZINK"); break;
    case 31: elements.push_back("GALIUM"); break;
    case 32: elements.push_back("GERMANIUM"); break;
    case 33: elements.push_back("ARSENIC"); break;
    case 34: elements.push_back("SELENIUM"); break;
    case 35: elements.push_back("BROMINE"); break;
    case 36: elements.push_back("KRYPTON"); break;
    case 37: elements.push_back("RUBIDIUM"); break;
    case 38: elements.push_back("STRONTIUM"); break;
    case 39: elements.push_back("YTTRIUM"); break;
    case 40: elements.push_back("ZIRCONIUM"); break;
    case 41: elements.push_back("NIOBIUM"); break;
    case 42: elements.push_back("MOLYBDENUM"); break;
    case 43: elements.push_back("TECHNETIUM"); break;
    case 44: elements.push_back("RUTENIUM"); break;
    case 45: elements.push_back("RHODIUM"); break;
    case 46: elements.push_back("PALADIUM"); break;
    case 47: elements.push_back("SILVER"); break;
    case 48: elements.push_back("CADMIUM"); break;
    case 49: elements.push_back("INDIUM"); break;
    case 50: elements.push_back("TIN"); break;
    case 51: elements.push_back("ANTIMONY"); break;
    case 52: elements.push_back("TELLURIUM"); break;
    case 53: elements.push_back("IODINE"); break;
    case 54: elements.push_back("XENON"); break;
    case 55: elements.push_back("CAESIUM"); break;
      
    default:
      elements.push_back("UNKNWON");
    }
    
    istringstream coords(string (what[2].first,what[2].second));
    double coord1, coord2, coord3;
    coords >> coord1 >> coord2 >> coord3;

    ostringstream new_coords;
    new_coords.setf(ios::fixed,ios::floatfield);
    new_coords.precision(8);
    new_coords << "   " << coord1 << "   " << coord2 << "   " << coord3 << endl;

    ostringstream element_name;
    element_name.width(20);
    element_name << left << elements.back();

    xyz_lines.push_back(element_name.str() 
			+ string (what[1].first,what[1].second)
      			+ new_coords.str());
    
    
    start = what[0].second;
  }

  // Create xyz.txt
  
  ofstream xyz ("xyz.txt");
  xyz << "              " <<  xyz_lines.size() << endl;
  for(int i = 0 ; i < xyz_lines.size();i++)
    xyz  << xyz_lines[i];
  xyz.close();
  

  // Create capint.txt
  start = scfdata.begin(); 
  expr = string("(\\s+[spdfSPDF]\\s+\\d+((\\s+[-]?\\d+\\.\\d+E[+-]\\d{2}){2})+)+");   
  
  ofstream capint("capint.txt");
  int count = 0;
  int gto_count = 1;
  while (regex_search(start, end, what, expr, boost::format_perl)) {
    capint << endl;
    start = what[0].second;
    string gto(string(what[0].first, what[0].second));

    gto = boost::regex_replace(gto,boost::regex(string("([sdpf])")),
                               "\\u\\1", boost::format_perl);


    string gto_new;
    int done = 1;
    while(done) {
      ostringstream oss;
      oss << gto_count;
      gto_new = boost::regex_replace(gto,
                                     boost::regex(string("(\\s{2,})([-]?\\d+\\.\\d+E[+-]\\d{2}\\s+[-]?\\d+\\.\\d+E[+-]\\d{2}$)")),
                                     string(string("")+ "\\1"+ oss.str()+" \\2").c_str(),
                                     boost::format_first_only | boost::format_perl);

      done=gto_new.compare(gto);
      if (done) gto_count++;
      gto = gto_new;

    }

    capint << elements[count++]
	   << gto << "\nE 1 2\n";
  }

  capint.close();

  /* Read the active space, frozen and deleted orbitals from OUTPUT
     Three arrays - frozen, deleted and usedorb are created. The total
     number of orbitals (a sum of the three) is used to determine the
     "bounds" (vector bounds) of the irreducible representations in
     the matrix with MO coefficients - and, therefore, delete and freeze
     orbitals*/

  ifstream infile11("OUTPUT");
  if (!infile11.is_open())
    throw string("Molcas_scf : \"OUTPUT\" does not exist.");

  buffer.str("");buffer.clear();
  buffer << infile11.rdbuf();
  infile11.close();
  string output(buffer.str());

  start = output.begin();end=output.end();
  boost::regex expr1("Orbital specifications:\\s+\\-+\\s+.+?"
          "Frozen orbitals:\\s+((:?\\d+\\s+)+\\d)"
          "\\s+Deleted orbitals:\\s+((:?\\d+\\s+)+\\d)"
          "\\s+Number of orbitals used:\\s+((:?\\d+\\s+)+\\d)");
  unsigned int frozen[] = {0,0,0,0,0,0,0,0};
  unsigned int deleted[] = {0,0,0,0,0,0,0,0};
  unsigned int usedorb[] = {0,0,0,0,0,0,0,0};
  while (regex_search(start, end, what, expr1, boost::format_perl)) {
      start = what[0].second;
      istringstream frzn(string (what[1].first, what[1].second));
      int i = 0;
      while (!frzn.eof())
          frzn >> frozen[i++];
      istringstream del(string (what[3].first, what[3].second));
      i = 0;
      while (!del.eof())
          del >> deleted[i++];

      istringstream used(string (what[5].first, what[5].second));
      i = 0;
      while (!used.eof())
          used >> usedorb[i++];
  }

  vector<int> bounds;
  bounds.push_back(0);
  int sum_index = 0;
  for (int i = 0; i < 8; ++i)
  {
    sum_index += frozen[i] + deleted[i] + usedorb[i];
    bounds.push_back(sum_index);
  }


  // Create mocoef.txt
  start = scfdata.begin();end=scfdata.end(); 
  expr = string("Sym=\\s+\\d+\\w+.+?Occup=\\s+\\d\\.\\d+\\s+(\\d+\\s+\\-?\\d\\.\\d+\\s+)+");   
  vector<vector <double> > scfmat;
  while (regex_search(start, end, what, expr, boost::format_perl)) {
    start = what[0].second;
    string scfvec(what[0].first, what[0].second);

    vector <double> scfv;
    boost::match_results<std::string::const_iterator> what1;
    string::const_iterator vec_start = scfvec.begin(), vec_end =scfvec.end();
    while (regex_search(vec_start, vec_end, what1, 
			boost::regex(string("\\d+\\s+(\\-?\\d\\.\\d+)")),
			boost::format_perl)) {
      vec_start = what1[0].second;
      istringstream iss(string(what1[1].first,what1[1].second));
      double comp;
      iss >> comp;
      scfv.push_back(comp);
    }
    // Stores the scf vectors in scfmat
    scfmat.push_back(scfv);
  }
  
  // Delete the frozen and virtual orbitals as specified by MOTRA
  for (int i  = 0; i < 8; i++) {
    int num_deleted = 0;
    while (frozen[i]--) {
       scfmat[bounds[i]+num_deleted++].clear();
    }
    num_deleted = 0;
    while (deleted[i]--) {
       scfmat[bounds[i+1] - num_deleted++].clear();
    }
  }
  // Purge the empty scf vectors
  vector<vector<double> >::iterator iter = scfmat.begin();
  while ( iter !=  scfmat.end()) {
    if (iter->empty()) {
       iter = scfmat.erase(iter);
     } else {
       iter++;
     }
  }
  




/*The SCF vectors in Molcas are ordered by symmetry
  We need to reoreder them by energy.
  This reordering has been done already in get_scfvec.
  I simply copy the code here. Lazy, lazy. */
//  dim1 = nAO                , dim2 = nMO
  int dim1 = scfmat[0].size(), dim2 = scfmat.size();
  double* mat = new double[dim1 * dim2];
  for (int j = 0; j < dim2; j++) 
    for(int i = 0; i < dim1; i++) {
          mat[i + j*dim1] =  scfmat[j][i];
    } 
  
  double* temp_mat = new double[dim1*dim2];
  
  for (int j = 0; j < dim2; j++) 
    for (int i = 0; i < dim1; i++) {
      temp_mat[i +j * dim1] = mat[i + j*dim1];
    }
  
  // Reordering by energy
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
      mat[i+j*dim1] = temp_mat[i + ind[j]*dim1] ;
    } 
  }
    
  delete [] temp_mat;

   ofstream mocoef("mocoef.txt");
   mocoef << dim1 << ' ' << dim2 << endl;
   for(int j = 0; j < dim2; j ++)
     for(int i = 0; i < dim1; i ++)
       mocoef << j+1 << " " << i+1 << " " << mat[i+j*dim1] << endl;
   mocoef.close();


}




Phis_molcas::Phis_molcas(const char *inp_backend)
{
  string backend_name(inp_backend);
  int backend;
  
  /* choose backend for PHIS */
  if(!backend_name.compare("molcas"))
    backend=MOLCAS | (SYM_BLOCKED << 16);
  else if(!backend_name.compare("guk"))
    backend=GUK;
  else 
    throw string("SCF_data_reader::get_info(): Uknown Backend \
Available: guk, molcas.\n");

  int cap;
  /* initialize PHIS */
  cap=phis_init(&backend);
  
  /* check if backend has the needed functions */
  if(~(cap | ~(HAVE_INFO | HAVE_LOA | HAVE_EPSI | HAVE_SYM | HAVE_OCC))) {
    
    throw string("SCF_data_reader::get_info: Backend incomplete. \
Terminating.\n");
  }

  sort();
  
}


// Produces two array (ind and inv_ind) members used in the reordering
// of the MO (from symmetry ordered to energy ordered)
void Phis_molcas::sort()
{


  unsigned number_irreps_, number_orbitals_;
  int atoms;
  get_info((int*) &number_irreps_, (int*) &number_orbitals_, &atoms);

  

  int test = number_orbitals_;
  // Obtain orbital energy information
  
  double *epsi = new double[number_orbitals_];
  double e_hf;
  vector<double> orbital_energy_;

  phis_get_epsi(&e_hf, epsi, &test);
  if(test != number_orbitals_)
    throw string("get_epsi(): Not enough space to save the orbital energies\n");
  
  orbital_energy_.assign(epsi, epsi + number_orbitals_);
  delete [] epsi;



  vector<double> temp_e(orbital_energy_);
  for (int i = 0; i < number_orbitals_; i++) {
    vector<double>::iterator min_e 
      = min_element(temp_e.begin(), temp_e.end());
    ind.push_back(min_e - temp_e.begin());
    *min_e = numeric_limits<double>::max();
  }


  inv_ind.assign(number_orbitals_,0);
  for (int i = 0; i < number_orbitals_; i++) 
    inv_ind[ind[i]] = i;

}



void Phis_molcas::get_info (int *nSym, int *nBas, int *nCenters)
{
  phis_get_info(nSym, nBas, nCenters);
}
void Phis_molcas::get_epsi (double *E_hf, double *e, int *n)
{
  phis_get_epsi (E_hf, e, n);
  if (*n > 0)
    reorder<double>(e,*n,1);
}
void Phis_molcas::get_sym (int *s, int *n)
{
  phis_get_sym (s, n);
  if (*n > 0)
    reorder<int>(s,*n,1);

}
void Phis_molcas::get_occ (double *o, int *n)
{
  phis_get_occ (o, n);
  if (*n > 0)
    reorder<double>(o,*n,1);

}
void Phis_molcas::get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint)
{
  phis_get_next_Vpqrs(p, q, r, s, vint);


  if (*s == 0) return;
  //  cout << *p-1 << ' ' << *q-1 << ' ' << *r-1 << ' ' << *s-1 << ' ' << *vint<< endl;
  *p = inv_ind[*p-1]+1;
  *q = inv_ind[*q-1]+1;
  *r = inv_ind[*r-1]+1;
  *s = inv_ind[*s-1]+1;

  //  cout << *p << ' ' << *q << ' ' << *r << ' ' << *s << ' ' << *vint<< endl;

}

void Phis_molcas::get_scfvec (double *C, int *n, int *len)
{
  //THIS CALL PRODUCES SEGMENTATION FAULT!!!
  phis_get_scfvec (C, n, len);
  // Here the symmetry ordering holds only for mo (coulumns)




  if ((*n > 0) && (*len > 0)) {
    double *mat = C;
    int dim1 = *len, dim2 = *n;
    double* temp_mat = new double[dim1*dim2];
   
    
    for (int j = 0; j < dim2; j++) 
      for (int i = 0; i < dim1; i++) {
	temp_mat[i +j * dim1] = mat[i + j*dim1];
      }
  
    // Reordering by energy
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	mat[i+j*dim1] = temp_mat[i + ind[j]*dim1] ;
      } 
    }
    
    delete [] temp_mat;

  }
}

// Seems to have incorrect behavior 
void Phis_molcas::get_overlap (double *S, int *n)
{
  phis_get_overlap (S, n);
}

void Phis_molcas::get_dip (double *x, double *y, double *z, int *n)
{
  phis_get_dip (x, y, z, n);
  if (*n > 0) {
    reorder<double>(x,*n,*n);
    reorder<double>(y,*n,*n);
    reorder<double>(z,*n,*n);
  }

}
void Phis_molcas::get_vel (double *x, double *y, double *z, int *n)
{
  phis_get_vel (x, y, z, n);
  if (*n > 0) {
    reorder<double>(x,*n,*n);
    reorder<double>(y,*n,*n);
    reorder<double>(z,*n,*n);
  }
}
void Phis_molcas::get_quad(double *xx, double *xy, double *xz,
			     double *yy, double *yz, double *zz, int *n)
{
  phis_get_quad(xx, xy, xz, yy, yz, zz, n);
  if (*n > 0 ) {
    reorder<double>(xx,*n,*n);
    reorder<double>(xy,*n,*n);
    reorder<double>(xz,*n,*n);
    
    reorder<double>(yy,*n,*n);
    reorder<double>(yz,*n,*n);
    
    reorder<double>(zz,*n,*n);
  }
}
void Phis_molcas::get_oneel(double *h, double *t, int *n)
{
  phis_get_oneel(h, t, n);
  if (*n > 0) {
    reorder<double>(h,*n,*n);
    reorder<double>(t,*n,*n);
  }
}


void Phis_molcas::get_geometry (int *n, double *geometry, double *Z_nuc, double *E_nuc)
{
  phis_get_geometry(n, geometry, Z_nuc, E_nuc);
}


void  Phis_molcas::get_mo_cap (double *boxx, double *boxy, double *boxz, int* nmo, double* capmo) 
{

  molcas2capinput();

  // Source copied from gus_reader.cpp
  static int ISPC=3000;
  
  // !    open(unit=8,file='mocoef.txt',STATUS='UNKNOWN',FORM='FORMATTED')
  // !    read (8, *) NB
  // !    close(8)
  // Get the number of basis functions as done above
  ifstream mocoef("mocoef.txt");
  
  if (!mocoef.is_open())
    throw string("GUS_reader::get_mo_cap : \"mocoef.txt\" does not exist.");
  int nb;
  mocoef >> nb;
  mocoef.close();
  
  //     !COMPUTES CAP MATRIX ELEMENTS IN MO BASIS
  //     !AND PRINTS THEM IN capmat.dat
  // !    call Compute_CAP_MO(NB,ISPC)

  // Creates a capmat.txt file containing the CAP matrix
  compute_cap_mo_(boxx,boxy,boxz,&nb,&ISPC);
 
  ifstream capmat("capmat.dat");
  if (!capmat.is_open())
    throw string("Phis_molcas::get_mo_cap : \"capmat.txt\" does not exist.");
 
  int dim;
  string line;
  istringstream is;
  getline (capmat,line);
  is.str(line);
  is >> dim;
  // Check the input parameter
  if (dim != *nmo) {
    *nmo = -dim;
    return;
  }
  is.clear();
  //Read data 

  // Triangular matrix
  unsigned row = 0, col = 0;
  while (row < dim ) {
    getline (capmat,line);
    is.str(line);
    is >> capmo[row+dim*col];
    capmo[col+dim*row] = capmo[row+dim*col];
    is.clear();
    // Its a triangular matrix
    if (row == col) {
      col = 0;
      row++;
    } else
	col++;
  }
  capmat.close();
  //all elements of occupied orbitals must be zero
  //get number of occupied first
  double* occ = new double[dim];
  get_occ(occ,&dim);
  int nocc = 0;
  for (int  i = 0; i < dim; i++)
    if(occ[i]) nocc++;
  delete [] occ;
  
  for(int i = 0; i < nocc; i++)
    for(int j = 0; j < dim; j++) 
     capmo[i+dim*j] = capmo[j+dim*i] =0.;
}




#endif
