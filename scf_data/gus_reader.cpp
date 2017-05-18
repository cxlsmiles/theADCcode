#ifdef GUS
#include "gus_reader.hpp"
#include <boost/regex.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include <fstream>
#include <vector>
#include "blas_matrix.hpp"

using namespace std;


#define RECORD_SIZE 791

// Abelian group numbering in Gamess-US:
static const char* irrep_table[][8] =
  { {"A" , ""   , ""   , ""   , ""  , ""   , ""   , ""   }, //C1
    {"Ag", "Au" , ""   , ""   , ""  , ""   , ""   , ""   }, //Ci
    {"A'", "A''", ""   , ""   , ""  , ""   , ""   , ""   }, //Cs
    {"A" , "B"  , ""   , ""   , ""  , ""   , ""   , ""   }, //C2
    {"A1", "A2" , "B1" , "B2" , ""  , ""   , ""   , ""   }, //C2v
    {"Ag", "Bg" , "Au" , "Bu" , ""  , ""   , ""   , ""   }, //C2h
    {"A" , "B1" , "B2" , "B3" , ""  , ""   , ""   , ""   }, //D2
    //{"Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"}  //D2h
    {"Ag", "B2u", "B3u", "B1g", "B1u", "B3g", "B2g", "Au" }  //D2h
  };



void GUS_reader::get_info (int *nSym, int *nBas, int *nCenters)
{
  int& sym  = *nSym;
  int& bas = *nBas;
  int& cent  = *nCenters;
  
  ifstream orbfile("orbital.txt");

  if (!orbfile.is_open())
    throw string("GUS_reader::get_info : \"orbital.txt\" does not exist.");

  ostringstream oss;  
  string input;
  while (getline(orbfile, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");


  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 
  istringstream iss;
  

  if (regex_search(start, end, what,
		   boost::regex("^\\s+THE POINT GROUP IS [\\w']+\\s*, NAXIS= \\d, ORDER= (\\d)$"), 
		   boost::format_perl)) {
    
    iss.str(string(what[1].first, what[1].second));
    iss >> sym;
    iss.clear();
      
    
  } else 
    throw string("GUS_reader::get_info : Could not read the number of irreps.");

  start = input.begin(); end = input.end();
  if (regex_search(start, end, what,
		   boost::regex("^ TOTAL NUMBER OF ATOMS                        =\\s+(\\d+)$"), 
		   boost::format_perl)) {

    iss.str(string(what[1].first, what[1].second));
    iss >> cent;
    iss.clear();
  } else 
    throw string("GUS_reader::get_info : Could not read the number of atoms.");

  
  bas = 0;
  start = input.begin(); end = input.end();
  while (regex_search(start, end, what,
		      boost::regex("^\\s*\\d+\\s+[\\w']+\\s+\\-?\\d+\\.\\d+$"), 
		      boost::format_perl)) {
    bas++;
    
    start = what[0].second; 
  }

  if (!bas) throw string("GUS_reader::get_info : Could not find MO orbitals.");

  orbfile.close();

}
void GUS_reader::get_epsi (double *E_hf, double *e, int *n)
{

  ifstream orbfile("orbital.txt");
  
  if (!orbfile.is_open())
    throw string("GUS_reader::get_info : \"orbital.txt\" does not exist.");
  
  ostringstream oss;  
  string input;

  while (getline(orbfile, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");


  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 
  istringstream iss;

  
  int bas = 0;
  while (regex_search(start, end, what,
		      boost::regex("^\\s*\\d+\\s+[\\w']+\\s+(\\-?\\d+\\.\\d+)$"), 
		      boost::format_perl)) {
   
    start = what[0].second; 
    iss.str(string(what[1].first, what[1].second));

    if (bas  < *n)
      iss >> e[bas];
    
    bas++;
    iss.clear();

  }

  if (!bas) throw string("GUS_reader::get_info : Could not find MO orbitals.");
  if (bas > *n) *n = -bas;
  else *n = bas;
  
 

  orbfile.close();
  
  
}
void GUS_reader::get_sym (int *s, int *n)
{

  ifstream orbfile("orbital.txt");
  
  if (!orbfile.is_open())
    throw string("GUS_reader::get_info : \"orbital.txt\" does not exist.");
  
  ostringstream oss;  
  string input;

  while (getline(orbfile, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");



  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 

  istringstream iss;
  string sym;
  int naxis;


  if (regex_search(start, end, what,
		   boost::regex("^\\s+THE POINT GROUP IS ([\\w']+)\\s*, NAXIS= (\\d), ORDER= \\d$"), 
		   boost::format_perl)) {
    
    sym.assign(what[1].first, what[1].second);
    iss.str(string(what[2].first, what[2].second));
    iss >> naxis;
    
  } else 
    throw string("GUS_reader::get_sym : Could not read the symmetry of the molecule.");


  int sym_num;
  if((naxis == 1) && regex_match(sym, boost::regex("Cn", boost::regex::icase)))
    
    sym_num = 0;
  
  else if(regex_match(sym, boost::regex("Ci", boost::regex::icase)))
    
    sym_num = 1;
  
  else if(regex_match(sym, boost::regex("Cs", boost::regex::icase)))
    
    sym_num = 2;
  
  else if((naxis == 2) && regex_match(sym, boost::regex("Cn", boost::regex::icase)))
    
    sym_num = 3;
  
  else  if(regex_match(sym, boost::regex("Cnv", boost::regex::icase)))
    
    sym_num = 4;
  
  else if(regex_match(sym, boost::regex("Cnh", boost::regex::icase)))
    
    sym_num = 5;
  
  else if(regex_match(sym, boost::regex("Dn", boost::regex::icase)))
    
    sym_num = 6;
  
  else if(regex_match(sym, boost::regex("Dnh", boost::regex::icase)))
    
    sym_num = 7;
  
  else
    
    throw string("GUS_reader::get_sym : Unknown symmetry label.");
  
  
  
  int bas = 0;
  start = input.begin(); end = input.end();
  while (regex_search(start, end, what,
		      boost::regex("^\\s*\\d+\\s+([\\w']+)\\s+\\-?\\d+\\.\\d+$"), 
		      boost::format_perl)) {
   
    start = what[0].second; 
    sym.assign(what[1].first, what[1].second);
    
    
    int i;
    if (bas  < *n) 
      for(i = 0; i < 8; i++) 
	if(regex_match(sym, boost::regex(irrep_table[sym_num][i], boost::regex::icase))) {
	  s[bas] = i+1;
	  break;
	}
      
    if (i == 8) throw string("GUS_reader::get_sym : Uknown irrep label.");
    bas++;
  }
  
  if (!bas) throw string("GUS_reader::get_info : Could not find MO orbitals.");
  

  if (bas > *n) *n = -bas;
  else *n = bas;
  
  orbfile.close();


}
void GUS_reader::get_occ (double *o, int *n)
{



  ifstream orbfile("orbital.txt");
  
  if (!orbfile.is_open())
    throw string("GUS_reader::get_info : \"orbital.txt\" does not exist.");
  
  ostringstream oss;  
  string input;

  while (getline(orbfile, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");
  
  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 
  istringstream iss;
  

  int alpha, beta;
  if (regex_search(start, end, what,
		   boost::regex("^ NUMBER OF OCCUPIED ORBITALS \\(ALPHA\\)          =\\s+(\\d+)$"), 
		   boost::format_perl)) {
    
    iss.str(string(what[1].first, what[1].second));
    iss >> alpha;
    iss.clear();
    
    
  } else 
    throw string("GUS_reader::get_occ : Could not read the number of alpha occupied orbitals.");

  start = input.begin(); end = input.end();
  if (regex_search(start, end, what,
		   boost::regex("^ NUMBER OF OCCUPIED ORBITALS \\(BETA \\)          =\\s+(\\d+)$"), 
		   boost::format_perl)) {
    
    iss.str(string(what[1].first, what[1].second));
    iss >> beta;
    iss.clear();
    
    
  } else 
    throw string("GUS_reader::get_occ : Could not read the number of beta occupied orbitals.");
  

  orbfile.close();

  int nSym, nBas, nCenters;
  get_info(&nSym, &nBas, &nCenters);

  int num_occ = (alpha == beta) ? alpha : alpha+beta;
  double el_per_orb = (alpha == beta) ? 2. : 1. ;

  for (int i = 0; i < nBas; i++) {
    if (i < *n)
      o[i] = (num_occ > i) ? el_per_orb : 0.;
  }

  if (nBas > *n) *n = -nBas;
  else *n = nBas;

}


void GUS_reader::get_next_Vpqrs(int *p, int *q, int *r, int *s, double *vint)
{

  static ifstream mofile;

  if (!mofile.is_open()){
    mofile.open("moint.txt");
    if (!mofile.is_open())
      throw string("GUS_reader::get_next_Vpqrs : Could find \"moint.txt\".");
  }

  // TODO: fix the bug here
  // NOTE: maybe a pgi/intel difference, this line is for pgi
  // attempts to read and then raises the EOF flag
  mofile >> *p >> *q >> *r >> *s >> *vint;  


  if (mofile.eof()) {
     mofile.close();
    *s = 0;
    return;
  }

  //intel line
  //mofile >> *p >> *q >> *r >> *s >> *vint;
  //  if (*p  == 46)
  //cout << *p << ' ' << *q << ' ' << *r << ' ' << * s << '=' << *vint << endl;


}

void GUS_reader::get_scfvec (double *C, int *n, int *len)
{

  ifstream mocoef("mocoef.txt");
  
  if (!mocoef.is_open())
    throw string("GUS_reader::get_scfvec : \"mocoef.txt\" does not exist.");

  string line;
  istringstream iss;
  int mo_num, ao_num;


  getline(mocoef, line);
  iss.str(line);
  iss >> ao_num >> mo_num;
  iss.clear();

  if ( (*n != mo_num) || (*len != ao_num) ) {
    *len = -ao_num;
    return;
  }


  while (getline(mocoef, line)) {
    iss.str(line);
    int i, j;
    double coef;
    iss >> i >> j >> coef;
    i--; j--;
    C[ i * ao_num + j] = coef;

    iss.clear();
    
  }

  mocoef.close();

}
void GUS_reader::get_overlap (double *S, int *n)
{
  ifstream mocoef("overlap.txt");
  
  if (!mocoef.is_open())
    throw string("GUS_reader::get_overlap : \"overlap.txt\" does not exist.");

  istringstream is;
  string line;

  unsigned row = 0, col = 0;
  while (getline (mocoef,line) ) {

    if ( (row < *n)  && (col < *n)) {
      is.str(line);
      is >> S[row * (*n) + col];
      S[col * (*n) + row] =  S[row * (*n) + col];
      is.clear();
      // Its a triangular matrix
    }
      
    if (row == col) {
      col = 0;
      row++;
    } else
      col++;
    
  }


  if (!((row != *n)  && (col != *n)))
    *n = col > row ? col : row;
  else
    *n = col > row ? -col : -row;
    
  mocoef.close();
  

}


void GUS_reader::get_dip (double *x, double *y, double *z, int *n)
{


  int nSym, nBas = 0, nCent;
  get_info(&nSym, &nBas, &nCent);
  if (*n != nBas) {*n = -nBas; return;}

  ifstream dip_file("dipAO.txt");
  
  if (!dip_file.is_open())
    throw string("GUS_reader::get_info : \"dipAO.txt\" does not exist.");
  
  ostringstream oss;  
  string input;

  while (getline(dip_file, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");

  dip_file.close();
  
  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 
  istringstream iss;
  
  vector<vector<double> > data;



  int rows;
  while (regex_search(start, end, what,
		      boost::regex("^\\s+(\\d+)\\s+\\w+\\s+\\d+\\s+\\w+\\s+((:?-?\\d+\\.\\d+(:?\\s+)?)+)$"), 
		      boost::format_perl)) 
    {
      start = what[0].second;
     
      iss.str(string(what[1].first, what[1].second));
      iss >> rows;
      iss.clear();

      if (data.size() < rows) 
	data.push_back(vector<double>(0));
      
     
      
      string vals(what[0].first, what[0].second);
      
      double num;
      string::const_iterator list_start=vals.begin(), list_end=vals.end();
      while (regex_search(list_start, list_end, what,
			  boost::regex("-?\\d+\\.\\d+"), 
			  boost::format_perl)) 
	{
	  
	  iss.str(string(what[0].first, what[0].second));
	  iss >> num;
	  iss.clear();

	  data[rows-1].push_back(num);
	  list_start = what[0].second;
	}
      
    }

  int size = data.size();

  
  Blas_matrix xAO(size,size), yAO(size,size), zAO(size,size);
  
  // REad the AO dipole matrices
  for(int i = 0; i < data.size(); i++) {
    int col_size = data[i].size()/3;
    for(int j =0; j < col_size;j++) {

      xAO(i,j) = xAO(j,i) = data[i][j];
      yAO(i,j) = yAO(j,i) = data[i][j+col_size];
      zAO(i,j) = zAO(j,i) = data[i][j+2*col_size];
    }
  }



  //MO Transformation
  int c_len = 0;
  get_scfvec (0, &nBas, &c_len);
  c_len *= -1;
  Blas_matrix C(c_len,nBas);
  // Load the SCF vectors
  get_scfvec (&C(0,0), &nBas, &c_len);

  Blas_matrix xMO(nBas,nBas), yMO(nBas,nBas), zMO(nBas,nBas);
  xMO = 0; yMO = 0; zMO = 0;

  Blas_matrix T(nBas,c_len);

  T = 0;
  T.dgemm('T','N',C,xAO);
  xMO.dgemm('N','N',T,C);

  T = 0;
  T.dgemm('T','N',C,yAO);
  yMO.dgemm('N','N',T,C);


  T = 0;
  T.dgemm('T','N',C,zAO);
  zMO.dgemm('N','N',T,C);


  // Load the MO matrices
  for(int i = 0; i < nBas; i++) 
    for(int j = 0; j < nBas;j++) {
      z[i+j*nBas] = zMO(i,j);
      x[i+j*nBas] = xMO(i,j);
      y[i+j*nBas] = yMO(i,j);
      
    }
  
   
}


extern "C" void compute_cap_mo_(double*, double*, double*, int*, int*);

void GUS_reader::get_mo_cap(double *boxx, double *boxy, double *boxz, int* nmo, double* capmo)
{
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
    throw string("GUS_reader::get_mo_cap : \"capmat.txt\" does not exist.");
 
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
