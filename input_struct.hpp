#ifndef __INPUT_STRUCT_HPP__
#define __INPUT_STRUCT_HPP__

#include <set>

#define UNDEF 0

enum Propags {ADC2DIP = 1, ADC2PP, ADC2IPx, NDADC3IP, NDADC3AP};
enum Caps {FULL = 1, SUBSPACE};
enum StepTypes {EXP = 1, LINEAR, POWER};
enum DiagTypes {DiagFull = 1, DiagLanc};

struct Prop_info {
  std::set<unsigned> symms, holes, spin;
  int  method;
};

struct Diag_info {
  int type; 
  int iter;
  std::set<unsigned> vecs;
};


struct Eigen_info {
  double ps; 
  double thresh;
};


struct Cap_info {
  static const int filesize = 64;
  char output[filesize];
  double deltaeta, srcap, sicap;
  double eadcmin, eadcmax;
  int nicap;
  int ioffset;
  int increase;
  int type;
  double boxx, boxy, boxz;
};

// Test only
struct Dip_info {
  double nucl;
};


#endif // #ifndef __INPUT_STRUCT_HPP__
