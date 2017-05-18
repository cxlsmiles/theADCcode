#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

// The file contains the declaration of the Config structure

// The Config structure defines an electronic configuration. The format
// is specific for the ADC2 double ionization propagator: |ii>,
// |ij>, |jiir>, |ijkr,T> , where T is the type of 
// the spin function (see Chem. Phys. 329, 11, Table A.1)

struct Config {
  unsigned short vir;       // A virtual orbital index for 3h1p configurations
  unsigned char occ[3];     // Occupied orbital indexes for 2h and 3h1p configurations
  unsigned char type;       // Spin function index for <ijkr| configurations
  Config(unsigned char i, unsigned char j, unsigned char k = 0, 
	 unsigned short r = 0, unsigned char t = 0)
  {occ[0] = i; occ[1] = j; occ[2] = k; vir = r; type = t;}    
};  

#endif //#ifndef __CONFIG_HPP__
