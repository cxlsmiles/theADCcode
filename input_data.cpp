#include "input_data.hpp"
#include "input_struct.hpp"
#include <boost/regex.hpp>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;

// A temporary input parser,
// replace with elparso later?
// No.

// Input structure:
// &frontend   
// &propagator
// &diagonalizer
// &eigen
// &popana --optional
// &cap    --optional
// &dip    --optional, test only

// Module options:
// &frontend: guk, gus, molcas
// &propagator:
// method: ndadc3ip, ndadc3ap, adc2ipx, adc2dip, adc2pol
// sym  <list> --optional
// spin <list> --optional
// &diagonalizer:
// full, lanczos, ... (Davidson?, some day)
// iter  <num>
// roots <list> --optional 
// &eigen
// ps = <num>   -- pole strength threshold
// thresh = <num> -- component size threshold
// &cap
// full, subspace
// Nicolas' input:
// 5 5 5 ! box size in au x,y,z
// 4 ! nel // not needed anymore
// 40 100 0.015 'exp' ! nIcap, ioffset, sIcap, increase
// 0.0 ! sRcap
// 0.0 85.0 ! eADCmin eADCmax
// 'test.dat'
// &popana
// <name_group> <list of atomic orbitals>
// &dip



Input_data::Input_data() :
  propag(0), diag(0), dip(0), eig(0), cap (0)
{ 



  string input;
  ostringstream oss;  
  // load the input line by line
  while (getline(cin, input))     
    oss << input << endl;    
  input = oss.str();
  oss.str("");

  input = // strip comments 
    boost::regex_replace(input, boost::regex("#.*?$"), "",boost::match_default | boost::format_perl);


  // map of keywords
  // Define the main sections of the input 
  vector<string> keys;
  keys.push_back("^\\s*\\&\\s*frontend\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*propagator\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*diagonalizer\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*eigen\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*cap\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*popana\\s*\\n(.*?)(?=(\\&|\\Z))");
  keys.push_back("^\\s*\\&\\s*dip\\s*\\n(.*?)(?=(\\&|\\Z))"); // test only 


  // get the input strings corresponding to the sections
  vector<string> sections = analyze_input(input, keys);

  keys.clear();
  input.clear();


  //1. Define the frontend
  keys.push_back("(?:^|;)\\s*(\\w+)\\s*(?=(?:;|$))"); // name of the front end
  vector<string> result = analyze_input(sections[0], keys);
  sections[0].clear();
  keys.clear();
  init_frontend(result);
  
  
  //2. Define the keyword patterns of the propagator section
  keys.push_back("(?:^|;)\\s*(\\w+)\\s*(?=(?:;|$))"); //name of the propagator
  keys.push_back("(?:^|;)\\s*spin\\s+((?:\\d\\s*\\.\\.\\s*\\d\\s*|\\d\\s*)+)");
  keys.push_back("(?:^|;)\\s*sym\\s+((?:\\d\\s*\\.\\.\\s*\\d\\s*|\\d\\s*)+)");
  //keys.push_back("^\\s*holes\\s+((?:\\d\\s*\\.\\.\\s*\\d\\s*|\\d\\s*)+)");
 
  result = analyze_input(sections[1], keys);
  sections[1].clear();
  keys.clear();
  init_propag(result);
  
  //3. Define the keywords of the diagonalizer section
  keys.push_back("(?:^|;)\\s*(\\w+)\\s*(?=(?:;|$))"); // name of the diagonalizer
  keys.push_back("(?:^|;)\\s*iter\\s+(\\d+)");
  keys.push_back("(?:^|;)\\s*roots\\s+((?:\\d\\s*\\.\\.\\s*\\d\\s*|\\d\\s*)+)");

  result = analyze_input(sections[2], keys);
  sections[2].clear();
  keys.clear();
  init_diag(result);

  //4. Define the keywords of the eigen section
  keys.push_back("(?:^|;)\\s*ps\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*thresh\\s+(\\d+\\.?\\d*|\\.\\d+)");

  result = analyze_input(sections[3], keys);
  sections[3].clear();
  keys.clear();
  init_eigen(result);

  //5. Define the keywords of the CAP section
  keys.push_back("(?:^|;)\\s*(\\w+)\\s*(?=(?:;|$))"); // name of the CAP
  keys.push_back("(?:^|;)\\s*boxx\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*boxy\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*boxz\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*nicap\\s+(\\d+)");
  keys.push_back("(?:^|;)\\s*ioffset\\s+(\\d+)");
  keys.push_back("(?:^|;)\\s*sicap\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*incr\\s+(\\w+)");
  keys.push_back("(?:^|;)\\s*deltaeta\\s+([+-]?\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*srcap\\s+([+-]?\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*eadcmin\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*eadcmax\\s+(\\d+\\.?\\d*|\\.\\d+)");
  keys.push_back("(?:^|;)\\s*output\\s+([\\w\\.\\d]+)");

  result = analyze_input(sections[4], keys);
  sections[4].clear();
  keys.clear();
  init_cap(result);

  //6. Define the popana section
  keys.push_back("(?:^|;)\\s*((?:(?:\\w+\\s+(?:\\d+\\s*\\.\\.\\s*\\d+\\s+|\\d+\\s+)*(?:\\d+\\s*\\.\\.\\s*\\d+|\\d+))(?:\\s|\\n|;|$)+)+)");
  result = analyze_input(sections[5], keys);
  sections[5].clear();
  keys.clear();
  init_popana(result[0]);
  //7. Define the dipole moment section, TEST ONLY
  keys.push_back("(?:^|;)\\s*nucl\\s+([+-]?\\d+\\.?\\d*|\\.\\d+)"); //any word
  result = analyze_input(sections[6], keys);
  sections[6].clear();
  keys.clear();
  init_dip(result);

}



vector<string> Input_data::analyze_input(string input, vector<string> input_sections)
{

  vector<string> sections_found;  
  
  // Variables of the regexp search.
  string::const_iterator start, end;
  start = input.begin(); end = input.end();
  boost::match_results<std::string::const_iterator> what; 

  for (unsigned int sec = 0; sec < input_sections.size(); sec++) {
    if (regex_search(start, end, what, boost::regex(input_sections[sec]), boost::format_perl)) {
      
      sections_found.push_back(string(what[1].first, what[1].second)); 
      input = boost::regex_replace(input,boost::regex(string(what[0].first, what[0].second)), 
				   "", boost::format_first_only  | boost::format_perl);

      start = input.begin(); end = input.end();
     
    } else 
      
      sections_found.push_back(""); 

  }
  if (start != end) { // Find if there is something  left in the string
    
    string nokey;
    if (regex_search(start, end, what, boost::regex("([^\\s;]+)"), boost::format_perl)) {
      nokey.assign(what[1].first, what[1].second);
      
      throw string("Input error found near:\n") + nokey;
    }
  }


  return sections_found;
  

}


void Input_data::init_frontend(std::vector<std::string> result)
{
  backend = result[0];
}


void Input_data::init_propag(vector<string> prop)
{
  if (is_empty_section(prop)) return;

  propag = new struct Prop_info;

  //Get the type of propagator
  if(regex_match(prop[0], boost::regex("ADC2DIP", boost::regex::icase))) 
    
    propag->method = ADC2DIP;
  
  else if(regex_match(prop[0], boost::regex("ADC2POL", boost::regex::icase)))
    
    propag->method = ADC2PP;
  
  else if(regex_match(prop[0], boost::regex("NDADC3IP", boost::regex::icase)))
    
    propag->method = NDADC3IP;

  else if(regex_match(prop[0], boost::regex("ADC2IPx", boost::regex::icase)))

    propag->method = ADC2IPx;
  
  else if(regex_match(prop[0], boost::regex("NDADC3AP", boost::regex::icase)))

    propag->method = NDADC3AP;
  
  else
    
    propag->method = UNDEF;
  
  //Get spin
  if (!prop[1].empty()) {
    get_list(prop[1],propag->spin);
  } 

  
  
  //Get symmetries
  if (!prop[2].empty()) {
    get_list(prop[2],propag->symms);
  }
  
//   //Get occupied list
//   if (!prop[5].empty()) {
//     get_list(prop[5],propag->holes);
//   } 


}


bool Input_data::is_empty_section(vector<string> keys)
{

  bool is_empty = true;
  
  for(int i = 0; i < keys.size(); i++) {
    if (!keys[i].empty()) { is_empty = false; break;}
  }
  return is_empty;

}



void Input_data::init_diag(vector<string> keys)
{
  
  if (is_empty_section(keys)) return;
  


  diag = new struct Diag_info;
  //Get the type of analysis
  if(regex_match(keys[0], boost::regex("FULL\\w*", boost::regex::icase)))
    
    diag->type = DiagFull;

  else if(regex_match(keys[0], boost::regex("LANC\\w*", boost::regex::icase)))

    diag->type = DiagLanc;

  else

    diag->type = UNDEF;


  if (!keys[1].empty())
    get_number(keys[1],diag->iter);
  else
    diag->iter = UNDEF;

  if (!keys[2].empty()) {
    get_list(keys[2],diag->vecs);
  } 

}


void Input_data::init_eigen(vector<string> keys)
{

  if (is_empty_section(keys)) return;
  
  eig = new struct Eigen_info;
    

  //Get the type of analysis
  if (!keys[0].empty())
    get_number(keys[0],eig->ps);
  else
    eig->ps = 0.;

  if (!keys[1].empty()) {
    get_number(keys[1],eig->thresh);
  } else
    eig->thresh = 0.;
  

}




void Input_data::init_cap(vector<string> keys)
{


  if (is_empty_section(keys)) return;


  cap = new struct Cap_info;

  //Get the type of analysis
  if(regex_match(keys[0], boost::regex("FULL\\w*", boost::regex::icase)))
    
    cap->type = FULL;

  else if(regex_match(keys[0], boost::regex("subspace\\w*", boost::regex::icase)))
    
    cap->type = SUBSPACE;
  
  else
    
    cap->type = UNDEF;

  if (!keys[1].empty())
    get_number(keys[1],cap->boxx);
  else
    cap->boxx = 0.;

  if (!keys[2].empty())
    get_number(keys[2],cap->boxy);
  else
    cap->boxy = 0.;

  if (!keys[3].empty())
    get_number(keys[3],cap->boxz);
  else
    cap->boxz = 0.;

  if (!keys[4].empty())
    get_number(keys[4],cap->nicap);
  else
    cap->nicap = 0;

  if (!keys[5].empty())
    get_number(keys[5],cap->ioffset);
  else
    cap->ioffset = 0;

  if (!keys[6].empty())
    get_number(keys[6],cap->sicap);
  else
    cap->sicap = 0.;


  if(regex_match(keys[7], boost::regex("exp\\w*", boost::regex::icase)))
    
    cap->increase = EXP;
  
  else if(regex_match(keys[7], boost::regex("lin\\w*", boost::regex::icase)))
    
    cap->increase = LINEAR;
  
  else if(regex_match(keys[7], boost::regex("pow\\w*", boost::regex::icase)))

    cap->increase = POWER;

  else
    
    cap->increase = EXP;//default

  if (!keys[8].empty())
    get_number(keys[8],cap->deltaeta);
  else
    cap->deltaeta = 0.;
    ;

  if (!keys[9].empty())
    get_number(keys[9],cap->srcap);
  else
    cap->srcap = 0.;

  if (!keys[10].empty())
    get_number(keys[10],cap->eadcmin);
  else
    cap->eadcmin = 0.;

  if (!keys[11].empty())
    get_number(keys[11],cap->eadcmax);
  else
    cap->eadcmax = 0.;

  
  if (keys[12].empty()) 
    keys[12] = "cap_traj.dat";
 
  // Create a null terminated string
  int slen = keys[12].size() > cap->filesize-1 ? cap->filesize-1 : keys[12].size();
  memcpy(cap->output, keys[12].c_str(),slen);
  cap->output[slen] = '\0';


  
}

void Input_data::init_dip(vector<string> keys)
{


  if (is_empty_section(keys)) return;



  dip = new struct Dip_info;


  if (!keys[0].empty())
    get_number(keys[0],dip->nucl);
  else
    dip->nucl = 0.;


}



Input_data::~Input_data() {
  delete propag;
  delete diag;
  delete dip;
  delete eig;
  delete cap;
}


// Extracts the information for the population analysis in the form of a list of named sets
// using get_list
void Input_data::init_popana(string result) {
  
  boost::match_results<std::string::const_iterator> res_what;
  string::const_iterator res_start = result.begin(),
    res_end  = result.end();

  // Extract one or more lists each having a name tag
  while (regex_search(res_start, res_end, res_what, 
		      boost::regex("^\\s*(\\w+)\\s+((?:\\d+\\s*\\.\\.\\s*\\d+\\s*|\\d+\\s*)+)(?:\\s|\\n|;|$)+"),
		      boost::format_perl))
    
    {
      
      // A pair that contains the name the set and its name.
      pair<string, set<unsigned> > named_set;
      named_set.first = string(res_what[1].first, res_what[1].second);       // First is the name
      get_list(string(res_what[2].first, res_what[2].second), named_set.second); // Second is the list
      el_groups.push_back(named_set);
      
      res_start = res_what[0].second;
    }
}





// Auxiliary methods

// Extracts a number from a regular expression
void Input_data::get_number(string result, int& number) 
{
  istringstream iss(result);
  iss >> number;  
}

void Input_data::get_number(string result, double& number) 
{
  istringstream iss(result);
  iss >> number;
  
}


// Extracts a set of numbers in the form 1 2 4 6..8 4..1, etc.
void Input_data::get_list(string result, set<unsigned>& numlist) {
	
  istringstream iss;
  
  string::const_iterator res_start = result.begin(),
    res_end = result.end();

  boost::match_results<std::string::const_iterator> res_what;
  
  while (regex_search(res_start, res_end, res_what, 
		      boost::regex("(\\d+)\\s*\\.\\.\\s*(\\d+)|(\\d+)"), 
		      boost::format_perl))
    {
      int list_start, list_end;
      
      // If just a number is found
      if (res_what[3].matched) {
	iss.clear();
	iss.str(string(res_what[3].first, res_what[1].second));	      
	
	iss >> list_start;
	list_end = list_start;
	
      }
      
      // If a region is found
      if (res_what[2].matched) {
	iss.clear();
	iss.str(string(res_what[1].first, res_what[1].second));	      
	iss >> list_start;
	
	iss.clear();
	iss.str(string(res_what[2].first, res_what[2].second));	      
	iss >> list_end;
	
	int temp;
	
	if (list_start > list_end) {
	  
	  temp = list_start;
	  list_start = list_end;
	  list_end = temp;
	  
	}
	
      }
      // Insert the numbers in the list
      for (unsigned ii = list_start; ii <= list_end; ii++) 
	numlist.insert(ii-1);
      
      res_start = res_what[0].second;
      
    }
}


