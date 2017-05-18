#ifndef __INPUT_DATA_HPP__
#define __INPUT_DATA_HPP__


#include <list>
#include <set>
#include <utility>
#include <string>
#include <vector>

struct Prop_info;
struct Diag_info;
struct Eigen_info;
struct Cap_info;
struct Dip_info; // Test only


struct Input_data {

  std::vector<std::string> analyze_input(std::string input, std::vector<std::string> input_sections);
  void get_number(std::string result, int& number);
  void get_number(std::string result, double& number);
  void get_list(std::string result, std::set<unsigned>& numlist);

  void init_frontend(std::vector<std::string> result);
  void init_diag(std::vector<std::string> result);
  void init_eigen(std::vector<std::string> result);
  void init_popana(std::string result);
  void init_propag(std::vector<std::string> result);
  void init_cap(std::vector<std::string> result);
  void init_dip(std::vector<std::string> result);

  bool is_empty_section(std::vector<std::string> keys);

public:

  struct Prop_info* propag;
  struct Diag_info* diag;
  struct Dip_info* dip;
  struct Eigen_info* eig;
  struct Cap_info* cap;

  std::string backend;
  std::list<std::pair<std::string, std::set<unsigned> > > el_groups;
  Input_data();
  ~Input_data();
};


#endif // #ifndef __INPUT_DATA_HPP__
