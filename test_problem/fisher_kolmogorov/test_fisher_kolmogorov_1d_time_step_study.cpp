#include "fisher_kolmogorov_problem.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
constexpr int n_cells=200; constexpr double final_time=19.2, diffusion=2.0e-4, alpha=2.0;
std::vector<double> solve_case(const std::string &exe,const std::string &graph,double dt){
 char *args[]={const_cast<char *>(exe.c_str()),const_cast<char *>(graph.c_str())};
 femg::fisher_kolmogorov_problem p(final_time,dt); p.init(2,args);
 p.set_edge_diffusion_coefficients(std::vector<double>(p.number_of_edges(),diffusion));
 p.set_vertex_reaction_coefficients(std::vector<double>(p.number_of_vertices(),alpha));
 p.set_vertex_initial_condition(std::vector<double>(p.number_of_vertices(),0.0));
 p.set_time_scheme(femg::fisher_kolmogorov_problem::TimeScheme::backward_euler);
 p.set_newton_parameters(1.0e-11,30); p.set_output_enabled(false); p.set_verbose(false); p.set_coefficients();
 std::vector<double> initial(n_cells+1,0.0); initial[n_cells/2]=0.1;
 p.set_edge_initial_values(0,initial); p.assemble_matrices(); p.solve();
 auto values=p.edge_values(0);
 for(double value:values) if(!std::isfinite(value)||value < -1e-7||value > 1.0+1e-7) throw std::runtime_error("Concentration outside physical range.");
 return values;
}
double l2(const std::vector<double>&v,const std::vector<double>&r){
 double sum=0,h=2.0/n_cells; for(std::size_t i=0;i<v.size();++i){double w=(i==0||i+1==v.size())?.5:1.; sum+=w*(v[i]-r[i])*(v[i]-r[i]);} return std::sqrt(h*sum);
}
double front(const std::vector<double>&v){
 for(int i=n_cells/2;i<n_cells;++i) if(v[i]>=.5&&v[i+1]<.5){double x=-1.+2.*i/n_cells; return x+2./n_cells*(.5-v[i])/(v[i+1]-v[i]);} return 1.;
}
double mean(const std::vector<double>&v){
 double sum=0; for(std::size_t i=0;i<v.size();++i) sum+=((i==0||i+1==v.size())?.5:1.)*v[i]; return sum/n_cells;
}}
int main(int argc,char*argv[]){try{
 std::string graph=argc>=2?argv[1]:"data/interval_weickenmeier_200.txt";
 std::filesystem::path out=argc>=3?argv[2]:"benchmarks/18_fisher_kolmogorov_1d_sensitivity/results"; std::filesystem::create_directories(out);
 std::vector<double>dts={.025,.05,.1,.2,.3,.4}; std::vector<std::vector<double>> solutions;
 for(double dt:dts) solutions.push_back(solve_case(argv[0],graph,dt));
 std::ofstream profiles(out/"time_step_profiles.csv"),summary(out/"time_step_study.csv"); if(!profiles||!summary) throw std::runtime_error("Unable to write output.");
 profiles<<"dt,x,c\n"<<std::setprecision(16); summary<<"dt,l2_error,max_error,front_position,mean_concentration,min_c,max_c\n"<<std::setprecision(16);
 const auto&reference=solutions.front();
 for(std::size_t k=0;k<dts.size();++k){const auto&v=solutions[k]; double max_error=0;
  for(std::size_t i=0;i<v.size();++i){double x=-1.+2.*i/n_cells; profiles<<dts[k]<<","<<x<<","<<v[i]<<"\n"; max_error=std::max(max_error,std::abs(v[i]-reference[i]));}
  double min=*std::min_element(v.begin(),v.end()),max=*std::max_element(v.begin(),v.end());
  summary<<dts[k]<<","<<l2(v,reference)<<","<<max_error<<","<<front(v)<<","<<mean(v)<<","<<min<<","<<max<<"\n";
  std::cout<<"dt = "<<dts[k]<<", front = "<<front(v)<<", range = ["<<min<<", "<<max<<"]\n";
 }}catch(const std::exception&e){std::cerr<<"Error: "<<e.what()<<"\n";return 1;} return 0;}
