#ifndef _UTIL_
#define _UTIL_
#
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iterator>
#include <random>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

typedef double Real;
typedef int Smallint;
typedef int Bigint;
using std::vector;
using std::string;

const Real pi = 2.*asin(1);

void espic_error(const std::string& msg)
{
  std::cerr << "--> Error: " << msg << ". <--" << std::endl;
  exit(-1);
}

string unknown_cmd_info(const std::string& cmd, 
                             const std::string& file)
{
  std::ostringstream oss;
  oss << "Unknown command \"" << cmd << "\" in [" << file << "]";
  return oss.str();
}

string out_bound_info(const std::string& var,
                      const std::string& file)
{
    std::ostringstream oss;
    oss << "Out the boundary of "<< var << "define in file [" << file << "]";
    return oss.str(); 
}

Real random_G01()
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<Real> RDr(0.,1.);
    return RDr(g);
}

void random_sample_Maxwell(size_t n, vector<Real>& x, vector<Real>& v)
{
    Real sq2 = sqrt(2);
    Real vel;
    for (size_t i = 0; i < n; ++i) { 
        vel = sq2*sin(2*pi*random_G01())*sqrt(-log(random_G01()));
        v.emplace_back(vel);
        x.emplace_back(random_G01());
    }
}

void print(vector<Real>::iterator beg, vector<Real>::iterator end)
{
    std::copy(beg, end, std::ostream_iterator<Real>(std::cout, " "));
    std::cout << std::endl;
}

void print(vector<Real>::const_iterator beg, vector<Real>::const_iterator end)
{
    std::copy(beg, end, std::ostream_iterator<Real>(std::cout, " "));
    std::cout << std::endl;
}

void print(vector<string>::const_iterator beg, vector<string>::const_iterator end)
{
    std::copy(beg, end, std::ostream_iterator<string>(std::cout, " "));
    std::cout << std::endl;
}

void fprint(vector<Real>::const_iterator beg, vector<Real>::const_iterator end, std::ofstream& of)
{
    std::copy(beg, end, std::ostream_iterator<Real>(of, " "));
    of << std::endl;
}

inline Real toReal(const string& str) 
{
    return (Real)atof(str.c_str());
}

inline Real vsum(const vector<Real>& v)
{
    Real sum = 0;
    for(size_t i = 0; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

inline Real Pcoll(Real nu, Real dt) 
{
    return 1 - exp(-nu*dt);
} 

inline Real findmax(vector<Real>::const_iterator beg, vector<Real>::const_iterator end)
{
    Real max = 0;
    for(auto i = beg; i != end; ++i) 
        if (max < *beg) max = *beg;

    return max;
    
}

#endif