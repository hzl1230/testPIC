#ifndef _PARTICLES_H
#define _PARTICLES_H

#include <vector>
#include "utility.h"

class Particle {
  public:
    // constructors
    // default constructor
    Particle() :
      pos_ {0, 0, 0},
      vel_ {0, 0, 0}, 
      nu_(0)
      { }

    Particle(Real x, Real vx, Real y, Real vy, vector<Real> nu, Real z=0., Real vz=0.) :
      pos_ {x, y, z},
      vel_ {vx, vy, vz},
      nu_(nu)
      { }

    // copy constructor
    Particle(const Particle& other) :
      pos_{other.pos_[0], other.pos_[1], other.pos_[2]},
      vel_{other.vel_[0], other.vel_[1], other.vel_[2]},
      nu_(other.nu_)
      { }
// assignment operator
    Particle& operator=(const Particle& rhs) {
      pos_[0] = rhs.pos_[0]; pos_[1] = rhs.pos_[1]; pos_[2] = rhs.pos_[2];
      vel_[0] = rhs.vel_[0]; vel_[1] = rhs.vel_[1]; vel_[2] = rhs.vel_[2];
      nu_ = rhs.nu_;
      return *this;
    }

    Real& x()  { return pos_[0]; }
    Real& y()  { return pos_[1]; }
    Real& z()  { return pos_[2]; }
    Real& vx() { return vel_[0]; }
    Real& vy() { return vel_[1]; }
    Real& vz() { return vel_[2]; }
    const Real& x()  const { return pos_[0]; }
    const Real& y()  const { return pos_[1]; }
    const Real& z()  const { return pos_[2]; }
    const Real& vx() const { return vel_[0]; }
    const Real& vy() const { return vel_[1]; }
    const Real& vz() const { return vel_[2]; }

    const Real vel_tot() { return sqrt(vel_[0]*vel_[0] 
                      + vel_[1]*vel_[1] + vel_[2]*vel_[2]); }
    vector<Real>& nu() { return nu_; }
    const vector<Real>& nu() const { return nu_; }
    void update_nu(const vector<Real>& nu) { nu_ = nu; }

    Real* pos() { return pos_; }
    Real* vel() { return vel_; }

  private:
    Real pos_[3];
    Real vel_[3];
    vector<Real> nu_;
};

class Particles {
  friend class Particle;
  typedef std::vector<Particle>::const_iterator const_iterator;
  typedef std::vector<Particle>::iterator iterator;

  public:

    typedef std::size_t size_type;
    Particles():
        nparticles(0), 
        cluster(0), 
        scalar(0) { }

    Particles(const std::vector<Real>& x,
              const std::vector<Real>& y,
              const std::vector<Real>& z,
              const std::vector<Real>& vx,
              const std::vector<Real>& vy,
              const std::vector<Real>& vz):
    nparticles(x.size())
    {
        if (nparticles != y.size())
            espic_error("Failed to construct particles because list length is not matched");
        if (nparticles != z.size())
            espic_error("Failed to construct particles because list length is not matched");
        if (nparticles != vx.size())
            espic_error("Failed to construct particles because list length is not matched");
        if (nparticles != vy.size())
            espic_error("Failed to construct particles because list length is not matched");
        if (nparticles != vz.size())
            espic_error("Failed to construct particles because list length is not matched");

        cluster.resize(nparticles);
        for (size_type ip = 0; ip < nparticles; ip++) {
            cluster[ip].x() = x[ip];
            cluster[ip].y() = y[ip];
            cluster[ip].z() = z[ip];
            cluster[ip].vx() = vx[ip];
            cluster[ip].vy() = vy[ip];
            cluster[ip].vz() = vz[ip];
        }
    }

    // copy constructor
    Particles(const Particles& other)
    : nparticles(other.nparticles),
    cluster(other.cluster) 
    { }

    // desctructor
    ~Particles() { 
      nparticles = 0;
    }

    // public methods
    size_type size() const { return nparticles; }

    const std::vector<Real>& x() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].x();
      return scalar;
    }
    const std::vector<Real>& y() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].y();
      return scalar;
    }
    const std::vector<Real>& z() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].z();
      return scalar;
    }
    const std::vector<Real>& vx() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].vx();
      return scalar;
    }
    const std::vector<Real>& vy() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].vy();
      return scalar;
    }
    const std::vector<Real>& vz() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = cluster[ip].vz();
      return scalar;
    }
    const std::vector<Real>& get_particles_energy() 
    {
      resize_scalar();
      for(size_type ip = 0; ip < size(); ip++) {
        Real vx2, vy2, vz2;
        vx2 = cluster[ip].vx() * cluster[ip].vx();
        vy2 = cluster[ip].vy() * cluster[ip].vy();
        vz2 = cluster[ip].vz() * cluster[ip].vz();
        scalar[ip] = 0.5 * (vx2 + vy2 + vz2);
      }
      return scalar;
    }

    // reserve space for storing n particles
    void reserve(size_type n) {
      cluster.reserve(n);
    }

    // append one particle to the end
    void append(const Particle& p) 
    {
        cluster.emplace_back(p);
        ++nparticles;
    }
    // append a list of particles
    void append(const std::vector<Particle>& p_arr)
    {
        for (size_type i = 0; i < p_arr.size(); i++) 
            append(p_arr[i]);
        nparticles += p_arr.size();
    }
    // append particles given bewteen two iterators
    void append(const_iterator beg, const_iterator end)
    {
        cluster.insert(cluster.end(), beg, end);
        nparticles += static_cast<size_type>(end - beg);
    }
    void append(iterator beg, iterator end) {
        cluster.insert(cluster.end(), beg, end);
        nparticles += static_cast<size_type>(end - beg);
    }
    // append a Particles instance
    void append(const Particles& others) 
    {
        cluster.insert(cluster.end(), others.cluster.begin(), others.cluster.end());
        nparticles += others.size();
    }
    void append(Particles* particles) 
    {
        append(*particles);
    }

    void overwrite(size_type id, const Particle& particle) 
    {
      cluster[id] = particle;
    }

    // particles[id1] = particles[id2]
    void overwrite(size_type id1, size_type id2) {
      cluster[id1] = cluster[id2];
    }

    // "particle = particles[id]"
    // fetch particles[id]'s property and store in particle
    void fetch(size_type id, Particle& particle) {
      particle = cluster[id];
    }

    Particle& operator[] (size_type i) { return cluster[i]; }
    const Particle& operator[] (size_type i) const { return cluster[i]; }
    Particle& at(size_type i) { return cluster[i]; }
    const Particle& at(size_type i) const { return cluster[i]; }

    // erase particle with id 
    void erase(size_type id) 
    {
        cluster[id] = cluster.back();
        cluster.pop_back();
        --nparticles;
    }
    // erase n particles starting with id 
    void erase(size_type id, size_type n)
    {
        n = id + n <= nparticles ? n : nparticles - id;
        for(size_type i = std::max(id+n, nparticles-n); id < nparticles; id++) 
            cluster[id] = cluster[i];
        pop_back(n);
    }

    // pop_back n particles from the end of array
    void pop_back(size_type n) 
    {
        cluster.erase(cluster.end()-n, cluster.end());
        nparticles -= n;
    }

    void particles_shuffle() 
    { 
        srand((unsigned int)time(NULL));
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(cluster.begin(), cluster.end(), g);
    }

    void get_sub_particles(size_type n, Particles*& sub)
    {
        particles_shuffle();
        if (sub->size() != 0) {
            sub->cluster.clear();
            sub->nparticles = 0;
        }   
        std::copy(cluster.begin(), cluster.begin() + n, back_inserter(sub->cluster));
        sub->nparticles += n;
    }

  private:
    size_type nparticles;

    std::vector<Particle> cluster;
    std::vector<Real> scalar;

    void resize_scalar() {
      if(scalar.size() != size()) scalar.resize(size());
    }
};
#endif