#ifndef _COLL_
#define _COLL_
#
#include "cross_section.h"
#include "particle.h"
#include "utility.h"
#include "espic_math.h"

typedef std::vector<std::pair<Particle, Particle>> CollPair;
using namespace ESPIC;

class Collisionpair {
public:    // In class all velocity except for the Update part are relative velocity 
    Collisionpair(Particle& particle, Real mass1, Real mass2, Real vth1, Real vth2)
    : pt(particle), 
    vx(pt.vxr()), vy(pt.vyr()), vz(pt.vzr()),
    vth(vth1), vtb(vth2), // inject(1) and background(2) thermo velocity
    F1(mass1/(mass1+mass2)), F2(mass2/(mass1+mass2))
    {
        if (vx == 0) { theta = 0.5 * ESPIC::PI; }
        else { theta = atan2(sqrt(vy*vy+vz*vz), vx); }
        if (vy == 0) { 
            if (vz > 0) phi = 0.5 * ESPIC::PI;
            else phi = -0.5 * ESPIC::PI;
        } else phi = atan2(vz, vy);
        st = sin(theta);
        ct = cos(theta);
        sp = sin(phi);
        cp = cos(phi);
        wx = F1 * vx;
        wy = F1 * vy;
        wz = F1 * vz;
        
        energy = fabs(get_energy(vx, vy, vz));
        vel = sqrt(2.0 * energy);
    }

    void ParticleElasticCollision() 
    { 
        chi = acos(1.0 - 2.0*RG01());
        eta = ESPIC::PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleExcitatinCollision(Real th) 
    {
        energy = fabs(energy - th);
        vel = sqrt(2.0 * energy);
        chi = acos(1.0 - 2.0 * RG01());
        eta = ESPIC::PI2 * RG01();
        UpdateParticleVelInfo();
        
        pt.lostenergy() += th;
    }

    void ParticleIonizationCollision(Real th)
    {
        Real en_ej, en_sc;
        Real vel_ej, chi_ej, eta_ej;
        Real w = 10.3/kTe0;

        // energy = 0.5 * vel * vel;
        energy = fabs(energy - th);
        en_ej = w * tan(RG01() * atan(0.5*energy/w));
        en_sc = fabs(energy - en_ej);
        vel = sqrt(2.0 * en_sc);
        vel_ej = sqrt(2.0 * en_ej);
        chi = acos(sqrt(en_sc / energy));
        chi_ej = acos(sqrt(en_ej / energy));
        eta = ESPIC::PI2 * RG01();
        eta_ej = eta + ESPIC::PI;
    #ifdef DEBUG
        std::ofstream of("out/coll.dat", std::ofstream::app);
        of << "Eng sc: " << en_sc << " "
           << "Eng ej: " << en_ej << std::endl;
        of << "chi sc: " << chi << " "
           << "chi ej: " << chi_ej << std::endl;
        of << "eta sc: " << eta << " "
           << "eta ej: " << eta_ej << std::endl;
        of << "theta: " << theta << " "
           << "phi: " << phi << " "
           << "F1: " << F1 << " "
           << "F2: " << F2 << std::endl;
        of << "vx: " << vx << " "
           << "vy: " << vy << " "
           << "vz: " << vz << std::endl;
    #endif
        UpdateParticleVelInfo();
        Particle newelectron = pt;
        Particle newion = pt;
        EjectEletronReaction(chi_ej, eta_ej, vel_ej, newelectron);
        EjectIonReaction(newion);
    #ifdef DEBUG
        of << newelectron.vx() << " " 
           << newelectron.vy() << " "  
           << newelectron.vz() <<std::endl; 
        of << newion.vx() << " " 
           << newion.vy() << " "  
           << newion.vz() <<std::endl; 
        of << std::endl;
        of.close();
    #endif
        particle_pair = std::make_pair(std::move(newelectron), std::move(newion));
        pt.lostenergy() += th;

        
    }

    void ParticleIsotropicCollision()
    {
        chi = acos(1.0 - 2.0*RG01());
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleBackwardCollision()
    {
        chi = PI;
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    std::pair<Particle, Particle>& ion_products() { return particle_pair; }



private:
    Particle& pt;
    Real vx, vy, vz;
    const Real vth, vtb;
    const Real F1, F2;
    Real chi, eta, theta, phi;
    Real wx, wy, wz;
    Real st, ct, sp, cp;
    std::pair<Particle, Particle> particle_pair;
    Real vel, energy;
    // std::vector<Real> velbuffer;


void EjectEletronReaction(Real chi_, Real eta_, Real vel_, Particle& particle)
{
    Real sc(sin(chi_)), cc(cos(chi_));
    Real se(sin(eta_)), ce(cos(eta_));
    // Real vxb, vyb, vzb;

    vx = vel_ * (ct * cc - st * sc * ce);
    vy = vel_ * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    vz = vel_ * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    particle.vx() = wx + F2*vx;
    particle.vy() = wy + F2*vy;
    particle.vz() = wz + F2*vz;
    // VelBoltzDistr(vth, vxb, vyb, vzb);
    // RelativeVelocity(particle, vxb, vyb, vzb);
}

void UpdateParticleVelInfo()
{
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    vx = vel * (ct * cc - st * sc * ce);
    vy = vel * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    vz = vel * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    pt.vx() = wx + F2*vx;
    pt.vy() = wy + F2*vy;
    pt.vz() = wz + F2*vz;
}

void EjectIonReaction(Particle& particle)
{
    Real vx, vy, vz;
    VelBoltzDistr(vtb, vx, vy, vz);
    particle.vx() = vx; 
    particle.vy() = vy; 
    particle.vz() = vz;
    // VelBoltzDistr(vtb, vxb, vyb, vzb);
    // RelativeVelocity(particle, vxb, vyb, vzb);
}

};


#endif