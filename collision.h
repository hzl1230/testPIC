#ifndef _COLL_
#define _COLL_
#define DEBUG_VEL
#include "cross_section.h"
#include "particle.h"
#include "utility.h"
#include "espic_math.h"

typedef std::vector<std::pair<Particle, Particle>> CollPair;
using namespace ESPIC;

class Collisionpair {
public:    // In class all velocity except for the Update part are relative velocity 
    Collisionpair(Particle& particle, Real mass1, Real mass2, Real vth1, Real vth2)
    : pt(particle), mass(mass1),
    vxr(particle.vxr()), vyr(particle.vyr()), vzr(particle.vzr()),
    vth(vth1), vtb(vth2), // inject(1) and background(2) thermo velocity
    F1(mass1/(mass1+mass2)), F2(mass2/(mass1+mass2))
    {
        if (vxr == 0) { theta = 0.5 * ESPIC::PI; }
        else { theta = atan2(sqrt(vyr*vyr+vzr*vzr), vxr); }
        if (vyr == 0) { 
            if (vzr > 0) phi = 0.5 * ESPIC::PI;
            else phi = -0.5 * ESPIC::PI;
        } else phi = atan2(vzr, vyr);
        st = sin(theta);
        ct = cos(theta);
        sp = sin(phi);
        cp = cos(phi);

        Real vxb, vyb, vzb;
        vxb = pt.vx() - vxr;
        vyb = pt.vy() - vyr;
        vzb = pt.vz() - vzr;
        
        wx = F1 * (pt.vx() + vxb);
        wy = F1 * (pt.vy() + vyb);
        wz = F1 * (pt.vz() + vzb);
        
        energy = pt.rel_velsqr() * mass;
        // vel = sqrt(2.0 * pt.rel_velsqr());
    }

    void ParticleElasticCollision() 
    { 
        // std::ofstream of("out/ela.dat", std::ofstream::app);
        // of << "Ori_Eng: "<< energy << std::endl;
        vel = sqrt(2. * pt.rel_velsqr());
        chi = acos(1.0 - 2.0*RG01());
        eta = ESPIC::PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleExcitatinCollision(Real th) 
    {   
        // std::ofstream of("out/exc.dat", std::ofstream::app);
        // of << "Ori_Eng: "<< energy << std::endl;
        energy = fabs(energy - th);
        vel = sqrt(2.0 * energy/mass);
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
        std::ofstream of("out/ion.dat", std::ofstream::app);
        of << "Ori_Eng: "<< energy << " ";
        energy = fabs(energy - th);

        en_ej = w * tan(RG01() * atan(0.5*energy/w));
        en_sc = fabs(energy - en_ej);
        vel = sqrt(2.0 * en_sc/mass);
        vel_ej = sqrt(2.0 * en_ej/mass);
        chi = acos(sqrt(en_sc / energy));
        chi_ej = acos(sqrt(en_ej / energy));
        eta = ESPIC::PI2 * RG01();
        eta_ej = eta + ESPIC::PI;
    #ifdef DEBUG_VEL
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
        of << "vx: " << vxr << " "
           << "vy: " << vyr << " "
           << "vz: " << vzr << std::endl;
    #endif
        UpdateParticleVelInfo();
        Particle newelectron = Particle(pt.x(),pt.y(),pt.z());
        Particle newion = Particle(pt.x(),pt.y(),pt.z());
        EjectEletronReaction(chi_ej, eta_ej, vel_ej, newelectron);
        EjectIonReaction(newion);
    #ifdef DEBUG_VEL
        of << pt.vx() << " "
           << pt.vy() << " "
           << pt.vz() <<std::endl;
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
        vel = sqrt(2. * pt.rel_velsqr());
        chi = acos(1.0 - 2.0*RG01());
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleBackwardCollision()
    {
        vel = sqrt(2. * pt.rel_velsqr());
        chi = PI;
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    std::pair<Particle, Particle>& ion_products() { return particle_pair; }



private:
    Particle& pt;
    const Real mass;
    Real vxr, vyr, vzr;
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
    Real gx, gy, gz;

    gx = vel_ * (ct * cc - st * sc * ce);
    gy = vel_ * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = vel_ * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    particle.vx() = wx + F2*gx;
    particle.vy() = wy + F2*gy;
    particle.vz() = wz + F2*gz;
}

void UpdateParticleVelInfo()
{
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    Real gx, gy, gz;

    gx = vel * (ct * cc - st * sc * ce);
    gy = vel * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = vel * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    pt.vx() = wx + F2*gx;
    pt.vy() = wy + F2*gy;
    pt.vz() = wz + F2*gz;
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