#include "tile.h"
// #include <iomanip>
using std::cout;
using std::endl;
#define DEBUG
Bigint ESPIC::Random::seed = (Bigint)time(NULL);
ESPIC::Random ranf;

int main(int argc, char **argv)
{
    size_type Np = 10000;
    Real dt = 0.1, time = 0.0;
    Real ntime = 10.;
    if (argc > 1) ntime = (Real)atof(argv[1]);
    Tile* tile = new Tile(Np);

    Real massb = tile->background->mass;
    Real vtb = tile->background->vth;
    

    std::ofstream of("out/out.dat"), cof("out/coll.dat");
    std::ofstream newpof("out/nele.dat");
    std::ofstream ooff("out/ion.dat");
    ooff.close();
    while(time < ntime) {
        of << time << " ";
        Real ten = 0;
        for(size_type ires = 0; ires < tile->num_coll_species(); ++ires)
        {
            int specid = tile->reaction_arr[ires].first;
            Species* & species = tile->spec_arr[specid];
            Particles & particles = *(species->particles);
            Reaction* & reaction = tile->reaction_arr[ires].second;

            int coll = 0, ela = 0, exc = 0, ion = 0; 
            Real numax = 0;
            numax = tile -> UpdateMaxCollFreq(ires);

            size_type ncoll = static_cast<size_type>(Pcoll(numax,dt)*particles.size()); 

            int ntype(reaction->isize());
            CollPair newpp;
            Real mass(species->mass), vth(sqrt(2.*species->temp/mass));
            particles.particles_shuffle();   // copy will not change pt.vel

            for(size_type ip = 0; ip < ncoll; ++ip)
            {
                Particle& pt = particles[ip];
                const vector<Real>& nu = pt.nu();
                Real nuj = 0, rnd = ranf();
                int itype = 0;

                while (itype != ntype)
                {
                    nuj += nu[itype];
                    if (rnd < nuj / numax) {
                        Real ori_en, aft_en;
                        Collisionpair collpair = Collisionpair(pt, mass, massb, vth, vtb);
                        if(ires == 0) {    
                            if (itype == 0) ++ela;
                            else if (itype == 1) ++exc;
                            else ++ion;}
                        ori_en = pt.velsqr() * species->mass;
                        tile->ParticleCollision(itype, reaction, collpair, newpp);
                        aft_en = pt.velsqr() * species->mass + pt.lostenergy();

                        cof << itype << " " << "Ori: " << ori_en << " "
                             << "Aft: " << aft_en << endl;

                        coll++; break;
                    }
                    ++itype;
                }  
            }
            if(!newpp.empty()) { 
                for(size_type ipair = 0; ipair < newpp.size(); ++ipair){
                    species->particles->append(newpp[ipair].first);
                    tile->spec_arr[specid+1]->particles->append(newpp[ipair].second);
                    newpof << newpp[ipair].first.vx() << " " 
                           << newpp[ipair].first.vy() << " "
                           << newpp[ipair].first.vz() << endl;
                }
            }
            species->update_tot_energy();
            of << species->particles->size() << " ";
            if (ires == 0)
            of << ela << " " << exc << " " 
               << ion << " " ;
            of << species->totenergy << " ";
            ten += species->totenergy;
        }
        cof << endl;
        of << endl;
        
        time += dt;
    }
    of.close();
    cof.close();
    delete tile;
    return 0;
}