#include "tile.h"
using std::cout;
using std::endl;
Bigint ESPIC::Random::seed = (Bigint)time(NULL);
ESPIC::Random ranf;

int main(int argc, char **argv)
{
    size_type Np = 10000;
    Real dt = 0.1, time = 0.0;
    Real ntime = 100.;
    if (argc > 1) ntime = (Real)atof(argv[1]);
    Tile* tile = new Tile(Np);

    Real massb = tile->background->mass;
    Real vtb = tile->background->vth;
    std::ofstream of("out/out.dat"), rof("out/particle_num.dat"), onu("out/nu.dat");
    // std::ofstream ooof("out/coll.dat");
    // ooof.close();
    // Reaction* & reaction = tile->reaction_arr[0].second;
    // Particles* & particles = tile->spec_arr[0]->particles;
    while(time < ntime) {
        for(size_type ires = 0; ires < tile->num_coll_species(); ++ires)
        {
            int specid = tile->reaction_arr[ires].first;
            Species* & species = tile->spec_arr[specid];
            Particles& particles = *(species->particles);
            Reaction* & reaction = tile->reaction_arr[ires].second;

            int coll = 0;
            Real numax = 0;
            numax = tile -> UpdateMaxCollFreq(ires);
            onu << numax << " ";
            size_type ncoll = static_cast<size_type>(Pcoll(numax,dt)*particles.size());
            Particles subparticles = Particles();
            particles.get_sub_particles(ncoll,subparticles);
            int ntype(reaction->isize());
            CollPair newpp;
            Real mass(species->mass), vth(sqrt(2.*species->temp/mass));
            rof << particles.size() << " " << subparticles.size() << endl;

            for(size_type ip = 0; ip < subparticles.size(); ++ip)
            {
                Particle& pt = subparticles[ip];
                vector<Real>& nu = pt.nu();
                Real nuj(0), rnd(ranf());
                
                int itype = 0;

                while (itype != ntype)
                {
                    nuj += nu[itype];
                    if (rnd < nuj / numax) {
                        Collisionpair collpair = Collisionpair(pt, mass, massb, vth, vtb);
                        tile->ParticleCollision(itype, reaction, collpair, newpp);
                        coll++; break;
                    }
                    ++itype;
                }  
            }
            if(!newpp.empty()) { 
                // vector<Particle> electron, ion;
                // for(size_type ipair = 0; ipair < newpp.size(); ++ipair) {
                //     electron.emplace_back(newpp[ipair].first);
                //     ion.emplace_back(newpp[ipair].second);
                // }
                for(size_type ipair = 0; ipair < newpp.size(); ++ipair){
                    species->particles->append(newpp[ipair].first);
                    tile->spec_arr[specid+1]->particles->append(newpp[ipair].second);
                }
                newpp.clear();
                newpp.shrink_to_fit();
            }
            species->update_tot_energy();
            of << coll << " " << species->totenergy << " ";
        }
        onu << endl;
        of << endl;
        time += dt;
    }
    delete tile;
    return 0;
}