#include "SPH_handle.h"
#include "json.hpp" 

#include <iostream>

void SPH_handle::set_export_file( const char * path )
{
    if(fp == nullptr) fp = fopen( path , "wb" ) ;
    
    fwrite(&particle_number, sizeof(particle_number), 1, fp); 

    for(Particle &p : particles)
    {
        fwrite(&(p.pos), sizeof(p.pos), 1, fp); 
    }
    
    return ;
}

void SPH_handle::run()
{
    stopwatch timer ; 
    timer.start() ; 
    
    for( int i = 0 ; i < total_step ; i++) 
    {
        if(i % 100 == 0) 
        {
            std::cout << "Step: " << i << "\t" << "Elapse:" << timer.elapse() << "ms" << std::endl;
        }

        std::vector< std::set<unsigned int> >  particle_collection ; 

        for(unsigned int j = 0 ; j < particle_number ; j++)
        {
            std::set<unsigned int> surrounding_particles;

            cubic_zone cubic[7];
            cubic[0] = particles[j].get_grid(para.h) ; 
            for(int k = 1 ; k < 7; k++) cubic[k] = cubic[0];
            cubic[1].x += 1 ; 
            cubic[2].x -= 1 ; 
            cubic[3].y += 1 ; 
            cubic[4].y -= 1 ; 
            cubic[5].z += 1 ; 
            cubic[6].z -= 1 ;    
            for(int k = 0 ; k < 7 ; k++)
            {
                for(unsigned int p : particle_list[cubic[k]] )
                {
                    if( particles[j].distance_squre(particles[p].pos) < para.h * para.h ) 
                    {
                        surrounding_particles.insert( p ) ; 
                    }
                }
            }   
            particle_collection.push_back( surrounding_particles ) ; 
            particles[j].rho = particles[j].get_rho(surrounding_particles, particles, para ); 
        }
        
        tem_par = particles ; 

        for(unsigned int j = 0 ; j < particle_number; j++ )
        {   
            // particle_collection[j].erase(j) ;
            particle_list[particles[j].get_grid(para.h)].erase(j) ;
            particles[j].update(particle_collection[j], particles, tem_par, dt,  para , volume); 
            particle_list[particles[j].get_grid(para.h)].insert(j) ; 
            fwrite( &(particles[j].pos), sizeof(particles[j].pos), 1, fp); 
        }
    }

    return ;
}

SPH_handle::SPH_handle( std::string filename) {

    double he;
    double wi;
    double le;

    std::ifstream input_file(filename);

    if (input_file.good() )
    {
        nlohmann::json j ; 
        input_file >> j ;  

        total_step      = j["step"] ; 
        particle_number = j["number"] ; 
        he              = j["container"]["height"] ; 
        wi              = j["container"]["width"] ;
        le              = j["container"]["length"] ; 
        dt              = j["dt"];
        para.h          = j["kernel_radius"];
        para.rho0       = j["rest_density"];
        para.miu        = j["viscosity"];
        para.k          = j["gas_constant"];
        para.sigma      = j["tension_coefficient"] ;
        para.g          = j["gravity"] ; 
    }

    fp          = nullptr;
    
    volume      = vector3d(le, wi, he) ; 
    srand(time(0));
    for(int it = 0; it < particle_number; it++)
    {
        double he_pos = he/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(he/3)));
        double wi_pos = wi/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(wi/3)));
        double le_pos = le/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(le/3)));

        Particle tmp( vector3d(le_pos, wi_pos, he_pos) );
        particles.push_back(tmp); // tmp -> it
        particle_list[tmp.get_grid(para.h)].insert(it);
    }
}

SPH_handle::~SPH_handle() {
    fclose(fp) ; 
}


void stopwatch::start()
{
    start_time = std::chrono::system_clock::now() ; 
}

float stopwatch::elapse()
{
    std::chrono::time_point<std::chrono::system_clock> current_time = std::chrono::system_clock::now() ;
    
    using namespace std::chrono_literals ;
    float elapse = (current_time - start_time )  / 1ms ; 

    return elapse;
}