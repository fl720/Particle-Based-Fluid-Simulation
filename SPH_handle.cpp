#include "SPH_handle.h"
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

void SPH_handle::run(unsigned int step )
{
    total_step += step ; 
    for( int i = 0 ; i < step ; i++) 
    {
        if(i % 100 == 0) std::cout << i << std::endl;

        // TODO: save surrounding_particles
        for(int j = 0 ; j < particle_number ; j++)
        {
            std::set<Particle> surrounding_particles;
            
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
                        surrounding_particles.insert(particles[p]) ;   
                }
            }     
                        
            particles[j].rho = particles[j].get_rho(surrounding_particles, para.h ); 
        }

        for(int j = 0 ; j < particle_number; j++ )
        {
            std::set<Particle> surrounding_particles;
            
            cubic_zone cubic[7];
            cubic[0] = particles[j].get_grid(para.h) ; 
            for(int k = 1 ; k < 7; k++) cubic[k] = cubic[0];
            cubic[1].x += 1 ; 
            cubic[2].x -= 1 ; 
            cubic[3].y += 1 ; 
            cubic[4].y -= 1 ; 
            cubic[5].z += 1 ; 
            cubic[6].z -= 1 ;  
            // 'rho's in particles[] has been updated nor in particle_list
            for(int k = 0 ; k < 7 ; k++)
            {
                for(unsigned int p : particle_list[cubic[k]] )
                {
                    if( particles[j].distance_squre(particles[p].pos) < para.h * para.h &&
                        j != p
                    ) 
                        surrounding_particles.insert(particles[p]) ;   
                }
            }     

            particle_list[particles[j].get_grid(para.h)].erase(j) ;
            particles[j].update(dt, surrounding_particles, para , volume); 
            particle_list[particles[j].get_grid(para.h)].insert(j) ; 
            fwrite( &(particles[j].pos), sizeof(particles[j].pos), 1, fp); 
        }
    }

    return ;
}

SPH_handle::SPH_handle( int number , double le, double wi , double he) {
    total_step  = 1 ; 
    dt          = 0.001 ; 

    fp          = nullptr;
        
    para.h           = 0.16; // kernel radius
    para.rho0        = 998 ; // rest density    
    para.miu         = 100; // viscosity
    para.k           = 13.8; // ideal gas constant
    para.sigma       = 72.75; // tension coefficient
    para.g           = 9.81 ; // gravitational constant 
    
    volume      = vector3d(le, wi, he) ; 

    particle_number = number ;
    srand(time(0));
    for(int it = 0; it < number; it++)
    {
        double he_pos = he/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(he/3)));
        double wi_pos = wi/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(wi/3)));
        double le_pos = le/3 + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(le/3)));

        Particle tmp( vector3d(le_pos, wi_pos, he_pos) );
        particles.push_back(tmp);
        particle_list[tmp.get_grid(para.h)].insert(it);
    }
}

SPH_handle::~SPH_handle() {
    fwrite(&total_step, sizeof(total_step), 1, fp); 
    fclose(fp) ; 
}
