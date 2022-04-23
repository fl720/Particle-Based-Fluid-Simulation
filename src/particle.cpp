#include "particle.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#define _USE_MATH_DEFINES

#ifndef M_PI
    #define M_PI 3.1415926
#endif

cubic_zone Particle::get_grid(double h) 
{
    int a = floor(pos.x / h) ; // floor()
    int b = floor(pos.y / h) ; 
    int c = floor(pos.z / h) ; 
    
    cubic_zone cubic_no(a,b,c);
    return cubic_no; 
}

Particle::Particle(vector3d input_pos)
{
    pos = input_pos;

    v = vector3d(0.1,0,0) ; 
}

void Particle::update( std::set<unsigned int> &surrounding_particles,std::vector<Particle> &particles, std::vector<Particle> &tem_par, double dt, parameter &p , vector3d volume) 
{  
    // accelerate
    vector3d g(0, 0, -p.g);

    vector3d fp = get_pressure(surrounding_particles, tem_par , p);
    vector3d fv = get_viscosity(surrounding_particles, tem_par , p);
    vector3d ft = get_tension(surrounding_particles, tem_par , p) ;
    vector3d fg = (g * rho) ;
    
    vector3d a  = (fp + fv + ft + fg) / rho ;
    
    // move
    pos         = pos + v * dt + a * dt * dt / 2; 
    v           = v + a * dt ; 

    // bouce back

    if ( pos.x > volume.x )
    {
        v.x = - v.x ;
        pos.x = 2 * volume.x - pos.x;
    } 
    if ( pos.x < 0 )
    {
        v.x = - v.x ;
        pos.x = - pos.x;
    } 

    if ( pos.y > volume.y  )
    {
        v.y = - v.y ;
        pos.y = 2 * volume.y - pos.y;
    }
    if ( pos.y < 0 )
    {
        v.y = - v.y ;
        pos.y = - pos.y;
    }

    if ( pos.z > volume.z  ) 
    {
        v.z = - v.z ;
        pos.z = 2 * volume.z - pos.z;
    }
    if ( pos.z < 0 ) 
    {
        v.z = - v.z ;
        pos.z = - pos.z;
    }
}

double Particle::distance_squre(vector3d s) 
{
    return (s - pos) * (s - pos);
}

double Particle::get_rho(std::set<unsigned int> &surrounding_particles, std::vector<Particle> &particles , parameter &p )
{
    double ans = 0 ; 

    for(int pa : surrounding_particles)
    {
        ans += p.m * kernel_poly6(pos - particles[pa].pos , p.h );
    }

    return ans ;
}

vector3d Particle::get_pressure(std::set<unsigned int> &surrounding_particles, std::vector<Particle> &particles, parameter &p )
{
    vector3d ans ;

    for(int pa : surrounding_particles)
    {
        ans += (kernel_spiky_gradient( pos - particles[pa].pos , p.h ) * p.m * (p.k * (particles[pa].rho - p.rho0) + p.k*(rho - p.rho0)) / (2*particles[pa].rho) ); 
    }

    ans = ans * -1 ; 

    return ans ; 
}

vector3d Particle::get_viscosity(std::set<unsigned int> &surrounding_particles, std::vector<Particle> &particles, parameter &p )
{
    vector3d f;

    for(int pa : surrounding_particles)
    {
        f += ((particles[pa].v - v) / particles[pa].rho) * kernel_viscosity_laplacian(pos - particles[pa].pos, p.h);
    }

    f *= p.miu * p.m;
    return f;
}

vector3d Particle::get_tension(std::set<unsigned int> &surrounding_particles, std::vector<Particle> &particles , parameter &p)
{
    vector3d f;

    vector3d cs_grad; // n
    double cs_lap = 0;



    for (auto &pa : surrounding_particles)
    {

        if( (particles[pa].pos - pos) * (particles[pa].pos - pos) > p.h * p.h)
        {
            std::cout << "error" << std::endl ;
        }

        cs_grad += kernal_poly6_gradient( pos - particles[pa].pos , p.h ) * p.m / particles[pa].rho ;
        cs_lap  += kernal_poly6_laplacian( pos - particles[pa].pos , p.h ) * p.m / particles[pa].rho ; 
    }

    if(cs_grad.abs() > 0) f = (cs_grad / cs_grad.abs()) * (- p.sigma) * cs_lap ;


    return f;
}

double Particle::kernel_poly6(vector3d r , double h )
{
    if( h * h < r * r ) return 0 ; 
    return 315*(h*h - r * r)*(h * h - r * r)*(h * h - r * r) / ( 64*M_PI*pow(h,9) )  ;
}

vector3d Particle::kernal_poly6_gradient(vector3d r , double h )
{
    if( h * h < r * r ) return vector3d() ; 
    return r * -945.0 / ( 32.0 * M_PI * pow(h,9) ) * pow((h*h - r * r ), 2) ;
}

double Particle::kernal_poly6_laplacian(vector3d r , double h )
{
    if( h * h < r * r ) return 0 ; 
    return 945.0 / ( 32.0 * M_PI * pow(h, 9) ) * (h * h - r * r) * ( r * r * 7 - 3 * h * h ); 
}

vector3d Particle::kernel_spiky_gradient(vector3d r , double h )
{
    double r_abs = r.abs();
    if ( r_abs <= 1e-7 || h < r_abs ) return vector3d();
    return r * (-45*(h - r_abs)*(h - r_abs) / (r_abs * M_PI * pow(h,6))) ;
}

double Particle::kernel_viscosity_laplacian(vector3d r , double h )
{
    if( h * h < r * r ) return 0 ; 
    return 45 * (h - r.abs()) / (M_PI * pow(h, 6));
}

