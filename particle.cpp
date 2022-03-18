#include "particle.h"
#include <cmath>
#define _USE_MATH_DEFINES

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

    v = vector3d() ; 
}

void Particle::update(double dt, std::set<Particle> &surrounding_particles,  parameter p) 
{  
    // accelerate
    vector3d g(0, 0, -p.g);
    
    vector3d f  = get_pressure(surrounding_particles , p) + (g * rho) + get_viscosity(surrounding_particles , p)  ;
    
    vector3d a  = f / rho ;
    
    // move
    pos         = pos + v * dt + a * dt * dt / 2; 
    v           = v + a * dt ; 
}

double Particle::distance_squre(vector3d s) 
{
    double dist = (s.x - pos.x) * (s.x - pos.x) + (s.y - pos.y) * (s.x - pos.y) + (s.z - pos.z) * (s.z - pos.z); 
    
    return dist ; 
}

double Particle::get_rho(std::set<Particle> &surrounding_particles, double h )
{
    double ans = 0 ; 

    for(auto &pa : surrounding_particles)
    {
        ans += pa.m * kernel_poly6(pos - pa.pos , h );
    }

    return ans ;
}

vector3d Particle::get_pressure(std::set<Particle> &surrounding_particle, parameter &p)
{
    vector3d ans ;

    for(auto &pa : surrounding_particle)
    {
        ans += (kernel_spiky_gradient( pos - pa.pos , p.h ) * pa.m * (p.k * (pa.rho - p.rho0) + p.k*(rho - p.rho0)) / (2*pa.rho) ); 
    }

    ans = ans * -1 ; 

    return ans ; 
}

vector3d Particle::get_viscosity(std::set<Particle> &surrounding_particle, parameter &p)
{
    vector3d f;

    for(auto &pa : surrounding_particle)
    {
        f += ((pa.v - v) / pa.rho) * kernel_viscosity_laplacian(pos - pa.pos, p.h);
    }

    f *= p.miu * m;
    return f;
}

vector3d Particle::get_tension(std::set<Particle> &surrounding_particle, parameter &p)
{
    vector3d f;

    return f;
}

double Particle::kernel_poly6(vector3d r , double h )
{
    double r_abs = r.abs();
    if( h < r_abs ) return 0 ; 
    return 315*(h*h - r_abs*r_abs)*(h*h - r_abs*r_abs)*(h*h - r_abs*r_abs) / ( 64*M_PI*pow(h,9) )  ;
}

vector3d Particle::kernel_spiky_gradient(vector3d r , double h )
{
    double r_abs = r.abs();
    if ( h < r_abs ) return vector3d();
    return r * (-45*(h - r_abs)*(h - r_abs) / (r_abs * M_PI * pow(h,6))) ;
}

double Particle::kernel_viscosity_laplacian(vector3d r , double h )
{
    if ( h < r.abs()) return 0;
    return 45 * (h - r.abs()) / (M_PI * pow(h, 6));
}

