#include <vector> 
#include <set>
#include "vector3d.h"

struct parameter
{
    double h ; // kernel radius
    double rho0 ; // rest density
    double miu ; // viscosity
    double k ; // ideal gas constant
    double sigma ; // tension coefficient
    double g ;  //gravitational constant
};

struct cubic_zone 
{
    cubic_zone():x(0),y(0),z(0){}
    cubic_zone(int x_, int y_, int z_):x(x_),y(y_),z(z_){}

    int x ; 
    int y ; 
    int z ; 

    bool operator < (const cubic_zone &p) const
    {
        return x < p.x;
    }
}; 

class Particle 
{
public:

    Particle(vector3d input_pos);
    ~Particle(){};

    cubic_zone get_grid(double h);
    double distance_squre(vector3d s);
    double get_rho(std::set<Particle> &surrounding_particles, double h );
    void update(double dt , std::set<Particle> &surrounding_particles, parameter p);

    vector3d pos;
    vector3d v; // velocity
    double rho; 
    constexpr static double m = 1; // mass

    bool operator < (const Particle &p) const
    {
        return pos.x < p.pos.x;
    }

private:

    vector3d get_pressure(std::set<Particle> &surrounding_particle, parameter &p);
    vector3d get_viscosity(std::set<Particle> &surrounding_particle, parameter &p);
    vector3d get_tension(std::set<Particle> &surrounding_particle, parameter &p);

    double kernel_poly6(vector3d r ,double h);
    vector3d kernel_spiky_gradient(vector3d r, double h);
    double kernel_viscosity_laplacian(vector3d r, double h);

};
