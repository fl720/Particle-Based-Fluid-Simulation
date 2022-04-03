#include <vector> 
#include <set>
#include "vector3d.h"

struct parameter
{
    constexpr static double m = 1; // mass
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
    double get_rho(std::set<int> &surrounding_particles, std::vector<Particle> &particles, parameter &p);
    void update(std::set<int> &surrounding_particles, std::vector<Particle> &particles, double dt , parameter &p, vector3d volume);

    vector3d pos;
    vector3d v; // velocity
    double rho; 

    bool operator < (const Particle &p) const
    {
        return pos.x < p.pos.x;
    }

private:

    vector3d get_pressure(std::set<int> &surrounding_particles, std::vector<Particle> &particles, parameter &p);
    vector3d get_viscosity(std::set<int> &surrounding_particles, std::vector<Particle> &particles, parameter &p);
    vector3d get_tension(std::set<int> &surrounding_particles, std::vector<Particle> &particles, parameter &p);

    double kernel_poly6(vector3d r ,double h);
    vector3d kernal_poly6_gradient(vector3d r , double h );
    double kernal_poly6_laplacian(vector3d r , double h );
    vector3d kernel_spiky_gradient(vector3d r, double h);
    double kernel_viscosity_laplacian(vector3d r, double h);

};
