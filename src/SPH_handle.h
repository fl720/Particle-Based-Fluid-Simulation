#include <stdio.h>
#include <vector>
#include <map> 
#include <set>
#include <random>
#include <fstream>
#include <string>
#include <chrono>

#include "particle.h"

class SPH_handle
{
    public:

        SPH_handle(std::string filename) ; 
        //SPH_handle(int number, double he = 100, double wi = 100, double le = 100);
        ~SPH_handle();

        void set_export_file(const char * path);
        void run(unsigned int step );

    private:

        FILE * fp;
        unsigned int total_step; 

        double dt ; 
        vector3d volume;

        parameter para;

        unsigned int particle_number;

        std::vector<Particle> particles;
        std::vector<Particle> tem_par  ;
        std::map<cubic_zone , std::set<unsigned int> > particle_list; // spatial hash 

};

class stopwatch
{
    public: 
    
        void start() ; 
        float elapse() ; 

    private:

        std::chrono::time_point<std::chrono::system_clock> start_time ;
};

    