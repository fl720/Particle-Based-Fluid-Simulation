#include "SPH_handle.h"

int main(int argc, char ** argv)
{  
    std::string infile("input.txt");
    if(argc > 1)
        infile = argv[1];
    SPH_handle simulator(infile);
    simulator.set_export_file("data.bin");
    simulator.run(1000);

    return 0;
}

