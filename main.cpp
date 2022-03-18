#include "SPH_handle.h"

int main()
{   
    SPH_handle simulator(10, 100, 100, 100);
    simulator.set_export_file("data.bin");
    simulator.run(10);

    return 0;
}

