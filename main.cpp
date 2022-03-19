#include "SPH_handle.h"

int main()
{   
    SPH_handle simulator(1000, 100, 100, 12);
    simulator.set_export_file("data.bin");
    simulator.run(1000);

    return 0;
}

