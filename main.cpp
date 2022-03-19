#include "SPH_handle.h"

int main()
{   
    SPH_handle simulator(10, 12, 9, 9);
    simulator.set_export_file("data.bin");
    simulator.run(1000);

    return 0;
}

