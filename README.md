# Particle-Based-Fluid-Simulation

## Usage

You could compile code with `sh compile.sh` and modify parameters in `input.txt`, then execute following code:

```
./main input.json data.bin
python3 SPH_painter.py data.bin
```
Start of the simulation:
![Simulation of 5000 particles](sim_pic/simulation_begin.png "simulation")

End of the simulation: 
![Simulation of 5000 particles](sim_pic/simulation_result.png "simulation")

## TODO
- Adding 2D support.
- Adding animation drawer.
- Adding multi thread support.

![Frame](sim_pic/Frame.png "Frame")
