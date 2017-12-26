# N-body-Simulation
Implementation of N-Body simulation <br />
in n-body-serial you can find serial implementation <br />
in n-body-parallel you can find parallel implementation <br />
### Compilation <br />
You can specify number of bodies, simulation tim and deltaTime through command-line parameters, respectively. <br />
by default numberOfBodies = 100, simulationTime = 1.0 and dt = 0.1 <br />
Serial compilation "gcc <file_name>" <br />
Serial running "./<application_name> 0.1 100" <br />
Parallel compilation "mpicc <file_name>.c -o <file_name> -lm" <br />
Parallel running "mpiexec -n 10 ./<application_name> 100 1.0 0.1"
10 - number of processes
### Notes 
number of bodies should be divisible evenly by number of processes, <br />
in order for each process to work with same amount of tasks
