# N-body-Simulation
Implementation of N-Body simulation <br />
in n-body-serial you can find serial implementation <br />
in n-body-parallel you can find parallel implementation <br />
Compilation <br />
You can specify deltaTime and number of bodies through command-line parameters, respectively.  <br />
by default dt = 0.01 and numberOfBodies = 100 <br />
Serial compilation "gcc <file_name>" <br />
Serial running "./<application_name> 0.1 100" <br />
Parallel compilation "mpicc <file_name>.c -o <file_name> -lm" <br />
Parallel running "mpiexec -n 10 ./<application_name> 0.1 100"