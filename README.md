# parallel_full_gyro
This is a particle-in-cell code to compare full kinetic and gyrokinetic  simulation with hybrid parallel scheme

0. By the comparison between the gyrokinetic simulation  and full-orbit simulation with the presence of the perturbation which has large amplitude and short wavelength (under the periodic and natural boundary conditions (only the periodic one is available now.)), the purpose of this code is to verify the correctness of the gyrokinetic theory and simulations when such kind of perturbation exists. 

1. Only the electrostatic potential is considered in the current version.

2. Since the perpendicular drift of charged particles magnetized in the magnetic field is the focus of gyrokinetic theory,
the simulation domain is confined within the surface perpendicular to the magnetic field, and the parallel dimension is 
not included. 

3. The linked list is used to store the sampling particles and each process is responsible for the particles, the spatial position of which locates at that process. 

4. This code includes the functionalities to trace the orbit of test particles driven by the self-consitent electrostatic field. The number of test particles can be as much as you want and they are also stored by the test particles. 

5. The full orbit can be advanced by boris algorithm, 2nd-order and 4th-order runge-kutta algorihtm. The gyrocenter orbit can be driven by the 2nd-order and 4th-order runge-kutta algorihtm.

6. The double-gyroaveage operation in the gyrokinetic quasi-neutrality equation is computed by the double-gyroaveage interpolation algorithm, given by the paper : https://www.mdpi.com/2571-6182/2/2/9 

7. The sampling of particles along the magnetic moment "mu" dimension can be random sampling along that dimansion, or that by dividing the "mu" domain into a mesh, the number of particles at each node is obtain according to the distribution, and then all the particles with the same "mu" are randomly sampled.

8. The cubic-b spline is used for the finite-element-method.

9. This code can do integrated simulations as well as tracing the single particle orbit.
