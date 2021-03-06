This is the README file for MiniDFT 1.06.

A.) Description

Mini-DFT is a plane-wave denstity functional theory (DFT) mini-app for modeling materials.  Given an set of atomic coordinates and pseudopotentials, mini-DFT computes self-consistent solutions of the Kohn-Sham equations using either the LDA, PBE or PBE0 exchange-correlation functionals. For each iteration of the self-consistent field cycle, the Hamiltonian matrix is constructed and then diagonalized. To build the Hamiltonian matrix, Fast Fourier Transforms are used to tranform orbitals from the plane wave basis ( where the kinetic energy is most readily compted ) to real space ( where the potential is evaluated ) and back. Davidson diagonalization is used to compute the orbital energies and update the orbital coefficients.

The mini-DFT mini-app was excised from the general-purpose Quantum Espresso (QE) code.
Quantum Espresso is licensed per the GNU General Public License (GPL).
A copy of the GPL is provided in the 'License' file in this directory.

The mini-DFT distribution consists for four top-level directories.
src: source code and makefiles for mini-DFT
test: small input files for build validation
benchmark: input files for performance measurement
espresso: mirrors the test and benchmark directories, but with QE-compatible input files

B.) How to Build mini-DFT:

B.0) Dependencies: mini-DFT requires ScaLAPACK, BLAS and FFT libraries. Mini-DFT calls the FFT routines through the FFTW3 interface.

B.1) Define the mini-DFT root directory 

> export MDFT_ROOT=/path/to/mini-DFT

B.2) Move to the src directory.

> cd $MDFT_ROOT/src

B.3) Configure manually...
     a) Select the Makefile.system.compiler that most closely resembles your own.
     b) Edit the chosen makefile to your liking.
        NB: The first line of the makefile can be uncommented to enable OpenMP support.
     c) Create a link to the selected makefile.

> ln -s Makefile.system.compiler Makefile

B.4) Initiate the build

> make

Note: parallel make (e.g., -j 4) will not work.

B.5) Validate the build

You can test the build using the small input (small.pbs and small.in). The total energy should match the total energy to within 10^-5 Ry:

!    total energy              =   -2888.59081781 Ry

C.) How to run mini-DFT:

Mini-DFT accepts three mutally-compatible command line options for PBE0 calculations. 

The -in flag (required) is used to identify the input file read by mini-DFT.

The -nbgrp flag (optional) enables band-group parallelism to improve the parallel scaling exact exchange calculation.  The number of task groups must be a divisor of the number of MPI ranks. The default value for bgrp is 1.  The chosen value must evenly divide the total number of MPI tasks used in the calculation.

The -ndiag flag (optional) sets the number of MPI ranks used for diagonalizaton Ndiag must be a square integer. The default value for Ndiag is the largest square integer less than half the total number of MPI ranks per group. The cannot exceed the tasks per group.

For PBE calculations you may use the following flags in lieue of -nbgrp:

The -ntg flag (optional) enables task-group parallelism to improve the parallel scaling of the FFTs.  The number of task groups must be a divisor of the number of MPI ranks. The default 
value for ntg is 1.

The -npool flag (optional) divides the MPI ranks into "pools".  K-points are distributed among pools, and individual k-points can be evaluated with very little communication between pools. The default value for npool is 1.

To run the code, do someting like:
export OMP_NUM_THREADS=6
mpirun -np 32 ./mini_dft -in small.in -nbgrp 4 > small.out

D.) Required Problems:
Input decks for the large problem size are provided in the $MDFT_ROOT/benchmark directory - large.in and large.pbs example pbs script. The Ti.pbe.nml and O.pbe.nml must be in the run 
directory for a successful run. 

E) Verification:
The benchmark runs are validated based on the total energy after one SCF cycle with exact exchange. The total energy should agree 
with the reference value to within 1e-5 Ry.

Reference output listings for the small and large case are in the ../sample_outputs directory.

F.) Timing. 

The time measured by the miniDFT benchmark excludes initialization and finalization stages. It is labeled "Benchmark_Time" (without quotes) and is on the final line of output.

G.) Reporting

For the electronic submission, include all the source and the makefiles used to build on the target platform.  Include all standard output files.

H.) Authorship

The mini-DFT mini-app was excised from the general-purpose Quantum Espresso http://www.quantum-espresso.org/ (QE) code. QE is an open-source program licensed per the GNU General Public License (GPL). A copy of the GPL is provided in the 'License' file in this directory.
