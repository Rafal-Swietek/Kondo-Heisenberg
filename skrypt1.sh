#PBS -S /bin/bash
#PBS -q main
#PBS -l walltime=95:59:59
#PBS -l select=1:mem=32GB

module load Armadillo

cd $PBS_O_WORKDIR/

g++ main.cpp Hamiltonian.cpp -o program.o -larmadillo -std=c++17

./program.o >& log.txt