#PBS -S /bin/bash
#PBS -q main
#PBS -l walltime=95:59:59
#PBS -l select=1:mem=32GB

module load Armadillo
module load mkl

cd $PBS_O_WORKDIR/

g++ -g main.cpp Hamiltonian.cpp -o program.o -larmadillo -std=c++17

#valgrind --leak-check=full --show-reachable=yes --track-origins=yes ./program.o >& log.txt
./program.o >& log.txt