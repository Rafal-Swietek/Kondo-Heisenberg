g++ -g main.cpp Hamiltonian.cpp -o program.o -fopenmp -larmadillo -std=c++17
./program.o >& log.txt
