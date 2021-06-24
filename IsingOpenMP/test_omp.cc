#include <omp.h>
#ifndef OMP_H
#define OMP_H
#endif
#include "spin.h"

int main()
{
    IsingModel ising(3, 1.0);
    ising.SetSize(4);
    ising.SetTemp(0.22165);
    ising.OMP_init();
    
    //ising.NeighborList.print();
    omp_set_num_threads(8);
    for (int iter=0; iter<ising.N; iter++)
    {
        CheckerboardMonteCarlo(ising, 0);
        CheckerboardMonteCarlo(ising, 1);
        std::cout<<ising.Chi()<<" "<<ising.E()<<std::endl;
        std::cout<<ising.Chi_omp()<<" "<<ising.E_omp()<<std::endl;
    }
    exit(1);
    return 1;
}