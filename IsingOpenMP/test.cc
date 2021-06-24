#include "spin.h"

int main()
{
    IsingModel ising(3, 1.0);
    ising.SetSize(4);
    ising.SetTemp(0.22165);
    //ising.NeighborList.print();
    ising.OMP_init();

    for (int iter=0; iter<ising.N*ising.N; iter++)
    {
        for (int j=0; j<ising.N; j++)
        {
            std::cout<<LocalMonteCarlo(ising, rand_int(ising.N));
            std::cout<<" ";
        }
        std::cout<<std::endl;
        std::cout<<ising.Chi()<<" "<<ising.E()<<std::endl;
    }
    return 1;
}