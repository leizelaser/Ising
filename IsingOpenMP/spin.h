#include "ndarray.h"

#define DIM 3
#define DIM2 DIM*2


#ifdef OMP_H
#include "grid_field_virtual.h"
#endif

class IsingModel
{
    public:
    int L;
    int N;                       // N = L^DIM
    double beta;
    double Pcluster;                // = 1 - exp(-2.0*beta)
    int * s;    
    ndarray<int> NeighborList;
    
    IsingModel();
    IsingModel(int l, double bet);
    ~IsingModel();

    void SetSize(int l);
    void SetTemp(double b);

    // deprecated single-core CPU operation; can be done with multiple-core GPU
    double M();
    double Chi();
    double E();

    // OMP operations: multithread calculation, and multithread Monte Carlo
    #ifdef OMP_H
    double M_omp();
    double Chi_omp();
    double E_omp();

    ndarray<int> Checkerboard;      // two sets of sites
    int halfN;

    int OMP_init();
    #endif
};

IsingModel::IsingModel(): L(0), N(0), s(NULL)
{};

IsingModel::IsingModel(int l, double bet): L(l), beta(bet)
{
    N = pow(L, DIM);
    Pcluster = 1 - exp(-2.0*beta);
    s = new int[N];
    for (int i=0; i<N; i++)
        s[i] = 1;
    NeighborList.resize(N, 2*DIM);
    NeighborList.autolink();

    // NeighborList to be checked
    for (int j=0; j<DIM; j++)
    {
        int Lj = pow(L, j);
        for (int i=0; i<N; i++)
        {
            NeighborList[i][j] = i + (int(floor(i/Lj)+1)%L - int(floor(i/Lj))%L) * Lj;
            NeighborList[i][j+DIM] = i + (int(floor(i/Lj)-1+L)%L - int(floor(i/Lj))%L) * Lj;
        }
    }
}

void IsingModel::SetSize(int l)
{
    L = l; N = pow(L, DIM);
    if (s != NULL) delete[] s;
    s = new int[N];
    for (int i=0; i<N; i++)
        s[i] = 1;
    NeighborList.resize(N, 2*DIM);
    NeighborList.autolink();
    for (int j=0; j<DIM; j++)
    {
        int Lj = pow(L, j);
        for (int i=0; i<N; i++)
        {
            NeighborList[i][j] = i + (int(floor(i/Lj)+1)%L - int(floor(i/Lj))%L) * Lj;
            NeighborList[i][j+DIM] = i + (int(floor(i/Lj)-1+L)%L - int(floor(i/Lj))%L) * Lj;
        }
    }
}

void IsingModel::SetTemp(double b)
{
    beta = b;
    Pcluster = 1 - exp(-2.0*beta);
}


IsingModel::~IsingModel()
{
    if (s != NULL) delete[] s;
}

double IsingModel::M()
{
    double result(0);
    for (int i=0; i<N; i++)
        result += s[i];
    return result;
}

double IsingModel::Chi()
{
    double m = M();
    return beta*m*m/double(N);
}

double IsingModel::E()
{
    double energy = 0;
    for (int i=0; i<N; i++)
        for (int d=0; d<DIM; d++)
            energy -= s[i] * s[NeighborList[i][d]];
    return energy;
}

#include "random_func.h"
// Monte Carlo as external functions:
// ind = rand()%N
int LocalMonteCarlo(IsingModel& Ising, int ind)
{
    double deltaE = 0;
    for (int i=0; i<DIM2; i++)
        deltaE += Ising.s[Ising.NeighborList[ind][i]];
    deltaE *= Ising.beta * Ising.s[ind];
    //Flip or not
    if (deltaE <= 0)
    {
        Ising.s[ind] = -Ising.s[ind];
        return 1;
    }
    else if(rand_unif()<exp(-2 * Ising.beta * deltaE))
    {
        Ising.s[ind] = -Ising.s[ind];
        return 1;
    }    
    else
        return 0;
}

// Parallelization
#ifdef OMP_H

double IsingModel::M_omp()
{
    double result(0);
    #pragma omp parallel for reduction (+:result)
    for (int i=0; i<N; i++)
        result += s[i];
    return result;
}

double IsingModel::Chi_omp()
{
    double m = M_omp();
    return beta*m*m/double(N);
}

double IsingModel::E_omp()
{
    double energy = 0;
    #pragma omp parallel for reduction (+:energy)
    for (int i=0; i<N; i++)
    {
        double tmp = 0;
        for (int d=0; d<DIM; d++)
            tmp -= s[i] * s[NeighborList[i][d]];
        energy += tmp;
    }
    return energy;
}


int IsingModel::OMP_init()
{
    if (N%2 !=0)
    {
        std::cout << "N must be an even integer."<<std::endl;
        return 0;
    }
    halfN = N/2;
    Checkerboard.resize(2, halfN);
    Checkerboard.autolink();
    grid_field_virtual<DIM> indexing;
    indexing.set_size(L);
    int n[2]; n[0] = 0; n[1] = 0;
    vector<DIM, unsigned> tempvec;
    for(int ind=0; ind<N; ind++)
    {
        tempvec = indexing.IndtoVec(ind);
        unsigned j = tempvec.sum()%2;
        Checkerboard[j][n[j]] = ind; n[j]++;
    }
    return 1;
}

// Monte Carlo with omp
int CheckerboardMonteCarlo(IsingModel& Ising, int CheckerboardInd)
{    
    #pragma omp parallel for
    for (int i=0; i<Ising.halfN; i++)
    {
        //std::cout<<omp_get_num_threads()<<std::endl;
        LocalMonteCarlo(Ising, Ising.Checkerboard[CheckerboardInd][i]);
    }
    return 1;
}

#endif