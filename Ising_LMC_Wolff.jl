# Single-core Local MC, Swendsen-Wang, Wolff, Wang-Landau on D-dim lattice

using DelimitedFiles
using Statistics
using Random
using SpecialFunctions

Random.seed!(1234)

const L = 4
const DIM = 2
const N = L^DIM
const DIM1 = DIM + 1

### Critical Temperature: 2D Ising:  0.440686   3D Ising:  0.22165
const beta = 0.440686
const Pcluster = 1 - exp(-2 * beta)
#####           Initial settings

#const theta = Array{Float64, 1}(undef, N)
const s = Array{Int64, 1}(undef, N)
s[:] .= 1

const NeighborList = Array{Int64, 2}(undef, 2*DIM, N)
for j in 1 : DIM
    Lj = L^(j-1)
    for i = 0 : N -1
        NeighborList[2*j-1, i+1] = i + ( mod( fld(i, Lj) + 1, L)  - mod(fld(i, Lj), L) ) * Lj + 1
        NeighborList[2*j, i+1] = i + ( mod( fld(i, Lj) - 1, L)  - mod(fld(i, Lj), L) ) * Lj + 1
    end
end

const EnergyList = Array{Float64, 1}(undef, 2*DIM)
const ProbaList = Array{Float64, 1}(undef, 2*DIM)

EnergyList[:] = 1:(2*DIM)
ProbaList[:] = exp.(-2 * beta .* EnergyList[:])

#####           Generate a list for recording
struct RecordList
    List::Array{Int64,1}
    Iter::Int64
    Pointer::Int64

    function IfRecord(Iter, Pointer)
        if Iter == List[Pointer]
            Pointer += 1
            return true
        else
            return false
        end
    end
end

#####           Functions of Evaluations

Correlation(n1, n2) = s[n1] * s[n2]
M() = mean(s)

function chi()
    m = sum(s)
    return m*m/N
end


#####           Monte Carlo

### Local Metropolis Method

function MCMove(currentID::Int64)
    #currentID = rand(1:N)
    deltaEInd = sum(s[NeighborList[:, currentID]]) * s[currentID]
    # deltaE = -2 * sum(s[NeighborList[:, currentID]]) * s[currentID]
    # Paccept = exp(-2 * beta * sum(s[NeighborList[:, currentID]]) * s[currentID])
    if (deltaEInd <= 0)
        # accept
        s[currentID] = -s[currentID]
        return 1 
    elseif (rand() < ProbaList[deltaEInd])
        # accept 
        s[currentID] = -s[currentID]
        return 1 
    else
        # reject
        return 0 
    end
end
#=

##### Swendsen-Wang

const edges = Array{Bool, 2}(undef, DIM, N)

# Sites with the same number in the same cluster
const ClusterNumeber = Array{Int64, 1}(undef, N)

function SwendsenWang()
    # Assign the edges
    for ind = 1:N
        edges[:, ind] = (s[ind] .== s[NeighborList[1:2:2DIM, ind]])
    end
    edges .*= (rand(DIM, N) .< Pcluster)
    # Merging Clusters/Percolation

    # return TotalClusterNum
end

=#


##### Wolff

const InPocket = Array{Int64, 1}(undef, N)
# A stack of indices of spins inside the pocket
const Pocket = Array{Int64, 1}(undef, N)

function Wolff(currentID::Int64)
    SizePocket = 1    
    # Default 1, in the pocket: -1
    InPocket .= 1
    InPocket[currentID] = -1
    # Counter: the pointer in enlarging the pocket
    Pocket[SizePocket] = currentID
    Counter = 1

    while (true)
        for ind in NeighborList[:, currentID]
            if (InPocket[ind] != -1)
                if (s[ind] == s[currentID])
                    if (rand()<Pcluster)
                        # insert i in pocket
                        SizePocket += 1
                        Pocket[SizePocket] = ind
                        InPocket[ind] = -1
                    end
                end
            end
        end
        # move to next
        if (SizePocket == Counter)
            break
        end
        Counter += 1; currentID = Pocket[Counter]        
    end
    # Flip 
    s[:] .*= InPocket
end

# test
Cmean = 0
NIter = N*N*N*N

for i = 1:NIter
    global Cmean
    for iter = 1:N
        Wolff(rand(1:N))
    end
    Cmean += Correlation(1,2)
end

Cmean /= NIter
println(Cmean)
# value ~ 0.782




#=

Cmean = 0
NIter = N*N*N*N*N*N
#####                   Test
for i = 1:NIter
    global Cmean
    for iter = 1:N
        MCMove(rand(1:N))
    end
    Cmean += Correlation(1,2)
end
Cmean /= NIter
println(Cmean)
=#



#=


=#





