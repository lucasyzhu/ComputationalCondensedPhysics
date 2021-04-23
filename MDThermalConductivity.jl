# 2021/04/23 TODO: fuction initializeNeighbor
using LinearAlgebra
const KB = 8.617e-5
const TIME_UNIT = 1.018e+1
const KAPPA_UNIT = 1.574e+5

function main()
    nx = 4
    ny = 4
    nz = 4
    n0 = 4
    N = n0*nx*ny*nz
    Ne = 2e+4
    Np = 2e+4
    Ns = 10
    # number of heat current list
    Nd = Np/Ns
    Nc = Nd/10
    maxNeighbors = 100
    
    T0 = 60.0
    ax = 5.4
    ay = ax
    az = ax
    lx = ax*nx
    ly = ax*ny
    lz = az*nz
    # cutoff distance for neighbor list
    cutoff = 10.0
    timeStep = 10.0

    NN = zeros(Int32, N, 1)
    NL = zeros(Int32, N*maxNeighbors, 1)
    m = zeros(Float16, N, 1)
    x = zeros(Float32, N, 1)
    y = zeros(Float32, N, 1)
    z = zeros(Float32, N, 1)
    vx = zeros(Float32, N, 1)
    vy = zeros(Float32, N, 1)
    vz = zeros(Float32, N, 1)
    fx = zeros(Float32, N, 1)
    fy = zeros(Float32, N, 1)
    fz = zeros(Float32, N, 1)
    hx = zeros(Float32, Nd, 1)
    hy = zeros(Float32, Nd, 1)
    hz = zeros(Float32, Nd, 1)

    # initialize mass
    for n in 1:N
        m[n] = 40.0
    end
    # initialize position
    x, y, z = initializePosition(n0, nx, ny, nz, ax, ay, az, x, y, z)
    # initialize velocity
    vx, vy, vz = initializeVelocity(N, T0, m, vx, vy, vz)
    # initialize neighbor list
    initializeNeighbor(N, NN, NL, x, y, z, lx, ly, lz, maxNeighbors, cutoff)
end

function initializePosition(n0, nx, ny, nz, ax, ay, az, x, y, z)
    # FCC lattice
    x0 = [0.0, 0.0, 0.5, 0.5]
    y0 = [0.0, 0.5, 0.0, 0.5]
    z0 = [0.0, 0.5, 0.5, 0.0]
    n = 1
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                for i = 1:n0
                    x[n] = ((ix - 1) + x0[i])*ax
                    y[n] = ((iy - 1) + y0[i])*ay
                    z[n] = ((iz - 1) + z0[i])*az
                    n += 1
                end
            end
        end
    end
    return x, y, z
end

function initializeVelocity(N, T0, m, vx, vy, vz)
    averageMomentum = [0.0, 0.0, 0.0]
    temperature = 0.0
    for n = 1:N
        # Generate a random number between (0, 1)
        # It can be also given based on Maxwell's distribution
        vx[n] = rand()*2 - 1 
        vy[n] = rand()*2 - 1
        vz[n] = rand()*2 - 1
        averageMomentum[1] += m[n]*vx[n]/N
        averageMomentum[2] += m[n]*vy[n]/N
        averageMomentum[3] += m[n]*vz[n]/N
    end
    
    for n = 1:N
       # ???
        vx[n] -= averageMomentum[1]/m[n]
        vy[n] -= averageMomentum[2]/m[n]
        vz[n] -= averageMomentum[3]/m[n]
    end
    # scale velocity
    for n = 1:N
        vsquare = vx[n]*vx[n] + vy[n]*vy[n] + vz[n]*vz[n]
        temperature += m[n]*vsquare
    end
    temperature /= 3.0*KB*N
    scaleFactor = sqrt(T0/temperature)
    for n = 1:N
        vx[n] *= scaleFactor
        vy[n] *= scaleFactor
        vz[n] *= scaleFactor
    end
    return vx, vy, vz
end

function initializeVelocityFrenkel(N, T0, vx, vy, vz)
    sumVelocity = [0.0, 0.0, 0.0] 
    sumVelocitysquare = [0.0, 0.0, 0.0]
    for n = 1:N
        vx[n] = rand()*2 - 1
        vy[n] = rand()*2 - 1
        vz[n] = rand()*2 - 1
        sumVelocity[1] +=  vx[n]
        sumVelocity[2] +=  vy[n]
        sumVelocity[3] +=  vz[n]
        sumVelocitysquare[1] += vx[n]^2
        sumVelocitysquare[2] += vy[n]^2
        sumVelocitysquare[3] += vz[n]^2
    end
    sumVelocity /= N
    sumVelocitysquare /= N
    scaleFactor = sqrt.(3.0*T0*ones(size(sumVelocitysquare))./sumVelocitysquare)
    for n = 1:N
        vx[n] = (vx[n] - sumVelocity[1])*scaleFactor[1]
        vy[n] = (vy[n] - sumVelocity[2])*scaleFactor[2]
        vz[n] = (vz[n] - sumVelocity[3])*scaleFactor[3]
    end
    return vx, vy, vz
end

    

function initializeNeighbor(N, NN, NL, x, y, z, lx, ly, lz, maxNeighbors, cutoff)
    halflx = lx*0.5
    halfly = ly*0.5
    halflz = lz*0.5
    cutoffsquare = cutoff*cutoff
    for n1 = 1:(N - 1)
        for n2 = (n1 + 1):N
            x12 = x[n2] - x[n1]
            y12 = y[n2] - y[n1]
            z12 = z[n2] - z[n1]
            x12, y12, z12 = applyMic(lx, ly, lz, halflx, halfly, halflz, x12, y12, z12)
            distancesquare = x12^2 + y12^2 + z12^2
            if distancesquare < cutoffsquare
                NL[n1*maxNeighbors + NN[n1]] = n2
                NN[n1] += 1
                NL[n2*maxNeighbors + NN[n2]] = n1
                NN[n2] += 1
            end
            if NN[n1] > maxNeighbors
                throw(DomainError(cutoff, "cutoff for neighbor list is too large!"))
            end
        end
    end
end

function applyMic(lx, ly, lz, halflx, halfly, halflz, x12, y12, z12)
    # ???
    if x12  < -halflx
        x12 += lx
    elseif x12 > halflx
        x12 -= lx
    end
    if y12 < -halfly
        y12 += ly
    elseif y12 > halfly
        y12 -= ly
    end
    if z12 < -halflz
        z12 += lz
    elseif z12 > halflz
        z12 -= lz
    end
    return x12, y12, z12
end

