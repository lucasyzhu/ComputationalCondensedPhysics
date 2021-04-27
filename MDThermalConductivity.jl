# 2021/04/27 fixed: The bug of force field computation
# The result of heat current is 0.2361 \pm 0.0268 (Run five times in total)

using LinearAlgebra, DelimitedFiles
const KB = 8.617e-5
const TIME_UNIT = 1.018e+1
const KAPPA_UNIT = 1.574e+5

function main()
    nx = 4
    ny = 4
    nz = 4
    n0 = 4
    N = n0*nx*ny*nz
    Ne = 20000
    Np = 20000
    Ns = 10
    # number of heat current list
    Nd = Int(Np/Ns)
    Nc = Int(Nd/10)
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
    timeStep = 10.0 / TIME_UNIT

    NN = zeros(Int32, N, 1)
    NL = zeros(Int32, N, maxNeighbors)
    m = zeros(Float32, N, 1)
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
    hc = zeros(Float32, 3, 1)
    # initialize mass
    for n in 1:N
        m[n] = 40.0
    end
    # Initialize position
    x, y, z = initializePosition(n0, nx, ny, nz, ax, ay, az, x, y, z)
    # Initialize velocity
    vx, vy, vz = initializeVelocity(N, T0, m, vx, vy, vz)
    # Initialize neighbor list
    NN, NL = initializeNeighbor(N, NN, NL, x, y, z, lx, ly, lz, maxNeighbors, cutoff)
    # Find force and heat current
    fx, fy, fz, hc = findForce(N, NN, NL, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc)
    # equilibration
    @time begin
    
        for step in 1:Ne
            if step % 2000 == 0
                println("The steps of equilibration is ", step)
            end
            vx, vy, vz, x, y, z = integrate(N, timeStep, m, fx, fy, fz, vx, vy, vz, x, y, z, 1)
            fx, fy, fz, hc= findForce(N, NN, NL, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc)
            vx, vy, vz, x, y, z = integrate(N, timeStep, m, fx, fy, fz, vx, vy, vz, x, y, z, 2)
            vx, vy, vz = scaleVelocity(N, T0, m, vx, vy, vz)
        
        end
    end

    @time begin
        count::Int32 = 1
        for step in 1:Np
            if step % 2000 == 0
                println("The steps of production is ", step)
            end
            vx, vy, vz, x, y, z = integrate(N, timeStep, m, fx, fy, fz, vx, vy, vz, x, y, z, 1)
            fx, fy, fz, hc= findForce(N, NN, NL, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc)
            vx, vy, vz, x, y, z = integrate(N, timeStep, m, fx, fy, fz, vx, vy, vz, x, y, z, 2)
            if step % Ns == 0
                hx[count] = hc[1]
                hy[count] = hc[2]
                hz[count] = hc[3]
                count += 1
            end
        end
    end

    findHacKappa(Nd, Nc, timeStep*Ns, T0, lx*ly*lz, hx, hy, hz)
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

function scaleVelocity(N, T0, m, vx, vy, vz)
    temperature = 0.0
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
                NN[n1] += 1
                NL[n1, NN[n1]] = n2
                NN[n2] += 1
                NL[n2, NN[n2]] = n1
            end
            if NN[n1] > maxNeighbors
                throw(DomainError(cutoff, "cutoff for neighbor list is too large!"))
            end
        end
    end
    return NN, NL
end

function applyMic(lx, ly, lz, halflx, halfly, halflz, x12, y12, z12)
    # Apply periodic boundary condition and minimum image convention 
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

function findForce(N, NN, NL, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc)
    epsilon = 1.032e-2
    sigma = 3.405
    cutoff = sigma*3.0;
    cutoffsquare = cutoff*cutoff
    sigma3 = sigma*sigma*sigma
    sigma6 = sigma3*sigma3
    sigma12 = sigma6*sigma6
    factor1 = 24.0*epsilon*sigma6
    factor2 = 48.0*epsilon*sigma12
    for i in 1:3
        hc[i] = 0
    end
    for j in 1:N
        fx[j] = 0
        fy[j] = 0
        fz[j] = 0
    end
    halflx = lx*0.5
    halfly = ly*0.5
    halflz = lz*0.5
    for n = 1:N
        for m = 1:NN[n]
            nn = NL[n, m]
            if nn < n
               continue
            end
            xij = x[nn] - x[n]
            yij = y[nn] - y[n]
            zij = z[nn] - z[n]
            xij, yij, zij = applyMic(lx, ly, lz, halflx, halfly, halflz, xij, yij, zij)
            rsquare = xij*xij + yij*yij + zij*zij
            if rsquare > cutoffsquare
               continue
            end
            r4 = rsquare*rsquare
            r8 = r4*r4
            r14 = rsquare*r4*r8
            fij = factor1/r8 - factor2/r14
            fx[n] += fij*xij
            fx[nn] -= fij*xij
            fy[n] += fij*yij
            fy[nn] -= fij*yij
            fz[n] += fij*zij
            fz[nn] -= fij*zij
            # ???
            fdotv = xij*(vx[n] + vx[nn]) + yij*(vy[n] + vy[nn]) + zij*(vz[n] + vz[nn])
            fdotv *= fij*0.5
            hc[1] -= xij*fdotv
            hc[2] -= yij*fdotv
            hc[3] -= zij*fdotv
        end
    end
    return fx, fy, fz, hc
end

function integrate(N, timeStep, m, fx, fy, fz, vx, vy, vz, x, y, z, switch)
    timeStepHalf = timeStep*0.5
    for n = 1:N
        massInv = 1.0/m[n]
        ax = fx[n]*massInv
        ay = fy[n]*massInv
        az = fz[n]*massInv
        vx[n] += ax*timeStepHalf
        vy[n] += ay*timeStepHalf
        vz[n] += az*timeStepHalf
        if switch == 1
            x[n] += vx[n]*timeStep
            y[n] += vy[n]*timeStep
            z[n] += vz[n]*timeStep
        elseif switch == 2
            continue
        end
    end
    return vx, vy, vz, x, y, z
end

function findHacKappa(Nd, Nc, dt, T0, V, hx, hy, hz)
    dtinPs = dt*TIME_UNIT/1000.0
    M = Nd - Nc
    hacx = zeros(Float32, Nc, 1)
    hacy = zeros(Float32, Nc, 1)
    hacz = zeros(Float32, Nc, 1)
    rtcx = zeros(Float32, Nc, 1)
    rtcy = zeros(Float32, Nc, 1)
    rtcz = zeros(Float32, Nc, 1)
    hacx, hacy, hacz = findHac(Nc, M, hx, hy, hz, hacx, hacy, hacz)
    factor = dt*0.5*KAPPA_UNIT / (KB*T0*T0*V)
    rtcx, rtcy, rtcz = findRtc(Nc, factor, hacx, hacy, hacz, rtcx, rtcy, rtcz)
    datafile = open("kappa.dat", "a+")
    for nc = 1:Nc
        writedlm(datafile, [nc*dtinPs hacx[nc] hacy[nc] hacz[nc] rtcx[nc] rtcy[nc] rtcz[nc]])
    end
    close(datafile)
end

function findHac(Nc, M, hx, hy, hz, hacx, hacy, hacz)
    for nc = 1:Nc
        for m = 1:M
            hacx[nc] += hx[m]*hx[m + nc]
            hacy[nc] += hy[m]*hy[m + nc]
            hacz[nc] += hz[m]*hz[m + nc]
        end
        hacx[nc] /= M
        hacy[nc] /= M
        hacz[nc] /= M
    end
    return hacx, hacy, hacz
end

function findRtc(Nc, factor, hacx, hacy, hacz, rtcx, rtcy, rtcz)
    for nc = 2:Nc
        rtcx[nc] = rtcx[nc - 1] + (hacx[nc - 1] + hacx[nc])*factor
        rtcy[nc] = rtcy[nc - 1] + (hacy[nc - 1] + hacy[nc])*factor
        rtcz[nc] = rtcz[nc - 1] + (hacz[nc - 1] + hacz[nc])*factor
    end
    return rtcx, rtcy, rtcz
end


@time main()