using LinearAlgebra, PyPlot

function testGraphene()
    nu = 0.2:1:52
    w = nu*2*pi*10.18/1000 #??
    r1 = 2
    r2 = 2.6
    mass = 1.2
    r0, L, layerSize = findr([3 24 1])
    EnergyParameter = [1393.6 430 3.5333 2.2407 1.5724e-7 
                       0.72751 38049 4.3484 -0.93 1.8 2.1]
    Transmission = findTransmission(w, r0, r1, r2, L, layerSize, 
                                    mass, EnergyParameter)
    
end

function findr(nxyz)
    r0 = [1/2 0 0;0 1/6 0;0 1/2 0;1/2 2/3 0]
    n0 = size(r0, 1)
    layerSize = nxyz[2]*n0
    N = layerSize*nxyz[1]
    a = [1.42*sqrt(3), 1.42*3, 20]
    L = a.*nxyz
    r = zeros(N, 3)
    n = 0
    for nx = 0:nxyz[1] - 1
        for ny = 0:nxyz[2] - 1
            for m = 1:n0
                n += 1
                r[n, :] = a.*([nx, ny, 0] + r0[m, :])
            end
        end
    end
    plotfigure1D(r, "1")
    return r0, L, layerSize
end

function findTransmission(w, r0, r1, r2, L, layerSize, mass, EnergyParameter)
    numofw = length(w)
    Trans = zeros(1, numofw)
    H00, H01R, H01L = findHamiltonian(r0, r1, r2, L, layerSize, mass, EnergyParameter)
    for n = 1:numofw
        wn = w[n]
        sigmaL, sigmaR = findsigma(H00, H01R, H01L, wn)
        Trans[n] = findTransmissionvsw(H00, sigmaL, sigmaR, wn, layerSize*3)
    end
    return Trans
end

function findHamiltonian(r0, r1, r2, L, layerSize, mass, EnergyParameter)
    NN1, NL1 = findNeighbor(r0, 1, r1, L)
    NN2, NL2 = findNeighbor(r0, -1, r2, L)
    H00 = zeros(layerSize*3, layerSize*3)
    H01R = zeros(layerSize*3, layerSize*3)
    for n1 = 1:layerSize
        for l = 1:NN2[n1]
            n2 = NL2[n1, l]
            index1 = ((n1 - 1)*3 + 1):(n1*3)
            index2 = ((n2 - 1)*3 + 1):(n2*3)
            Htemp = findH1(EnergyParameter, n1, n2, r0, NN1, NL1, L)
            if n2 <= layerSize
                H00[index1, index2] = Htemp
            elseif n2 <= layerSize*2
                H01R[index1, index2 - layerSize*3] = Htemp
            end
        end
    end
    H00 = H00/mass
    H01R = H01R/mass
    H01L = H01R'
    return H00, H01R, H01L
end 
function findNeighbor()
    TODO
end
function findH1()
    TODO
end

function plotfigure1D(X, mode)
    figure(figsize = (6, 4.8))
    if mode == "1"
        scatter(X[:, 1], X[:, 2], color = "blue", linewidth = 1.0)
        xlabel("x (Angstrom)")
        ylabel("y (Angstrom)")
        savefig("graphene.png", bbox_inches = "tight", dpi = 300)
    elseif mode == "2"
        plot(X[:, 1], X[:, 2], color = "red", linewidth = 1)
        xlabel("Frequency (THz)")
        ylabel("Transmission")
        savefig("grapheneTransmission.png", bbox_inches = "tight", dpi = 300)
    end
end

@time testGraphene()