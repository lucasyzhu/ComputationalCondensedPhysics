# 2021-04-13 authored by Lucas Hu
# This script computes the edge state of BHZ model by surface Green's function 

using LinearAlgebra, DelimitedFiles, PyPlot

"""
This function constructs four four by four Gamma matrices
Output: gamma1, gamma2, gamma3, gamma4
gamma1 = s_x * sigma_y
gamma2 = s_0 * sigma_z
gamma3 = s_0 * sigma_y
gamma4 = s_z * sigma_x
Note that "*" stands for direct product
"""
function GammaMatrix()
    numofDim = (4, 4)
    gamma1 = zeros(ComplexF64, numofDim)
    gamma2 = zeros(ComplexF64, numofDim)
    gamma3 = zeros(ComplexF64, numofDim)
    gamma4 = zeros(ComplexF64, numofDim)
    # Gamma1
    gamma1[1, 4] = -1im
    gamma1[2, 3] = 1im
    gamma1[3, 2] = -1im
    gamma1[4, 1] = 1im
    # Gamma2
    gamma2[1, 1] = 1
    gamma2[2, 2] = -1
    gamma2[3, 3] = 1
    gamma2[4, 4] = -1
    # Gamma3
    gamma3[1, 2] = -1im
    gamma3[2, 1] = 1im
    gamma3[3, 4] = -1im
    gamma3[4, 3] = 1im
    # Gamma4
    gamma4[1, 2] = 1
    gamma4[2, 1] = 1
    gamma4[3, 4] = -1
    gamma4[4, 3] = -1
    return gamma1, gamma2, gamma3, gamma4
end

"""
This function constructs the partial Bloch Hamiltonian for k_y
Input: ky, float64 number
Output: H00, H01, 4 by 4 matrices
"""
function PartialBlochHamiltonian(ky::Float64)
    numofDim::Int64 = 4
    # Digonal terms
    H00 = zeros(ComplexF64, 4, 4)
    # Sub-digonal terms
    H01 = zeros(ComplexF64, 4, 4)
    # Gamma matrices
    gamma1 = zeros(ComplexF64, numofDim, numofDim)
    gamma2 = zeros(ComplexF64, numofDim, numofDim)
    gamma3 = zeros(ComplexF64, numofDim, numofDim)
    gamma4 = zeros(ComplexF64, numofDim, numofDim)
    # Parameters
    u::Float64 = -1.2
    m::Float64 = 0.3
    gamma1, gamma2, gamma3, gamma4 = GammaMatrix()
    # Construct the Bloch Hamiltonian
    for i in 1:numofDim
        for j in 1:numofDim
            H00[i, j] = m*gamma1[i, j] + (u + cos(ky))*gamma2[i, j] + sin(ky)*gamma3[i, j]
            H01[i, j] = 0.5*gamma2[i, j] + 0.5im*gamma4[i, j]
        end
    end
    return H00, H01
end

"""
This function constructs the iteration of surface green's function
Input: omega, ky, float64 number
Output:
"""
function IterationSGF(omega::ComplexF64, ky::Float64)
    # Parameters
    numofDim::Int64 = 4
    epsilon::Float64 = 1e-8
    error::Float64 = 1.0
    counter::Int64 = 0
    # ???
    g0 = zeros(ComplexF64, numofDim, numofDim)
    giter = zeros(ComplexF64, numofDim, numofDim)
    gdiff = zeros(ComplexF64, numofDim, numofDim)
    h00 = zeros(ComplexF64, numofDim, numofDim)
    h01 = zeros(ComplexF64, numofDim, numofDim)
    h00, h01 = PartialBlochHamiltonian(ky)
    # Surface green's function
    g0 = inv(omega*I(numofDim) - h00)
    while(error > epsilon)
        counter += 1
        giter = inv(omega*I(numofDim) - h00 - h01*g0*transpose(conj(h01)))
        gdiff = giter - g0
        g0 = giter
        error = abs(sum(gdiff))
    end
    return g0
end

function IterationSGF1985(omega::ComplexF64, ky::Float64)
    numofDim::Int32 = 4
    epsilon::Float64 = 1e-8
    h00, h01 = PartialBlochHamiltonian(ky)
    h10 = h01'
    e = omega*I(numofDim) - h00
    es = e
    a0 = h01
    b0 = h10
    counter = 0
    while maximum(abs.(b0)) > epsilon || maximum(abs.(a0)) > epsilon
        e0 = e
        es = es - b0*inv(e0)*a0
        e = e0 - a0*inv(e0)*b0 - b0*inv(e0)*a0
        a0 = a0*inv(e0)*a0
        b0 = b0*inv(e0)*b0
        counter += 1
    end
    g0 = inv(es)
    gbulk = inv(e)
    return g0, gbulk
end




function PlotDOS(filename::String, sizeofky::Int64, sizeofomega::Int64)
    figure(figsize = (6, 4.8))
    dataFile = open(filename)
    SDOS = readdlm(dataFile)
    ky = LinRange(-pi, pi, sizeofky) 
    omega = LinRange(-4, 4, sizeofomega)
    SDOSFig = pcolormesh(ky, omega, SDOS)
    xlabel("Wavevector")
    ylabel("Frequency")
    title("Edgestate")
    savefig("BHZsurfdos.png", bbox_inches = "tight", dpi = 300)
    close(dataFile)
end

function PlotDOS1985(filename::String, sizeofky::Int32, sizeofomega::Int32)
    x = LinRange(-1, 1, sizeofky)
    y = LinRange(-4, 4, sizeofomega)
    dataFile = open(filename)
    dosdata = readdlm(dataFile)
    dosdata1 = dosdata[:, 3]
    dosdata2 = dosdata[:, 4]
    surfdos = reshape(dosdata1, (500, 500))
    bulkdos = reshape(dosdata2, (500, 500))
    fig = figure(figsize = (13, 4.8))
    subplot(121)
    ax1 = pcolormesh(x, y, surfdos)
    colorbar(ax1)
    xlabel("Wavevector")
    ylabel("Frequency") 
    title("Edge DOS")
    subplot(122)
    ax2 = pcolormesh(x, y, bulkdos)
    colorbar(ax2)
    xlabel("Wavevector")
    ylabel("Frequency") 
    title("Bulk DOS")
    savefig("BHZsurfdos1985.png", bbox_inches = "tight", dpi = 300)
end

function testsurfold()
    numofDim::Int64 = 4
    numofomega::Int64 = 0
    sizeofomega::Int64 = 500
    omegaArray = LinRange(-4, 4, sizeofomega)
    sizeofky::Int64 = 500
    kyArray = LinRange(-pi, pi, sizeofky)
    g0 = zeros(ComplexF64, numofDim, numofDim)
    SDOS::ComplexF64 = 0.0
    SDOSdata = zeros(Float64, sizeofomega, sizeofky)
    eta::Float64 = 0.01
    dataFile = open("BHZsurfdos.dat", "w")
    for omega in omegaArray
        numofomega += 1
        numofky::Int64 = 0
        for ky in kyArray
            numofky += 1
            g0 = IterationSGF(omega + 1im*eta, ky)
            SDOS = -sum(g0)/pi
            SDOSdata[numofomega, numofky] = imag(SDOS)
        end
    end
    writedlm(dataFile, SDOSdata)
    close(dataFile)
    PlotDOS("BHZsurfdos.dat", sizeofky, sizeofomega)   
end


function testsurf1985()
    eta::Float32 = 0.01
    numofDim::Int32 = 4
    sizeofomega::Int32 = 500
    omegaArray = LinRange(-4, 4, sizeofomega)
    sizeofky::Int32 = 500
    kyArray = LinRange(-pi, pi, sizeofky)
    dataFile = open("BHZsurfdos1985.dat", "w")
    for ky in kyArray
        for omega in omegaArray
            g0, gbulk = IterationSGF1985(omega + 1im*eta, ky)
            surfdos = -imag(sum(g0))/pi
            bulkdos = -imag(sum(gbulk))/pi
            writedlm(dataFile, [ky/pi omega surfdos bulkdos])
        end
    end
    close(dataFile)
    PlotDOS1985("BHZsurfdos1985.dat", sizeofky, sizeofomega)
end

@time testsurf1985()


