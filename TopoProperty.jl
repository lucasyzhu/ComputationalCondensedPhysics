# This script computes Chern number and Berry Curvature of tight-binding model 
# using an effective method from ...
# authored by Lucas Hu
# 2020-04-21 v0.0.0 The case of Haldane model is TODO 

using LinearAlgebra, PyPlot

function testChernInsulator()
    n = 100
    dk = 2*pi/n
    chernNum::ComplexF64 = 0
    BerryCurvature = Vector()
    for kx in -pi:dk:(pi - dk)
        for ky in -pi:dk:(pi - dk)
            H = findHamiltonianChern(kx, ky)
            eigenvectors = eigvecs(H)
            vector = eigenvectors[:, 2]
            Hdkx = findHamiltonianChern(kx + dk, ky)
            eigenvectors = eigvecs(Hdkx)
            vectordkx = eigenvectors[:, 2]
            Hdky = findHamiltonianChern(kx, ky + dk)
            eigenvectors = eigvecs(Hdky)
            vectordky = eigenvectors[:, 2]
            Hdkxdky = findHamiltonianChern(kx + dk, ky + dk)
            eigenvectors = eigvecs(Hdkxdky)
            vectordkxdky = eigenvectors[:, 2]
            Ux = vector'*vectordkx/abs(vector'*vectordkx) # ?? dot(conj(vector),vectordkx)
            Uy = vector'*vectordky/abs(vector'*vectordky)
            Uxy = vectordky'*vectordkxdky/abs(vectordky'*vectordkxdky)
            Uyx = vectordkx'*vectordkxdky/abs(vectordkx'*vectordkxdky)
            Fxy = log(Ux*Uyx*(1/Uxy)*(1/Uy))
            append!(BerryCurvature, real(Fxy/1im))
            chernNum += Fxy
        end
    end
    chernNum /= (2*pi*1im)
    BerryCurvatureZ = reshape(BerryCurvature, (n - 1, n - 1))
    plotFigure2D(BerryCurvatureZ, dk)
    println("Chern number of the model is ", chernNum)
    
end



function findHamiltonian(kx, ky)
    t1 = 1.0
    t2 = 1.0
    t3 = 0.5
    m = -1.0
    hx = 2*t1*cos(kx)
    hy = 2*t1*cos(ky)
    hz = m + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky)
    sigmax, sigmay, sigmaz = constructPauli()
    H = hx*sigmax + hy*sigmay + hz*sigmaz
    return H
end

function findHamiltonianChern(kx, ky)
    t = 1
    B = 1
    M = -1
    hx = t*sin(kx)
    hy = t*sin(ky)
    hz = B*(2 + M - cos(kx) -cos(ky))
    sigmax, sigmay, sigmaz = constructPauli()
    H = hx*sigmax + hy*sigmay + hz*sigmaz
    return H
end

function findHamiltonianHaldane(kx, ky)
    t1 = 1
    t2 = 1
    M = 0
    phi = pi/2
    hx = t1*(cos(sqrt(3)/3*ky) + 2*cos(1/2*kx)*cos(sqrt(3)/6*ky))
    hy = t1*(sin(sqrt(3)/3*ky) - 2*cos(1/2*kx)*sin(sqrt(3)/6*ky))
    hz = M + 2*t2*sin(phi)*(sin(kx) - 2*sin(1/2*kx)*cos(sqrt(3)/2*ky))
    sigmax, sigmay, sigmaz = constructPauli()
    H = hx*sigmax + hy*sigmay + hz*sigmaz
    return H
end


function constructPauli()
    sigmax = zeros(ComplexF64, 2, 2)
    sigmay = zeros(ComplexF64, 2, 2)
    sigmaz = zeros(ComplexF64, 2, 2)
    sigmax[1, 2] = 1.0
    sigmax[2, 1] = 1.0
    sigmay[1, 2] = -1im
    sigmay[2, 1] = 1im
    sigmaz[1, 1] = 1.0
    sigmaz[2, 2] = -1.0
    return sigmax, sigmay, sigmaz
end

function plotFigure2D(Z, dk)
    x = -pi:dk:(pi - dk)
    y = -pi:dk:(pi - dk)
    fig = pcolormesh(x, y, Z)
    colorbar(fig)
    xlabel("kx")
    ylabel("ky")
    title("BerryCurvature")
    savefig("BerryCurvature.png", bbox_inches = "tight", dpi = 300)
end

@time testChernInsulator()