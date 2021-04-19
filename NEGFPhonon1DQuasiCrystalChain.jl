using LinearAlgebra, PyPlot

function main()
    m = 1
    C = 1
    b = (sqrt(5) - 1)/2
    phi = 0
    lambda = 0.5
    T = 300
    w0 = 0.01
    wmax = 2.5
    nw = 310
    eta = 1e-10
    omega = w0:nw:wmax
    nofN = 100
end

function findTransmissionvsw(N, w, m, delta)
    Dcc = findHcc()
    Dcl = zeros(N, 1)
    Dcl[1] = -1/m
    Dlc = Dcl'
    Dcr = zeros(N, 1)
    Dcr[N] = -1/m
    Drc = Dcr'
    gsurface = ((w*w - 2*1/m) - w*sqrt(w*w - 4*1/m))/(2/m^2)
    surfDos = -2*w*tr(imag(gsurface)/pi)
    selfEnergyL = Dcl*gsurface*Dlc
    selfEnergyR = Dcr*gsurface*Drc
    gammaL = -2*imag(selfEnergyL)
    gammaR = -2*imag(selfEnergyR)
    Gcc = inv((w^2 + 1im*delta)*I(N) - Dcc - selfEnergyL - selfEnergyR)
    dos = -2*w*diag(imag(Gcc)/pi)
    Transvsw = tr(gammaL*G*gammaR*conj(Gcc))
    return Transvsw
end

function findHcc(N, m, C, lambda, b, phi)
    forceConstant(n) = C*(1 + lambda*cos(2*pi*b*n + phi))



    