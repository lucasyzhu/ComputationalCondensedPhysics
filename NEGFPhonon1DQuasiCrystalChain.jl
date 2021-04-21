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
    dw = 0.005
    eta = 1e-10
    omega = w0:dw:wmax
    N = 500
    Transmission = findTransmission(omega, N, m, C, lambda, b, phi, eta)
    plotLineFigure(omega, Transmission)
end

function findTransmission(omega, N, m, C, lambda, b, phi, delta)
    nofomega = length(omega)
    Trans = zeros(ComplexF64, nofomega, 1)
    for i in 1:nofomega
        Trans[i] = findTransmissionvsw(N, omega[i], m, C, lambda, b, phi, delta)
    end
    return Trans
end

function findTransmissionvsw(N, w, m, C, lambda, b, phi, delta)
    Dcc = findHcc(N, m, C, lambda, b, phi)
    Dcl = zeros(N, 1)
    Dcl[1] = -1/m
    Dlc = Dcl'
    Dcr = zeros(N, 1)
    Dcr[N] = -1/m
    Drc = Dcr'
    gsurface = ((w*w - 2*1/m) - w*sqrt(complex(w*w - 4*1/m)))/(2/m^2)
    surfDos = -2*w*tr(imag(gsurface)/pi)
    selfEnergyL = Dcl*gsurface*Dlc
    selfEnergyR = Dcr*gsurface*Drc
    gammaL = -2*imag(selfEnergyL)
    gammaR = -2*imag(selfEnergyR)
    Gcc = inv((w^2 + 1im*delta)*I(N) - Dcc - selfEnergyL - selfEnergyR)
    dos = -2*w*diag(imag(Gcc)/pi)
    Transvsw = sum(diag(gammaL*Gcc*gammaR*conj(Gcc)))
    return Transvsw
end

function findHcc(N, m, C, lambda, b, phi)
    forceConstant(n) = C*(1 + lambda*cos(2*pi*b*n + phi))
    Hcc = zeros(N, N)
    temp = zeros(N, 1)
    temp[N] = forceConstant(N - 1) + 1
    temp[1] = forceConstant(1) + 1
    for i in 2:N-1
        temp[i] = forceConstant(i - 1) + forceConstant(i)
    end
    for j in 1:N
        Hcc[j, j] = 1/m*temp[j]
        if j < N
            Hcc[j, j+1] = -1/m*forceConstant(j)
            Hcc[j+1, j] = -1/m*forceConstant(j)
        end
    end  
    return Hcc
end

function plotLineFigure(x, y)
    figure(figsize = (6, 4.8))
    plot(x, y, color = "red", linewidth = 1)
    xlabel("Frequency")
    ylabel("Transmission")
    savefig("QCTransmission.png", bbox_inches = "tight", dpi = 300)
end

@time main()
