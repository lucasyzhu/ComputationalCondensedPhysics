using LinearAlgebra, PyPlot

function main()
    C = 1
    N = 100
    lambda = 0.5
    phi = 0
    nofb = 200
    f = zeros(N - 1, 1)
    b = range(0, stop = 1, length = nofb)
    #eigenvalues = zeros(N, nofb)
    eigenvalues = real(findeigenvalues(N, C, lambda, b, nofb, phi))
    plotHofstadterFigure(b, eigenvalues, N)
end

function findeigenvalues(N, C, lambda, b, nofb, phi)
    eigenvalues = zeros(ComplexF64, N, nofb)
    eigenEnergies = zeros(ComplexF64, 1, N)
    for i = 1:nofb
        H = constructHamiltonian(N, C, lambda, b[i], phi)
        eigenEnergies = eigvals(H)
        #println(eigenEnergies)
        eigenvalues[:, i] = abs.(eigenEnergies).^0.5
    end
    return eigenvalues
end

function constructHamiltonian(N, C, lambda, b, phi)
    forceConstant(n) = C*(1 + lambda*cos(2*pi*b*n + phi))
    H = zeros(N, N)
    temp = zeros(N, 1)
    temp[N] = forceConstant(N - 1) 
    temp[1] = forceConstant(1)
    for i in 2:N-1
        temp[i] = forceConstant(i - 1) + forceConstant(i)
    end
    for j in 1:N
        H[j, j] = temp[j]
        if j < N
            H[j, j+1] = -forceConstant(j)
            H[j+1, j] = -forceConstant(j)
        end
    end  
    return H
end

function plotHofstadterFigure(x, y, N)
    figure(figsize = (6, 4.8))
    for i = 1:N
        scatter(x, y[i, :], s = 0.5, color = "black")
    end
    xlabel("b")
    ylabel("Frequency")
    savefig("HafstadterButterfly.png", bbox_inches = "tight", dpi = 300)
end


@time main()