using LinearAlgebra, DelimitedFiles, PyPlot

x = LinRange(-1, 1, 500)
y = LinRange(-4, 4, 500)
dataFile = open("BHZsurfdos1985.dat")
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
