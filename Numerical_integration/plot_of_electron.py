import numpy as np, matplotlib.pyplot as plt

N = 100
alpha = 2
start = 0
end = 5
r = np.linspace(start,end,N)
psi = np.exp(-alpha*r)
approx_zero = 0.01

plt.plot(r,psi, label="Electron")
plt.plot([start,end], [approx_zero,approx_zero], label="approx zero")
plt.grid(); plt.legend()
plt.xlabel("Distance r"); plt.ylabel("Psi function of electron")
plt.show()
