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

###################################################
#######Plot monte carlo error######################
###################################################
file = open("brute_monte_carlo_accuracy.txt", "r")
lines = file.readlines()
size = len(lines)
i = 0
N = np.zeros(size)
error = np.zeros(size)
for i in range(size):
    line = lines[i].split()
    N[i] = line[0]
    error[i] = line[2]

plt.plot(N[2:-1],error[2:-1]); plt.grid()
plt.title("Accuracy by number of simulations")
plt.xlabel("N")
plt.ylabel("Error")
plt.show()




"""K = 273.15
C_0 = 5567.11
T_0 = 0+K
T_1 = 30.5+K
T_2 = 13.36+K
m = 0.250
C_v = 4200

H_m = C_0*(T_1-T_2)/m - C_v*(T_2-T_0)
print(H_m)


gasskonstant = 8.31446261815324
trout = 43900.41/(gasskonstant*(100+K))
print(trout)"""
