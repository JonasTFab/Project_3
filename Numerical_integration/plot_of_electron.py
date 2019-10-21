import numpy as np, matplotlib.pyplot as plt

"""N = 100
alpha = 2
start = 0
end = 5
r = np.linspace(start,end,N)
psi = np.exp(-alpha*r)
approx_zero = 0.01

plt.plot(r,psi, label="Electron")
#plt.plot([start,end], [approx_zero,approx_zero], label="approx zero")
plt.plot(2.3,0.01,'X', label="approx zero", color = "r")
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

"""


#################################################################################
########################Plot error as function as time###########################
#################################################################################

brute_force_time = np.array([1.9837e-05,0.00022733,0.00205109,0.0120933,0.104968,1.02758,10.2804])#time in seconds to complete the algorithm
brute_force_error = np.array([0.173965,0.0568929,0.0860584,0.0679437,0.0849495,0.0101632,0.0042974])#error

laguerre_error = np.array([0.00630838,0.00300674,0.00168394,0.00102497,0.000652004])
laguerre_time = np.array([0.565461,5.45619,30.8584,113.847,341.753])

legendre_error = np.array([0.0807253,0.0161529,0.0207566,0.0024796,0.00881355])
legenre_time = np.array([0.106745,0.856705,4.09135,15.8739,47.7655])

improved_error = np.array([0.0256453, 0.0350268, 0.00274833, 0.000468265, 0.00141248, 0.000285713,5.96489e-05])
improved_time = np.array([6.673e-05,0.00051678,0.00352412,0.0265326, 0.247553, 2.47919, 24.7062])


iterations = np.array([100,1000, 10000,100000,1000000,10000000,100000000])

plt.plot(np.log10(laguerre_time),laguerre_error,'--^',color = 'b',label = 'Laguerre')
plt.plot(np.log10(legenre_time),legendre_error,'--v',color ='r',label = 'Legendre')
plt.plot(np.log10(brute_force_time[-3:]),brute_force_error[-3:],'--o',color = "g",label = 'Brute Force MC')
plt.legend()
plt.ylabel('Absolute error')
plt.xlabel('CPU time [log10(s)]')
plt.title('Comparison of precision versus time' )
plt.show()

"""plt.plot(np.log10(iterations[:]),np.log10(improved_error[:]),"v", label = "Importance sampling",color = "r")
plt.plot(np.log10(iterations[:]),np.log10(brute_force_error[:]),"o", label = "Brute Force",color = "g")
plt.title('Absolute Error by Number of Iterations')
plt.xlabel("Iterations log10(N)")
plt.ylabel("Absolute Error log10(AbsError)")
plt.grid()
plt.legend()
plt.show()
plt.close('all')"""
#plt.plot(improved_time,improved_error)
#plt.plot(brute_force_time,brute_force_error)

#plt.show()
#plt.close()
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
