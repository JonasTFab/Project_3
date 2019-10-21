import matplotlib.pyplot as plt
import numpy as np

n = np.array([100,1000,10000,100000,1000000,10000000,100000000])
improved_time = np.array([6.673e-05,0.00051678,0.00352412,0.0265326, 0.247553, 2.47919, 24.7062])
only_paralell_time = np.array([0.000104346,0.000219668,0.0012813,0.0122907,0.127039,1.22464,12.3421])/improved_time
compiler_flag_02_time = np.array([8.971e-05,0.000130941,0.0012098,0.0102255,0.0931483,0.852854,8.45107])/improved_time
compiler_flag_03_time = np.array([9.1195e-05,0.000152967,0.000951482,0.00849046,0.0880818,0.871797,8.39174])/improved_time

plt.plot(np.log10(n)[1:],(1-only_paralell_time[1:])*100,'--b*',color = 'red',label = 'Only parallelized')
plt.plot(np.log10(n)[1:],(1-compiler_flag_02_time[1:])*100,'--b^',color = 'blue', label = 'Comp flag -O2')
plt.plot(np.log10(n)[1:],(1-compiler_flag_03_time[1:])*100,'--bo',color = 'green',label = 'Comp flag -O3')
plt.xlabel('Simulations [log10(N)]',fontsize=12)
plt.ylabel('Improvement %',fontsize=12)
plt.title('Improvement through optimization and parallization',fontsize=14)
plt.grid()
plt.legend()
plt.show()
