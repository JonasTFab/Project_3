import numpy as np, matplotlib.pyplot as plt, random as r

def func(x1, y1, x2, y2):

    r1 = np.sqrt(x1*x1+y1*y1)
    r2 = np.sqrt(x2*x2+y2*y2)
    r_diff = np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

    return r1,r2,r_diff

N=10000
lamb = 3

plt.subplot(211)
x1 = np.random.uniform(-lamb,lamb,N)
x2 = np.random.uniform(-lamb,lamb,N)
y1 = np.random.uniform(-lamb,lamb,N)
y2 = np.random.uniform(-lamb,lamb,N)
plt.plot(x1,y1, ".", label="First")
plt.plot(x2,y2, ".", label="Second")



plt.subplot(212)
x = np.linspace(1,N,N)
r1, r2, r_diff = func(x1, y1, x2, y2)
plt.plot(x,r_diff,".")
plt.show()
