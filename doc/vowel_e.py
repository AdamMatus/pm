import scipy.io as sio
import matplotlib.pyplot as plt

mat = sio.loadmat('./e.mat')
e = mat['e']

plt.plot(e)
plt.show()
