import numpy as np
import matplotlib.pyplot as plt
import math

y, x = np.mgrid[-3:3:100j, -3:3:100j]
theta = math.atan(y/x)
U = 5*cos(theta)
V = -5*sin(theta)

plt.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)



plt.show()
