import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('log', skiprows=0)
data1 = np.loadtxt('log1', skiprows=0)


x = data[:,3];
y = data[:,4];
x1 = data1[:,3];
y1 = data1[:,4];


plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle tracking in vortex')
#plt.xlim([cores[0],cores[3]])
#plt.ylim([time[0],time[3]])
#plt.yticks(time)
#plt.xticks(cores)
plt.plot(x, y,'b-', clip_on=False, label = 'RK4 advection')
plt.plot(x1, y1,'r-', clip_on=False, label = 'Euler advection')

#legend
legend = plt.legend(loc='center', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width



plt.show()
