import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('path.dat', skiprows=0)
#data1 = np.loadtxt('path1.dat', skiprows=0)
#data2 = np.loadtxt('path2.dat', skiprows=0)
#data3 = np.loadtxt('euler5.dat', skiprows=0)
#data4 = np.loadtxt('euler4.dat', skiprows=0)
#data5 = np.loadtxt('euler3.dat', skiprows=0)
#data6 = np.loadtxt('euler7.dat', skiprows=0)


x = data[:,3];
y = data[:,4];
#x1 = data1[:,3];
#y1 = data1[:,4];
#x2 = data2[:,3];
#y2 = data2[:,4];
#x3 = data3[:,3];
#y3 = data3[:,4];
#x4 = data4[:,3];
#y4 = data4[:,4];
#x5 = data5[:,3];
#y5 = data5[:,4];
#x6 = data6[:,3];
#y6 = data6[:,4];


plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle tracking in vortex')
plt.plot(x, y,'k-', clip_on=False, label = 'path')
plt.xlim(-0.5,0.5)
plt.ylim(-0.5,0.5)
#plt.plot(x1, y1,'r-', clip_on=False, label = 'path1')
#plt.plot(x2, y2,'b-', clip_on=False, label = 'path2')
#plt.plot(x3, y3,'g-', clip_on=False, label = 'Euler 5 advection')
#plt.plot(x4, y4,'y-', clip_on=False, label = 'Euler 4 advection')
#plt.plot(x5, y5,'m-', clip_on=False, label = 'Euler 3 advection')
#plt.plot(x6, y6,'c-', clip_on=False, label = 'Euler 7 advection')

#legend
legend = plt.legend(loc='upper right', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width



plt.show()
