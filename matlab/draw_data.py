# analyze real robot data
import numpy as np
import matplotlib.pyplot as plt

# load joints data
dt = 0.002
# THIS IS FOR SIMULATION 
file_name = "/home/okan/basketball/data/joints.txt"

# THIS IS FOR REAL ROBOT
#file_name = "/home/okan/basketball/data/16.10.2017/joints_real_left.txt"
#file_name = "/home/okan/basketball/data/21.10.2017/joints_both.txt"
M = np.genfromtxt(file_name)
t = dt * np.linspace(1,np.size(M,0),np.size(M,0))

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
fig1, axes1 = plt.subplots(7, 2, sharex = True)
fig1.suptitle("LEFT ARM", fontsize=16)
q_des = M[:,0:7]
qd_des = M[:,14:21]
q_act = M[:,28:35]
qd_act = M[:,42:49]

for i in range(0,7):
	axes1[i,0].plot(t,q_des[:,i],'r',label='desired')
	axes1[i,0].plot(t,q_act[:,i],'b',label='observed')
	axes1[i,0].set_ylabel(r'$\theta_' + str(i+1) + r'$',fontsize='x-large')
	axes1[i,0].legend(loc='upper right', shadow=True, fontsize='x-small',borderpad=0.5)
	axes1[i,1].plot(t,qd_des[:,i],'r',label='desired')
	axes1[i,1].plot(t,qd_act[:,i],'b',label='filtered')
	axes1[i,1].set_ylabel(r'$\dot{\theta}_' + str(i+1) + r'$',fontsize='x-large')
	axes1[i,1].legend(loc='upper right', shadow=True, fontsize='x-small',borderpad=0.5)
# Put a nicer background color on the legend.
#legend1.get_frame().set_facecolor('#00FFCC')
axes1[-1,0].set_xlabel('time (sec)',fontsize='x-large')
axes1[-1,1].set_xlabel('time (sec)',fontsize='x-large')

fig2, axes2 = plt.subplots(7, 2, sharex = True)
fig2.suptitle("RIGHT ARM", fontsize=16)
q_des = M[:,7:14]
qd_des = M[:,21:28]
q_act = M[:,35:42]
qd_act = M[:,49:56]

for i in range(0,7):
	axes2[i,0].plot(t,q_des[:,i],'r',label='desired')
	axes2[i,0].plot(t,q_act[:,i],'b',label='observed')
	axes2[i,0].legend(loc='upper right', shadow=True, fontsize='x-small',borderpad=0.5)
	axes2[i,0].set_ylabel(r'$\theta_' + str(i+1) + r'$',fontsize='x-large')
	axes2[i,1].plot(t,qd_des[:,i],'r',label='desired')
	axes2[i,1].plot(t,qd_act[:,i],'b',label='filtered')
	axes2[i,1].legend(loc='upper right', shadow=True, fontsize='x-small',borderpad=0.5)
	axes2[i,1].set_ylabel(r'$\dot{\theta}_' + str(i+1) + r'$',fontsize='x-large')

fig1.text(t[np.size(M,0)/2], 5, 'time (sec)', ha='center')
fig2.text(t[np.size(M,0)/2], 5, 'time (sec)', ha='center')

# Put a nicer background color on the legend.
#legend2.get_frame().set_facecolor('#00FFCC')
axes2[-1,0].set_xlabel('time (sec)',fontsize='x-large')
axes2[-1,1].set_xlabel('time (sec)',fontsize='x-large')

plt.show()
