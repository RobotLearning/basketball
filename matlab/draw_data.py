# analyze real robot data
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# 3d plot for the cartesian positions
cart_file_name = "/home/okan/basketball/data/cartesian.txt"
C = np.genfromtxt(cart_file_name)
down = 5
vec = np.arange(0,np.size(C,0),down)
left_des_pos = C[vec,0:3]
right_des_pos = C[vec,3:6]
left_act_pos = C[vec,12:15]
right_act_pos = C[vec,15:18]
fig3 = plt.figure(figsize=(12,12), dpi=300)
ax = fig3.add_subplot(111, projection='3d')
#ax.set_aspect('equal')

# PLOTTING ONLY THE FIRST TRIAL!
#idx = np.where(left_des_pos[:,0] == left_des_pos[0,0])
#print 'First trial length:', idx[0]
#idx1 = np.arange(idx[0][0],idx[0][2]-1,1)
idx1 = np.arange(0,np.size(vec,0),1)
ax.plot(left_des_pos[idx1,0],left_des_pos[idx1,1],left_des_pos[idx1,2],'r',label='desired')
ax.plot(left_act_pos[idx1,0],left_act_pos[idx1,1],left_act_pos[idx1,2],'k',label='observed')
ax.plot(right_des_pos[idx1,0],right_des_pos[idx1,1],right_des_pos[idx1,2],'r',label='desired')
ax.plot(right_act_pos[idx1,0],right_act_pos[idx1,1],right_act_pos[idx1,2],'b',label='observed')
ax.scatter(left_des_pos[0,0],left_des_pos[0,1], left_des_pos[0,2],s=50,c='r')
ax.scatter(right_des_pos[0,0],right_des_pos[0,1], right_des_pos[0,2],s=50,c='r')

# Draw the ball and initial robot loc
string_len = 1.0
ball_loc = [0.0, 0.3, 1.5-string_len]
basketball_color = np.array([207,83,0])/256
ball_radius = 0.1213

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = ball_loc[0] + ball_radius * np.outer(np.cos(u), np.sin(v))
y = ball_loc[1] + ball_radius * np.outer(np.sin(u), np.sin(v))
z = ball_loc[2] + ball_radius * np.outer(np.ones(np.size(u)), np.cos(v))
#ax.plot_surface(x, y, z,  rstride=5, cstride=5, color='r', linewidth=1, shade = 0, alpha=0.5)

fig1.savefig('/home/okan/basketball/data/23.10.2017/left_arm.png')   # save the figure to file
fig2.savefig('/home/okan/basketball/data/23.10.2017/right_arm.png')   # save the figure to file
fig3.savefig('/home/okan/basketball/data/23.10.2017/cart_pos.png')   # save the figure to file
plt.show()
