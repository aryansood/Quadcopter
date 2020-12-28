import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.animation as animation
import random
import time
mass = 1
g = 9.81
al = 0.6 #arm lenght
omega = np.array([[0.],[0.],[0.]])
#omegadot = np.array([[]])
inertia = np.array([[0.0232,0,0],
                    [0,0.0232,0],
                    [0,0,0.0468]]) #tensor of inertia
phi = 0#roll
theta = math.pi/6 #pitch
psi =  math.pi/2 #yaw
dpsi = 0 #time derivative of roll pitch and yaw
dtheta = 0
dpsi = 0
kb = 0 #thrust coefficent
kd = 0 # drag coefficient
velocity = np.array([0.,0.,0.])
torque = np.array([[0.],[0.],[0.]]) #torque in the body fixed frame
input_vec = np.array([[0.],[0.],[0.],[0.]])
angle_dot = np.array([[0.],[0.],[0.]])
pos = np.array([[-0.4],[0.],[0.]])
pos1 = np.array([[0.4],[0.],[0.]])
pos11 = np.array([[0.],[-0.4],[0.]])
pos12 = np.array([[0.],[0.4],[0.]])
grav_acc = np.array([[0.],[0.],[-mass*g]])
thrust = np.array([[0.],[0.],[0.]])
acc = np.array([[0.],[0.],[0.]])
vel = np.array([[0.],[0.],[0.]])
z_d = 10
z_dot = 0
z_dot2 = 0
k_d_z = 200
k_p_z = 25
pitch_d = 0
pitch_dot = 0
k_p_p = 30#30 #10
k_d_p = 200#200
k_i_p = 0
roll_d = 0
roll_dot = 0
k_p_r = 30
k_d_r = 200
yaw_d= 0
yaw_dot = 0
k_p_y = 25
k_d_y = 70
tempod = 0
sommat = 0
x_d  = 0.1
x_dot = 0
x_d_g = 7.5
x_p_g = 0.055
y_d = 0.1
y_d_g = -7.5
y_p_g = -0.055
inv_rm = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
ol = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
def IF2B():
    #print(psi)
    global psi,theta,phi,inv_rm
    rot = np.array([[math.cos(psi)*math.cos(theta), math.sin(psi)*math.cos(theta), -math.sin(theta)],
      [-math.sin(psi)*math.cos(phi)+math.cos(psi)*math.sin(theta)*math.sin(phi),math.cos(psi)*math.cos(phi)+math.sin(psi)*math.sin(theta)*math.sin(phi),math.cos(theta)*math.sin(phi)],
      [math.sin(psi)*math.sin(phi)+math.cos(psi)*math.sin(theta)*math.cos(phi),-math.cos(psi)*math.sin(phi)+math.sin(psi)*math.sin(theta)*math.cos(phi), math.cos(theta)*math.cos(phi)]])
    inv_rm = np.linalg.inv(rot)
    return rot
def B2IF():
    global psi, theta, phi
    rot = np.array([[math.cos(-theta)*math.cos(psi),-math.cos(theta)*math.sin(psi)+math.sin(phi)*math.sin(-theta)*math.cos(psi),math.sin(-theta)*math.sin(psi)+math.cos(phi)*math.sin(-theta)*math.cos(psi)],
        [math.cos(-theta)*math.sin(psi),math.cos(phi)*math.cos(psi)+math.sin(-theta)*math.sin(phi)*math.sin(psi),-math.sin(phi)*math.cos(psi)+math.cos(phi)*math.sin(-theta)*math.sin(psi)],
        [-math.sin(-theta),math.sin(phi)*math.cos(-theta),math.cos(phi)*math.cos(-theta)]])
    return rot
#function RB21 gives the derivatives of pitch yaw and roll angle
def RB2I():
    global omega
    rotation_frame = np.array([
        [1,math.sin(phi)*math.tan(theta),math.cos(phi)*math.tan(theta)],
        [0,math.cos(phi),-math.sin(phi)],
        [0,math.sin(phi)/math.cos(theta),math.cos(phi)/math.cos(theta)]])
    return rotation_frame.dot(omega)
def omega_dot():
    global inertia
    global omega
    global torque
    global phi,theta,psi
    p1 = inertia.dot(omega)
    inertia_inv = np.linalg.inv(inertia)
    omega_1 = omega.transpose()
    p1_1 = p1.transpose()
    p2 = np.cross(omega_1,p1_1)
    p2_2 = p2.transpose()
    p3 = torque-p2_2
    p4 = inertia_inv.dot(p3)
    return p4
def control_altitude(angle_dot):
    global mass, thrust, z_d, z_dot, z_dot2, g , k_d_z , k_p_z, vel, pos, pos1
    global torque , inertia, theta, k_p_p, k_d_p, pitch_d, pitch_dot ,phi, tempod, k_i_p
    global roll_d, roll_dot, k_p_r, k_d_r, yaw_d, yaw_dot, k_p_y, k_d_y, psi, x_d, x_d_g, x_p_g
    global y_d_g, y_d, y_p_g, ol
    #rad = math.sqrt(x_d**2+y_d**2)
    #angle = math.asin(
    mid = (pos1+pos)/2
    dist = math.sqrt((x_d-mid[0][0])**2+(y_d-mid[1][0])**2)
    rot_y = np.array([[math.cos(psi),math.sin(psi)],[-math.sin(psi),math.cos(psi)]])
    l = np.array([[x_d-mid[0][0]],[y_d-mid[1][0]],[0]])
    angle = math.atan(l[1][0]/l[0][0])
    dist = l[0][0]/math.cos(angle)
    thrust[2][0] = mass*(g+z_dot2+k_d_z*(z_dot-vel[2][0])+k_p_z*(z_d-mid[2][0]))/(math.cos(theta)*math.cos(phi))
    pitch_d =(x_d_g*(0-vel[0][0])+x_p_g*(x_d-mid[0][0]))
    roll_d = (y_d_g*(0-vel[1][0])+y_p_g*(y_d-mid[1][0]))
    print(pitch_d)

    if(roll_d>=math.pi/4):
        roll_d = math.pi/4
    if(roll_d<=-math.pi/4):
        roll_d = -math.pi/4
    if(pitch_d>=math.pi/4):
        pitch_d = math.pi/4
    if(pitch_d<=-math.pi/4):
        pitch_d = -math.pi/4
    yaw_d = angle
    torque[1][0] = inertia[1][1]*(k_p_p*(pitch_d-theta)+k_d_p*(pitch_dot-angle_dot[1][0]))
    torque[0][0] = inertia[0][0]*(k_p_r*(roll_d-phi)+k_d_r*(roll_dot-angle_dot[0][0]))
    torque[2][0] = inertia[2][2]*(k_p_y*(yaw_d-psi)+k_d_y*(yaw_dot-angle_dot[2][0]))

fig = plt.figure()
ax = fig.add_axes([0,0,1,1], projection = '3d')
conta = 0
#fig1, ax1 = plt.subplots(2)
#l, = ax.plot([-0.6,0.6], [0,0],[0,0], '-', label='line 1', linewidth=2)
tempo  = time.process_time()
def PlotQuad(frame):
    #fig = plt.figure()
    #plt.clf()
    ax.clear()
    #ax.grid(False)
    #ax.xaxis.pane.fill = False
    #ax.yaxis.pane.fill = False
    #ax.zaxis.pane.fill = False
    global tempo
    global psi
    global conta
    global omega
    global phi,theta,psi,pos,pos1,pos11,pos12,acc,thrust,vel,tempod, pitch_d, inv_rm, ol
    conta+=0.1
    #print(l12)
    #print(conta)
    #IF2B()
    #ax = fig.add_axes([0,0,1,1], projection = '3d')
    ax.set_xlim3d(-5, 5)
    ax.set_ylim3d(-5,5)
    ax.set_zlim3d(0,15)
    #ax.plot([0], [0], [conta], '_', c='green', marker='o', linewidth = 0.01)
    ax.plot([pos[0][0],pos1[0][0]], [pos[1][0],pos1[1][0]],[pos[2][0],pos1[2][0]], '-', label='line 1', linewidth=2)
    ax.plot([pos11[0][0],pos12[0][0]], [pos11[1][0],pos12[1][0]], [pos11[2][0],pos12[2][0]],'-', label = 'line 2',   linewidth=2)
    l12 = omega_dot()
    tempo2 = time.process_time()
    omega = omega+l12*(tempo2-tempo)/500
    angle_dot = RB2I()
    phi = phi+angle_dot[0][0]*(tempo2-tempo)
    theta = theta+angle_dot[1][0]*(tempo2-tempo)
    psi = psi+angle_dot[2][0]*(tempo2-tempo)
    ol = IF2B()
    posm = (pos+pos1)/2
    posm1 = (pos11+pos12)/2
    pos= np.array([[-0.4],[0],[0]])
    pos1= np.array([[0.4],[0],[0]])
    pos11 = np.array([[0],[-0.4],[0]])
    pos12 = np.array([[0],[0.4],[0]])
    pos = np.dot(inv_rm,pos)
    pos1 = np.dot(inv_rm,pos1)
    pos11 = np.dot(inv_rm,pos11)
    pos12 = np.dot(inv_rm,pos12)
    pos+=posm
    pos1+=posm
    pos11+= posm1
    pos12 += posm1
    tempod = tempo2-tempo
    control_altitude(angle_dot)
    #print(ol.dot(thrust))
    #print(str(phi)+" "+str(theta))
    acc = (grav_acc+inv_rm.dot(thrust))/(500)
    vel = vel + acc*(tempo2-tempo)
    pos+=vel
    pos1+=vel
    pos11+=vel
    pos12+=vel
    tempo = tempo2
    ax.plot([posm[0][0]], [posm[1][0]], [0], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=0.6)
    ax.plot([25], [30], [0], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=0.6)

animation1 = animation.FuncAnimation(fig, PlotQuad,interval=20,repeat=True, blit = False)
plt.show()
