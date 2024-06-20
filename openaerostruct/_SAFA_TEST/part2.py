#Safa Bakhshi
#AE 525 Homework 3
#Part 2

import numpy as np
import matplotlib.pyplot as plt
import copy
#from openaerostruct.aerodynamics.eval_mtx import _compute_finite_vortex as VORTEX

'''
#VORTEX Subroutine: Computes velocity contribution for unit vortex
def VORTEX(G, R):
    TOL  = 1e-10

    if len(G) != 3:
         raise Exception('G must be a 1x3 vector. G was size {}'.format(len(G)))
    if len(R) != 3:
         raise Exception('R must be a 1x3 vector. R was size {}'.format(len(G)))

    R = -1*R #Kroo convention: Vortex Root minus Control Point
    Gsq = np.sum(G**2)
    Rsq = np.sum(R**2) 
    R2 = R-G
    RXR2 = np.cross(R,R2)
    GXR = np.cross(G,R)
    sig1 = (np.dot(G,R))/np.linalg.norm(R)
    sig2 = (np.dot(G,R2))/np.linalg.norm(R2)
    
    E1 = np.sum(GXR**2)
    E2 = np.sum((G+R)**2)

    #Set equal to 0 if within TOL 
    if ((Rsq <= Gsq*TOL) or (E2 <= Gsq*TOL) or (E1 <= Gsq*TOL*Rsq)): #Check for colinearity (core model)
        return np.zeros(3)
    else:
        V = (sig1-sig2) / (4.0*np.pi);
        return V*np.array([(RXR2[0]/np.sum(RXR2**2)),(RXR2[1]/np.sum(RXR2**2)),(RXR2[2]/np.sum(RXR2**2))])

'''
def VORTEX(G, R):
    TOL  = 1e-10

    if len(G) != 3:
         raise Exception('G must be a 1x3 vector. G was size {}'.format(len(G)))
    if len(R) != 3:
         raise Exception('R must be a 1x3 vector. R was size {}'.format(len(G)))

    R = -1*R #Kroo convention Vortex Root minus Control Point
    Gsq = np.sum(G**2)
    Rsq = np.sum(R**2) 
    R2 = R-G
    RXR2 = np.cross(R,R2)
    GXR = np.cross(G,R)
    sig1 = (np.dot(G,R))/np.linalg.norm(R)
    sig2 = (np.dot(G,R2))/np.linalg.norm(R2)

    kfid = np.dot((G-R),(G/np.linalg.norm(G) - R/np.linalg.norm(R)))

	
    
    E1 = np.sum(GXR**2)
    E2 = np.sum((G+R)**2)

    if ((Rsq <= Gsq*TOL) or (E2 <= Gsq*TOL) or (E1 <= Gsq*TOL*Rsq)): #Check for colinearity (core model)
        return np.zeros(3)
    else:
        V = (kfid) / (4.0*np.pi)
        return V*np.array([(GXR[0]/np.sum(GXR**2)),(GXR[1]/np.sum(GXR**2)),(GXR[2]/np.sum(GXR**2))]) 

#Compute velocity induced by hairpin vortex
def hairpin(gam,x,h=2,w=2,th=np.pi/3):
	#Setup

	#Semi-infinite vortex termination point in x
	tr = -1e4

	#Segment 1
	#Origin: (-infty,-w/2,0)
	#Term: (0,-w/2,0)

	#Segment 2
	#Origin: (0,-w/2,0)
	#Term: (h*cos(th),-w/2,h*sin(th))

	#Segment 3
	#Origin: (h*cos(th),-w/2,h*sin(th))
	#Term: (h*cos(th),w/2,h*sin(th))

	#Segment 4
	#Origin: (h*cos(th),w/2,h*sin(th))
	#Term: (0,w/2,0)

	#Segment 5
	#Origin: (0,w/2,0)
	#Term: (-infty,w/2,0)


	#Init Vortex seg
	G1  = np.array([-tr,0,0])
	G2  = np.array([h*np.cos(th),0,h*np.sin(th)])
	G3  = np.array([0,w,0])
	G4  = np.array([-h*np.cos(th),0,-h*np.sin(th)])
	G5  = np.array([tr,0,0])

	#Init R (Root - Control Point)
	R1 = np.array([tr-x[0],-w/2-x[1],-x[2]])
	R2 = np.array([-x[0],-w/2-x[1],-x[2]])
	R3 = np.array([h*np.cos(th)-x[0],-w/2-x[1],h*np.sin(th)-x[2]])
	R4 = np.array([h*np.cos(th)-x[0],w/2-x[1],h*np.sin(th)-x[2]])
	R5 = np.array([-x[0],w/2-x[1],-x[2]])


	#DEBUGGING
	#print('1 V : {}'.format(VORTEX(G1,R1)))
	#print('2 V : {}'.format(VORTEX(G2,R2)))
	#print('3 V : {}'.format(VORTEX(G3,R3)))
	#print('4 V : {}'.format(VORTEX(G4,R4)))
	#print('5 V : {}'.format(VORTEX(G5,R5)))

	#Return Induced Velocity
	return gam*(VORTEX(G1,R1) + VORTEX(G2,R2) + VORTEX(G3,R3) + VORTEX(G4,R4) + VORTEX(G5,R5))


#Compute velocity induced by hairpin vortex with image
def hairpin_ground(gam,x,h=2,w=2,th=np.pi/3):
	#Setup

	#Semi-infinite vortex termination point in x
	tr = -1e4

	#Segment 1
	#Origin: (-infty,-w/2,0)
	#Term: (0,-w/2,0)

	#Segment 2
	#Origin: (0,-w/2,0)
	#Term: (h*cos(th),-w/2,h*sin(th))

	#Segment 3
	#Origin: (h*cos(th),-w/2,h*sin(th))
	#Term: (h*cos(th),w/2,h*sin(th))

	#Segment 4
	#Origin: (h*cos(th),w/2,h*sin(th))
	#Term: (0,w/2,0)

	#Segment 5
	#Origin: (0,w/2,0)
	#Term: (-infty,w/2,0)

	#Reflected segments#

	#Segment 6
	#Origin: (0,-w/2,0)
	#Term: (h*cos(th),-w/2,-h*sin(th))

	#Segment 7
	#Origin: (h*cos(th),-w/2,-h*sin(th))
	#Term: (h*cos(th),w/2,-h*sin(th))

	#Segment 8
	#Origin: (h*cos(th),-w/2,-h*sin(th))
	#Term: (0,w/2,0)


	#Init Vortex seg
	G1  = np.array([-tr,0,0])
	G2  = np.array([h*np.cos(th),0,h*np.sin(th)])
	G3  = np.array([0,w,0])
	G4  = np.array([-h*np.cos(th),0,-h*np.sin(th)])
	G5  = np.array([tr,0,0])

	G6  = np.array([h*np.cos(th),0,-h*np.sin(th)])
	G7  = np.array([0,w,0])
	G8  = np.array([-h*np.cos(th),0,h*np.sin(th)])

	#Init R (Root - Control Point)
	R1 = np.array([tr-x[0],-w/2-x[1],-x[2]])
	R2 = np.array([-x[0],-w/2-x[1],-x[2]])
	R3 = np.array([h*np.cos(th)-x[0],-w/2-x[1],h*np.sin(th)-x[2]])
	R4 = np.array([h*np.cos(th)-x[0],w/2-x[1],h*np.sin(th)-x[2]])
	R5 = np.array([-x[0],w/2-x[1],-x[2]])


	R6 = np.array([-x[0],-w/2-x[1],-x[2]])
	R7 = np.array([h*np.cos(th)-x[0],-w/2-x[1],-h*np.sin(th)-x[2]])
	R8 = np.array([h*np.cos(th)-x[0],w/2-x[1],-h*np.sin(th)-x[2]])

	#Compute real vortex contribution

	mainVortex = gam*(VORTEX(G1,R1) + VORTEX(G2,R2) + VORTEX(G3,R3) + VORTEX(G4,R4) + VORTEX(G5,R5))

	reflectVortex = -gam*(VORTEX(G1,R1) + VORTEX(G6,R6) + VORTEX(G7,R7) + VORTEX(G8,R8) + VORTEX(G5,R5))
	#Return Induced Velocity
	return mainVortex + reflectVortex



#Part a: Unit Test

gamma = 2
xtest1 = np.zeros(3)
xtest2 = np.array([0,0,1])
xtest3 = np.array([1,2,0.5])


print('Induced velocity at {} is {}'.format(xtest1,hairpin(gamma,xtest1)))
print('Induced velocity at {} is {}'.format(xtest2,hairpin(gamma,xtest2)))
print('Induced velocity at {} is {}'.format(xtest3,hairpin(gamma,xtest3)))

print('Velocity at {} is {}'.format(xtest1,hairpin(gamma,xtest1)+np.array([1,0,0])))
print('Velocity at {} is {}'.format(xtest2,hairpin(gamma,xtest2)+np.array([1,0,0])))
print('Velocity at {} is {}'.format(xtest3,hairpin(gamma,xtest3)+np.array([1,0,0])))




#Part b and c

#Compute streamlines w/ RK4
def rollup(x0,U0,xterm,dt):
	gam = 2
	iters = 0
	n = np.shape(x0)[0]
	x = np.zeros_like(x0)
	x_stream = copy.deepcopy(x0)
	f0 = np.zeros_like(x0)
	f1 = np.zeros_like(x0)
	f2 = np.zeros_like(x0)
	f3 = np.zeros_like(x0)

	while (min(x0[:,0])<xterm):
		#RK 4 Step. Loop for all streamlines
		for i in range(n):
			f0[i,:] = U0 + hairpin(gam,x0[i,:])
			f1[i,:] = U0 + hairpin(gam,x0[i,:] + (1/2)*dt*f0[i,:])
			f2[i,:] = U0 + hairpin(gam,x0[i,:] + (1/2)*dt*f1[i,:])
			f3[i,:] = U0 + hairpin(gam,x0[i,:] + dt*f2[i,:])
			x[i,:] = x0[i,:] + (dt/6)*(f0[i,:] + 2*f1[i,:] + 2*f2[i,:] + f3[i,:])
		x0 = copy.deepcopy(x)

		#Stack point history for plotting
		x_stream = np.dstack((x_stream,x))

		iters += 1

	return x_stream, iters

#Compute streamlines w/ RK4. Including image vortex for ground.
def rollup_ground(x0,U0,xterm,dt):
	gam = 2
	iters = 0
	n = np.shape(x0)[0]
	x = np.zeros_like(x0)
	x_stream = copy.deepcopy(x0)
	f0 = np.zeros_like(x0)
	f1 = np.zeros_like(x0)
	f2 = np.zeros_like(x0)
	f3 = np.zeros_like(x0)

	while (min(x0[:,0])<xterm):
		#RK 4 Step
		for i in range(n):
			f0[i,:] = U0 + hairpin_ground(gam,x0[i,:])
			f1[i,:] = U0 + hairpin_ground(gam,x0[i,:] + (1/2)*dt*f0[i,:])
			f2[i,:] = U0 + hairpin_ground(gam,x0[i,:] + (1/2)*dt*f1[i,:])
			f3[i,:] = U0 + hairpin_ground(gam,x0[i,:] + dt*f2[i,:])
			x[i,:] = x0[i,:] + (dt/6)*(f0[i,:] + 2*f1[i,:] + 2*f2[i,:] + f3[i,:])
		x0 = copy.deepcopy(x)
		x_stream = np.dstack((x_stream,x))

		iters += 1

	return x_stream, iters


def plot_streamlines(l_stream,ground=False):
	n = np.shape(l_stream)[0]
	ax = plt.figure().add_subplot(projection='3d')

	#Loop and plot all streamlines
	for i in range(n):
		x_stream = l_stream[i,0,:]
		y_stream = l_stream[i,1,:]
		z_stream = l_stream[i,2,:]
		ax.plot(x_stream,y_stream,z_stream,color='blue',label='Streamlines')
		ax.scatter(x_stream[0],y_stream[0],z_stream[0],color='red',label='Starting Points')


	#Plot Hairpin Vortex for reference
	if ground:
		ax.plot([-6,0,1,1,0,-6],[-1,-1,-1,1,1,1],[0,0,np.sqrt(3),np.sqrt(3),0,0],color='magenta',label='Hairpin vortex')
		ax.plot([-6,0,1,1,0,-6],[-1,-1,-1,1,1,1],[0,0,-np.sqrt(3),-np.sqrt(3),0,0],color='black',label='Hairpin vortex(image)')
		plt.title('Streamlines around a hairpin vortex w/ ground')
		ax.set_zlim(-2,3)
	else:
		ax.plot([-6,0,1,1,0,-6],[-1,-1,-1,1,1,1],[0,0,np.sqrt(3),np.sqrt(3),0,0],color='magenta',label='Hairpin vortex')
		plt.title('Streamlines around a hairpin vortex')
		ax.set_zlim(-0.5,3)

	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.set_xlim(-5,5)
	ax.set_ylim(-3,3)
	#plt.legend()
	
	
	plt.show()



#Initialize
U0  = np.array([1,0,0])
x0 = np.zeros([21,3])
x0[:,0] = -5
x0[:,1] = np.linspace(-1.5,1.5,21)
x0[:,2] = 1


#Part b: Hairpin
x_result, iters = rollup(x0,U0,5,.01)
plot_streamlines(x_result,ground=False)


#Part c: Hairpin
x_result_ground, iters_ground = rollup_ground(x0,U0,5,.01)
plot_streamlines(x_result_ground,ground=True)





