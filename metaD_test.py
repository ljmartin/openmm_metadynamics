import numpy as np
from matplotlib import pyplot as plt
import simtk.openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


def gaussian(yx, mux, muy, sig=0.1,A=1):
    #only want to add potential in the neighborhood of the center, otherwise it takes
    #ages to evaluate thousands of zeroes
    xrange = range(int(mux)-5, int(mux)+5)
    yrange = range(int(muy)-5, int(muy)+5)

    #iterate over points, calculate gaussian-shaped potential, and add to yx
    for y in yrange:
        for x in xrange:
            #gaussian:
            val = A*np.exp(-( np.power(x - mux, 2.) / (2 * np.power(sig, 2.))+np.power(y - muy, 2.) / (2 * np.power(sig, 2.))))
            #update the tabulated function 
            yx[y][x] = yx[y][x]+val
    return yx


#system setup
n_steps=10000
system = mm.System()
system.addParticle(1)
force = mm.CustomExternalForce('10*(x-1)^2*(x+1)^2 + y^2 + z^2')
force.addParticle(0, [])
system.addForce(force)

##
##This is the tabulated function parameters
##based on potential energy function above
xsize = 30
xmin = -3
xmax = 3
ysize = 30
ymin =-7
ymax = 7


#This is the matrix storing the hills. 
#np.reshape goes in 'C' order, which means the last axis changes the fastest. Based on tabulated
#function documentation, the fastest varying dimension is 'x'.
yx = np.zeros((ysize, xsize))

tabfun = Continuous2DFunction(xsize, ysize, yx.reshape(-1), xmin,xmax, ymin,ymax)

#make and add the biasing force
cv = CustomCentroidBondForce(1, "myfunc(x1,y1)")
cv.addTabulatedFunction("myfunc", tabfun)
cv.addGroup([0])
cv.addBond([0])
system.addForce(cv)

#simulation setup
integrator = mm.LangevinIntegrator(300, 1, 0.02)
context = mm.Context(system, integrator)
context.setPositions([[0, 0, 0]])
context.setVelocitiesToTemperature(300)

#for tracking the trajectory:
traj = np.zeros((n_steps, 3))

print('progress:')
#run "metadynamics":
for i in range(n_steps):
    print('start|'+' '*round((i/n_steps)*50)+'*'+' '*round((n_steps-i)/n_steps*50)+'|end', end='\r')
    #print('step:', i, end='\r')
    
    coords = (context.getState(getPositions=True)
            .getPositions(asNumpy=True)
            ._value)
    traj[i] = coords
    
    #get x,y coords of the particle
    xcoord = round((coords[0][0]-xmin)/(xmax-xmin)*xsize)
    ycoord = round((coords[0][1]-ymin)/(ymax-ymin)*ysize)

    #make sure they're in the grid...
    if xcoord<100 and ycoord<100 and i%1==0:
        #...then get the new hills matrix 
        yx = gaussian(yx, xcoord, ycoord)
        #update the biasing force:
        tabfun.setFunctionParameters(xsize, ysize, yx.reshape(-1), xmin,xmax, ymin,ymax)
        cv.updateParametersInContext(context)
    integrator.step(5)

#import simtk.openmm as mm
def propagate(n_steps=10000):
    system = mm.System()
    system.addParticle(1)
    force = mm.CustomExternalForce('10*(x-1)^2*(x+1)^2 + y^2 + z^2')
    force.addParticle(0, [])
    system.addForce(force)
    integrator = mm.LangevinIntegrator(300, 1, 0.02)
    context = mm.Context(system, integrator)
    context.setPositions([[0, 0, 0]])
    context.setVelocitiesToTemperature(500)
    x = np.zeros((n_steps, 3))
    print('progress:')
    for i in range(n_steps):
        print('start|'+' '*round((i/n_steps)*50)+'*'+' '*round((n_steps-i)/n_steps*50)+'|end', end='\r')
        x[i] = (context.getState(getPositions=True)
                .getPositions(asNumpy=True)
                ._value)
        integrator.step(5)
    return x

traj2 = propagate()

plt.figure(figsize=(15,15))
plt.subplot(2,1, 1)
plt.imshow(yx)
plt.subplot(2, 1, 2)
plt.plot(traj[:,0])
plt.plot(traj2[:,0])
plt.show()
