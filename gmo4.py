#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 17:10:56 2018

@author: yashodhan
"""
from numba import jit
import numpy as np
#import matplotlib.pyplot as plt
import rebound #rebound is an n body simulator, not available on windows(June 2018), pip3 install rebound to install
import xlrd #Module  to import Excel sheets
#from scipy.optimize import curve_fit
import pygmo 
import time


#Importing data
book=xlrd.open_workbook("tabby.xlsx")
sheet=book.sheet_by_index(0)
xdata=[]
ydata=[]
for row in range(60000,65316):
        if row%100==0:
            xdata.append(sheet.cell_value(row,0)) #xdata is the time in days, soon to be converted into years
for row in range(60000,65316):
        if row%100==0:
            ydata.append(sheet.cell_value(row,1)) #ydata is the normalised flux from Tabby's star
       
n=len(xdata)
for i in range(n):
    xdata[i]=xdata[i]/(365.25) #Converting units from days to years.
    xdata[i]=xdata[i]-3.99978
    
xdata=np.around(xdata,3)
    #Now defining the dip function. If the stars are partially obstructing each other then the first function will be used and if the smaller star in completely engulfed then we will get the second fucntion.
    # The functions give the fractional area of the bigger star obstructed by the smaller star. Hence you will need to subtract this from 1.
    
def dip(rf,rb,d):
        if d<np.abs(rf+rb):
            if d>np.abs(rf-rb):
                return (rf**2*np.arccos((d**2+rf**2-rb**2)/(2*d*rf))+rb**2*np.arccos((d**2+rb**2-rf**2)/(2*d*rb))-0.5*((-d+rf+rb)*(d+rf-rb)*(d-rf+rb)*(d+rf+rb))**0.5)/(np.pi*rb**2) #fraction of area between two overlapping circles
            else:
                if rf>rb:
                    return 1 #Complete blockade
                else:
                    return (rf/rb)**2
        else:
            return 0

def hip(rf,rm,rb,dfm,dmb,dbf):#Correction function
    if dbf<rf+rb:
        if dbf>np.abs(rf-rb):
            if dfm>np.abs(rf-rm):
                return dip(rf,rm,dfm)
            else:
                return dip(rm,rb,dmb)
        else:
            return dip(rm,rb,dmb)
        
    else:
        return 0
#@jit  
class traj:
    def fitness(self,v):
        """Generates 9 arrays of 3 particles and 3 coordinates each at n positions from t=0 to t=tmax
        Input the radii in increasing order. l1,l2 and l3 needn't be actual figures, just ratios."""
        
        #tmax=4.02590639288 
        m1=v[0]
        m2=v[1]
        m3=v[2]
        x2=v[3]
        y2=v[4]
        z2=v[5]
        x3=v[6]
        y3=v[7]
        z3=v[8]
        g2=v[9]
        h2=v[10]
        g3=v[11]
        h3=v[12]
        r1=v[13]/215
        r2=v[14]/215
        r3=v[15]/215
        l1=v[16]
        l2=v[17]
        l3=3.6751*(l1+l2) #Constraint for the minimum dip to be 0.78
        vx2=g2*x2+h2*(x3-x2*x2*x3-x2*y2*y3-x2*z2*z3)
        vy2=g2*y2+h2*(y3-y2*x2*x3*-y2*y2*y3-y2*z2*z3)
        vz2=g2*z2+h2*(z3-z2*x2*x3-z2*y2*y3-z2*z2*z3)
        vx3=g3*x2+h3*(x3-x2*x2*x3-x2*y2*y3-x2*z2*z3)
        vy3=g3*y2+h3*(y3-y2*x2*x3*-y2*y2*y3-y2*z2*z3)
        vz3=g3*z2+h3*(z3-z2*x2*x3-z2*y2*y3-z2*z2*z3)
        #Setting up the simulation
        sim=rebound.Simulation() #Starting simulation
        sim.integrator = "ias15"
        sim.units=('yr','AU','Msun') #Other units available, visit rebound documentation for more info.
        sim.add(m=m1,x=0,y=0,z=0,vx=0,vy=0,vz=0) #Adding a particle
        sim.add(m=m2,x=x2,y=y2,z=z2,vx=vx2,vy=vy2,vz=vz2)
        sim.add(m=m3,x=x3,y=y3,z=z3,vx=vx3,vy=vy3,vz=vz3)
        particles=sim.particles#3d array of particles' information
            
        A1=np.zeros(n) #Array that stores x coordinates of particle 0 in succession of time.
        B1=np.zeros(n)
        C1=np.zeros(n)
        A2=np.zeros(n)
        B2=np.zeros(n)
        C2=np.zeros(n)
        A3=np.zeros(n)
        B3=np.zeros(n)
        C3=np.zeros(n)
        D12=np.zeros(n)
        D23=np.zeros(n)
        D31=np.zeros(n)
        l=np.zeros(n) #Lightcurve array
        
        #Running the simulation
        for i in range(len(xdata)):
            sim.integrate(xdata[i],exact_finish_time=0) #xdata[i] is the ith time
            A1[i]=particles[0].x
            B1[i]=particles[0].y
            C1[i]=particles[0].z
            A2[i]=particles[1].x
            B2[i]=particles[1].y
            C2[i]=particles[1].z
            A3[i]=particles[2].x
            B3[i]=particles[2].y
            C3[i]=particles[2].z
            D12[i]=((A1[i]-A2[i])**2+(C1[i]-C2[i])**2)**0.5 #Distance between the stars in xz plane. Since we're looking in +y direction, this is all that is needed. Stars are flat disks to us.
            D23[i]=((A3[i]-A2[i])**2+(C3[i]-C2[i])**2)**0.5
            D31[i]=((A1[i]-A3[i])**2+(C1[i]-C3[i])**2)**0.5
            
        for i in range(n):
            if B1[i]<B2[i]: #Checking the order
                if B2[i]<B3[i]:
                    if B3[i]<B1[i]:
                        l[i]=100
                    else:
                        l[i]=l1+l2*(1-dip(r1,r2,D12[i]))+l3*(1-dip(r1,r3,D31[i])-dip(r2,r3,D23[i])+hip(r1,r2,r3,D12[i],D23[i],D31[i]))
                else:
                    if B3[i]<B1[i]:
                        l[i]=l3+l1*(1-dip(r3,r1,D31[i]))+l2*(1-dip(r3,r2,D23[i])-dip(r1,r2,D12[i])+hip(r3,r1,r2,D31[i],D12[i],D23[i]))#biggest star in front, middle one far behind. The smaller ones will be on either sides
                    else:
                        l[i]=l1+l3*(1-dip(r1,r3,D31[i]))+l2*(1-dip(r1,r2,D12[i])-dip(r3,r2,D23[i])+hip(r1,r3,r2,D31[i],D23[i],D12[i]))
            else:
                if B2[i]<B3[i]:
                    if B3[i]<B1[i]:
                        l[i]=l2+l3*(1-dip(r2,r3,D23[i]))+l1*(1-dip(r2,r1,D12[i])-dip(r3,r1,D31[i])+hip(r2,r3,r1,D23[i],D31[i],D12[i]))
                    else:
                        l[i]=l2+l1*(1-dip(r2,r1,D12[i]))+l3*(1-dip(r2,r3,D23[i])-dip(r1,r3,D31[i])+hip(r2,r1,r3,D12[i],D31[i],D23[i]))
                else:
                    if B3[i]<B1[i]:
                        l[i]=l3+l2*(1-dip(r3,r2,D23[i]))+l1*(1-dip(r3,r1,D31[i])-dip(r2,r1,D12[i])+hip(r3,r2,r1,D23[i],D12[i],D31[i]))
                    else:
                        l[i]=100
        
        l=l/(l1+l2+l3)
        dl=np.zeros(n)
        for i in range(n):
            dl[i]=(ydata[i]-l[i])**2
        return (sum(dl),)
    
    def get_bounds(self):           
        return ([0.01,0.5,0.5,-3,-3,-0.5,-3,-3,-0.5,-0.2,-0.2,-0.2,-0.2,0.01,0.5,5,0.01,2],[0.1,4,8,3,3,0.5,3,3,0.5,0.2,0.2,0.2,0.2,0.1,5,20,0.1,7])
    
start_time = time.time()
prob=pygmo.problem(traj())
algo=pygmo.algorithm(pygmo.pso(gen = 30))
pop = pygmo.population(prob, 50)
pop=algo.evolve(pop)
print(list(pop.champion_x))
print("--- %s seconds ---" % (time.time() - start_time))