import math
import numpy as np
import pandas as pd
from os.path import exists
import os


N=234941	# NUmber of cells
E=116800	# Number of elements

param=np.zeros((E,6))	# ((rho0, vx0, vy0, vz0, T0, p0),...)
v=np.zeros((E,3))	# ((rho0, vx0, vy0, vz0, T0, p0),...)

directory="/home/tama3rdgen/MURI/Deepack/"
readFile="ansys_VImarple"



if(exists(directory+"20000")==0):
	os.system("mkdir 20000")

#data=np.loadtxt(directory + readFile,skiprows=1,delimiter=",")
data=pd.read_csv(directory + readFile)
print(data)

'''
print ("start to get data file: " + directory + readFile)
fr=open(directory + readFile,"r")
start=4+3*N
loopCount=-1
flag=0
for line in fr:
	loopCount+=1
	if(loopCount==start):
		flag=1
		f=open(directory+"20000/rho","w")
		f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\
		\nobject\trho;}\n\n\ndimensions [1 -3 0 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
		print ("START DENSITY")
	if(loopCount==start+E):
		flag=2
		f.write(")\n;\n\n")
		f.close()
		print ("START VELOCITY READING")
	if(loopCount==start+2*E):
		flag=3
	if(loopCount==start+3*E):
		flag=4
	if(loopCount==start+4*E):
		flag=5
		f=open(directory+"20000/T","w")
		f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\
		\nobject\tT;}\n\n\ndimensions [0 0 0 1 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
		print ("START TEMPERATURE")
	if(loopCount==start+5*E):
		flag=6
		f.write(")\n;\n\n")
		f.close()
		f=open(directory+"20000/p","w")
		f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolScalarField;\nlocation\t\"20000\";\
		\nobject\tp;}\n\n\ndimensions [1 -1 -2 0 0 0 0];\n\ninternalField nonuniform List<scalar>\n"+str(E)+"\n(\n")
		print ("START PRESSURE")
	if(loopCount==start+6*E):
		flag=0
	if(flag==1 or flag==5 or flag==6):
		f.write(line)
	if(flag==2 or flag==3 or flag==4):
		v[loopCount-start-flag*E][flag-2]=float(line)

f.write(")\n;\n\n")
f.close()

print ("WRITING VELOCITY")
f=open(directory+"20000/U","w")
f.write("FoamFile\n{\n\tversion\t2.0;\n\tformat\tascii;\nclass\tvolVectorField;\nlocation\t\"20000\";\
\nobject\tU;}\n\n\ndimensions [0 0 0 1 0 0 0];\n\ninternalField nonuniform List<vector>\n"+str(E)+"\n(\n")
for i in np.arange(E):
	f.write("("+str(v.T[0][i])+" "+str(v.T[1][i])+" "+str(v.T[2][i])+")\n")
f.write(")\n;\n\n")
f.close()
'''