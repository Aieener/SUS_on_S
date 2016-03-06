# analy.py
# A python program to analyze the SUS weighting function in order to reach the following goals:
# 1. plot the weight function
# 2. generate the normalized distribution for Z=1
# 3. extrapolate the N distribution for different Zs given by the user.
# Author: Yuding Ai
# Date: 2015 Oct 23

import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def PN():
	PN = [] # a list of number distribution
	with open("P_N.txt","r") as file:
		for line in file:
			words = line.split()
			n = float(words[0]) #take the value
			PN.append(n); #append value into my WF list

	return PN

def WN():
	WN = []
	with open("SUSWeight_function.txt","r") as file:
		for line in file:
			words = line.split()
			n = float(words[0]) #take the value
			WN.append(n); #append value into my WF list

	return WN


def Pplot(PN,z):
	fig = plt.figure()	
	plt.plot(PN,'+b',markersize=3)
	Z = str(z)
	ylabel = 'P(N;Z='+ Z + ')'
	plt.ylabel('W(index)')
	# plt.xlabel('S')
	plt.xlabel('index; S = -0.5 + index/2000')

	title = Z + ').png'
	fig.savefig(title, dpi=300, bbox_inches='tight')

def enlargePplot(PN,z):
	fig = plt.figure()	
	plt.plot(PN,'+b-',markersize=3,linewidth = 0.1)
	plt.xlim(15300,16000)
	plt.ylim(0,0.006)
	Z = str(z)
	ylabel = 'P(N;Z='+ Z + ')'
	plt.ylabel(ylabel)
	plt.xlabel('N')
	title = 'ENLP(N;Z='+ Z + ').png'
	fig.savefig(title, dpi=300, bbox_inches='tight')



def exploPN(P,z):
	P_new = [] # a list of number distribution
	W = [] # the log of P
	for i in range(len(P)):
		P_new.append(P[i]*z**(i-32450))
		# if P[i] != 0:
		# 	W.append(math.log(P[i]))
		# else:
		#  W.append(-999)

	# for i in range(len(W)):

	# 	W[i] = W[i] + i*math.log(z) # update it

	# maxi = max(W)
	# if maxi > 500:
	# 	for j in range(len(W)):
	# 		W[j] = W[j]-maxi +500
	# 		P_new.append(math.exp(W[j]));
	# P_new = [float(k)/sum(P) for k in P_new]
	# print P
	return P_new


def main():


	P = PN() # take the P(N;z=1)
	W = WN()

	# Pe = exploPN(P,9.33)
	# Pplot(P,"1")
	Pplot(W,"Weighting_function")
	# Pe = exploPN(P,1.33)
	# Pplot(Pe,1.33)

	# enlargePplot(Pe,"9.33")

	# for i in range(10):
	# 	W = PN()[0] # take the original weighting function 		
	# 	t = 9.77 + 0.01*i
	# 	Pe = exploPN(W,t)
	# 	Pplot(Pe,t)
	# 	enlargePplot(Pe,t)

main()
