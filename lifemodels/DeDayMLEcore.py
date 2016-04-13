#! /usr/bin/python
import cPickle
import os
import sys
from multiprocessing import Pool
from multiprocessing import cpu_count
from multiprocessing import freeze_support
from itertools import chain
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
warnings.filterwarnings('ignore', 'The number of calls to function has reached maxfev')
warnings.filterwarnings('ignore', 'The occurrence of roundoff error is detected') #(numerical integration)
warnings.simplefilter("ignore") #(ubuntu 12.4LTS)
import time
from numpy import *
from scipy import stats
from scipy import optimize
from scipy import integrate
from scipy import special
import copy
#from random import sample
#from random import random
import os
import sys
import itertools
from itertools import chain
from itertools import combinations
import math
#from matplotlib.backends.backend_pdf import PdfPages
import warnings
from scipy import polyfit
from scipy.interpolate import griddata
from scipy.stats import norm
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
from pyPdf import PdfFileWriter, PdfFileReader
import numpy as np
np.seterr(all='ignore')
#from scipy.sparse.csgraph import _validation

#test data1, interval censored data. Numbers in days. Collected in Drosophila melanogaster.
atest1=array([[9.0, 10.0], [13.5, 14.5], [13.5, 14.5], [18.0, 19.0], [36.0, 37.0], [41.0, 42.0], [41.0, 42.0], [43.0, 44.0], [43.0, 44.0], [46.5, 47.5], [46.5, 47.5], [46.5, 47.5], [50.0, 51.0], [50.0, 51.0], [50.0, 51.0], [50.0, 51.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [64.0, 65.0], [64.0, 65.0], [66.0, 67.0], [66.0, 67.0], [66.0, 67.0], [66.0, 67.0], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [71.0, 72.0], [71.0, 72.0], [71.0, 72.0], [71.0, 72.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [75.5, 76.5], [75.5, 76.5], [78.0, 79.0], [80.5, 81.5], [88.0, 89.0], [3.5004999999999997, 4.5005], [3.5004999999999997, 4.5005], [3.5004999999999997, 4.5005], [3.5004999999999997, 4.5005], [11.0, 12.0], [13.5, 14.5], [23.0, 24.0], [28.0, 29.0], [32.0, 33.0], [41.0, 42.0], [43.0, 44.0], [46.5, 47.5], [46.5, 47.5], [46.5, 47.5], [50.0, 51.0], [52.0, 53.0], [52.0, 53.0], [52.0, 53.0], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [57.0, 58.0], [57.0, 58.0], [57.0, 58.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [64.0, 65.0], [64.0, 65.0], [66.0, 67.0], [66.0, 67.0], [66.0, 67.0], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [71.0, 72.0], [71.0, 72.0], [71.0, 72.0], [71.0, 72.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [75.5, 76.5], [75.5, 76.5], [78.0, 79.0], [84.5, 85.5], [11.0, 12.0], [16.0, 17.0], [23.0, 24.0], [36.0, 37.0], [36.0, 37.0], [38.5, 39.5], [38.5, 39.5], [41.0, 42.0], [46.5, 47.5], [50.0, 51.0], [50.0, 51.0], [50.0, 51.0], [50.0, 51.0], [50.0, 51.0], [52.0, 53.0], [52.0, 53.0], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [54.5, 55.5], [57.0, 58.0], [57.0, 58.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [59.0, 60.0], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [61.5, 62.5], [64.0, 65.0], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [68.5, 69.5], [71.0, 72.0], [73.0, 74.0], [73.0, 74.0], [73.0, 74.0], [80.5, 81.5], [90.0, 91.0]]).T
# lg2 -690.758735849 59.2752827002 9.10569058449
# gp2 -675.395632435 0.000589478562637 0.0753397845892
# wb2 -706.08214662 62.3924965447 3.78626428193
# lg3 -666.726793014 63.3585736359 5.81933703074 0.00319746762276
# gp3 -670.480581203 0.000192756610821 0.0921428247321 0.00209894141433
# wb3 -666.795323208 67.0450025386 6.68529112465 0.00285598803213

#test data2, time-event data. Same as the right time points in test data1
atest2=array([[10.0, 14.5, 14.5, 19.0, 37.0, 42.0, 42.0, 44.0, 44.0, 47.5, 47.5, 47.5, 51.0, 51.0, 51.0, 51.0, 60.0, 60.0, 60.0, 60.0, 60.0, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 65.0, 65.0, 67.0, 67.0, 67.0, 67.0, 69.5, 69.5, 69.5, 69.5, 69.5, 72.0, 72.0, 72.0, 72.0, 74.0, 74.0, 74.0, 74.0, 76.5, 76.5, 79.0, 81.5, 89.0, 4.5005, 4.5005, 4.5005, 4.5005, 12.0, 14.5, 24.0, 29.0, 33.0, 42.0, 44.0, 47.5, 47.5, 47.5, 51.0, 53.0, 53.0, 53.0, 55.5, 55.5, 55.5, 55.5, 55.5, 58.0, 58.0, 58.0, 60.0, 60.0, 60.0, 60.0, 60.0, 62.5, 62.5, 62.5, 62.5, 65.0, 65.0, 67.0, 67.0, 67.0, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 72.0, 72.0, 72.0, 72.0, 74.0, 74.0, 74.0, 74.0, 74.0, 76.5, 76.5, 79.0, 85.5, 12.0, 17.0, 24.0, 37.0, 37.0, 39.5, 39.5, 42.0, 47.5, 51.0, 51.0, 51.0, 51.0, 51.0, 53.0, 53.0, 55.5, 55.5, 55.5, 55.5, 55.5, 58.0, 58.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 65.0, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 72.0, 74.0, 74.0, 74.0, 81.5, 91.0]])
# lg2 -690.755334561 59.7762122227 9.10685931539
# gp2 -675.442671292 0.00056616144262 0.0753773511044
# wb2 -704.959741 62.9376755683 3.84735027035
# lg3 -666.988039979 63.8347335297 5.83175913703 0.00314518581382
# gp3 -670.683312826 0.000187385083105 0.0918753577297 0.00204800894592
# wb3 -667.052104136 67.5301949626 6.72177306096 0.00280812160446

#test data3, mixed-sensored data consisting of exact time-event obs. (treated as [t,t]) and right-censored obs. (treated as [t:+inf). From Minitab dataset 'reliable' temp80.
atest3=array([[50.0, 50.0], [60.0, 60.0], [53.0, 53.0], [40.0, 40.0], [51.0, 51.0], [99.0, inf], [35.0, 35.0], [55.0, 55.0], [74.0, 74.0], [101.0, inf], [56.0, 56.0], [45.0, 45.0], [61.0, 61.0], [92.0, inf], [73.0, inf], [51.0, 51.0], [49.0, 49.0], [24.0, 24.0], [37.0, 37.0], [31.0, 31.0], [67.0, 67.0], [62.0, 62.0], [100.0, inf], [58.0, 58.0], [46.0, 46.0], [51.0, 51.0], [27.0, 27.0], [52.0, 52.0], [48.0, 48.0], [79.0, inf], [48.0, 48.0], [67.0, 67.0], [66.0, 66.0], [27.0, 27.0], [59.0, 59.0], [48.0, 48.0], [77.0, inf], [58.0, 58.0], [51.0, 51.0], [97.0, inf], [34.0, 34.0], [79.0, inf], [91.0, inf], [41.0, 41.0], [64.0, 64.0], [81.0, inf], [105.0, inf], [84.0, inf], [54.0, 54.0], [23.0, 23.0]]).T
# lg2 -188.718140693 60.5072934568 15.9751395767
# gp2 -191.730799673 0.00436713058766 0.0254054962004
# wb2 -186.128217188 73.3444858713 2.31750651368
# lg3 nan 60.5072934568 15.9751395767 0.0
# gp3 nan 0.00436713058766 0.0254054962004 0.0
# wb3 nan 73.3444858713 2.31750651368 0.0

#the search interval for the 3rd parameter initial value: the time independent mortality. In negative log scale.

def input_data(filename, no_time, scale_fct=1.0):
	'''
	Read: 1. an input file; 2. Number of time variable (2 means interval censored data and 1 means time-event data). 3. a scaling factor
	Convert the file to: 1. a list of array of failure time/time interval. 2. all treatment combination in the dataset. 3. the header line of input file containing the definition of treatments.
	'''
	f1=open(filename)
	lines=f1.readlines()
	f1.close()
	split_lines=[item[:-1].split(' ') for item in lines]
	header=split_lines[0]
	bodydata=split_lines[1:]
	#no_time=1
	#number of time columns
	treat_combi=[' '.join(item[:-no_time]) for item in bodydata]
	time_combi=[map(float,item[-no_time:]) for item in bodydata]
	treat_set=list(set(treat_combi))
	#===to sort===
	treat_2sort=[item.split(' ') for item in treat_set]
	for i in range(len(treat_2sort[0]))[::-1]:
		treat_2sort.sort(key=lambda x: x[i])
	treat_set=[' '.join(item) for item in treat_2sort]
	#===end sort==
	return_days=[[jtem for item, jtem in zip(treat_combi, time_combi) if item==ktem] for ktem in treat_set]
	if no_time==2:
		return_days=[array(item).T/scale_fct for item in return_days]
	else:
		return_days=[array(item).reshape((1,len(item))) for item in return_days]
	header='\t'.join(header[:-no_time])
	#treat_set=[item.replace(' ','\t') for item in treat_set]
	return return_days, treat_set, header

#==========<define fibonacci grid and convert it back to a linear grid>===================================
#fibonacci series, positive and negative.
fibp=[1, 2, 3, 4, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073L, 4807526976L, 7778742049L, 12586269025L]
fibm=[-12586269025L, -7778742049L, -4807526976L, -2971215073L, -1836311903, -1134903170, -701408733, -433494437, -267914296, -165580141, -102334155, -63245986, -39088169, -24157817, -14930352, -9227465, -5702887, -3524578, -2178309, -1346269, -832040, -514229, -317811, -196418, -121393, -75025, -46368, -28657, -17711, -10946, -6765, -4181, -2584, -1597, -987, -610, -377, -233, -144, -89, -55, -34, -21, -13, -8, -5, -4, -3, -2, -1]
def fibspace(centre, range, steps=10):
	'''Analog to np.linspace. But the the differentce between steps follows Fibonacci series. 
	Attention: returns a array of 2*steps+1 long. 
	'''
	fibgrid=fibm[-steps:]+[0,]+fibp[:steps]
	fibgrid=[item*range*1.0/fibp[steps-1] for item in fibgrid]
	return array(fibgrid)+centre
def wrap_griddata(fibspaceX, fibspaceY, Z, steps=10):
	'''Convert 3D grid data (where X and Y are defined by fibspace()) to a 3D grid data with evenly spaced X-Y grid. Using linear interpolation.
	Attention: matplotlib.pyplot.contourf can only handle evenly X-Y grid
	'''
	egX=linspace(fibspaceX.min(),fibspaceX.max(),steps)
	egY=linspace(fibspaceY.min(),fibspaceY.max(),steps)
	if isnan(Z).all():
		egZ=zeros((steps, steps))+float('nan')
	else:
		goodidx=(isfinite(Z)&(Z<60))
		egZ=griddata((fibspaceX[goodidx], fibspaceY[goodidx]), Z[goodidx], (egX[None,:], egY[:,None]), method='linear')
	return egX, egY, egZ
#==========<end define fibonacci grid>====================================================================

#==========<Inital search for MLE>========================================================================
#for debug only: OLS search for starting guess of parameter vector. Doesn't work very fast.
if 1==0:
	olsdich={}
	olsdicp={}
	try:
		for item in ['wb2', 'wb3', 'lg2', 'lg3', 'gp2', 'gp3']:
			olsdich[item]=np.load('h'+item+'.npy')
			olsdicp[item]=np.load('p'+item+'.npy')
	except IOError:
		for dis_name, x1, x2 in zip([gp2,lg2,wb2],\
									[(-9,1,100), (-2.5,-0.5,30), (-2.5,-0.5,30)],\
									[(0,5,30), (-10,0,100), (0,5,100)]):
			olsdicp[dis_name.__name__]=array(zip(*[item.flatten() for item in meshgrid(linspace(*x1), linspace(*x2))]))
			olsdich[dis_name.__name__]=log(array(map(lambda x: dis_name.hzd(exp(x), array(range(0,50))*0.01), olsdicp[dis_name.__name__])))
		for dis_name, x1, x2, x3 in zip([gp3,lg3,wb3],\
									[(-9,1,100), (-2.5,-0.5,30), (-2.5,-0.5,30)],\
									[(0,5,30), (-10,0,100), (0,5,100)],\
									[(-14,0,51), (-14,0,51), (-14,0,51)]):
			rx, ry, rz=[linspace(*item) for item in x1, x2, x3]
			_base_zeros=zeros((len(rx), len(ry), len(rz)))
			rz3D=_base_zeros+rz
			ry3D=_base_zeros+ry[...,newaxis]
			rx3D=_base_zeros+rx[...,newaxis,newaxis]
			olsdicp[dis_name.__name__]=array(zip(*[item.flatten() for item in [rx3D, ry3D, rz3D]]))
			olsdich[dis_name.__name__]=log(array(map(lambda x: dis_name.hzd(exp(x), array(range(0,50))*0.01), olsdicp[dis_name.__name__])))
		for item in olsdich.keys():
			np.save('h'+item, olsdich[item])
			np.save('p'+item, olsdicp[item])

#exp(array(p3grid)) gives a search points for the 3rd parameter (such as the C term in Gomerptz-Makehamm)
p3grid=ppp=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5, 4.0, 4.5, 5, 5.75, 6.5, 7.25, 8, 9.25, 10.5, 11.75, 12.5, 14.0]

#the search grip for initial values. Search grips for 3 parameter distribution are based on MLE of relevant 2p distribution.
#this and the following guess_paraX() functions are for brutal force inital search. Very slow for models with >2 parameters.
search_dic={
	'gp2':[exp(linspace(*item)) for item in [(-9,1,100), (0,5,30)]],
	'lg2':[exp(linspace(*item)) for item in [(-3,1,30), (-10,0,100)]],
	'wb2':[exp(linspace(*item)) for item in [(-3,1,30), (0,5,100)]],
	'llg2':[exp(linspace(*item)) for item in [(-3,1,30), (0,5,100)]],
	'gp3':lambda x: [exp(linspace(*item)) for item in [(x[0]-8,x[0],17), (x[1]-1,x[1]+1,23), (-14,1.5,81)]],
	'lg3':lambda x: [exp(linspace(*item)) for item in [(x[0]-1,x[0]+1,33), (x[1]-2,x[1]+0.5,33), (-14,1.5,81)]],
	'wb3':lambda x: [exp(linspace(*item)) for item in [(x[0]-1,x[0]+1,17), (x[1]-0.5,x[1]+2.,23), (-14,1.5,81)]],
	'llg3':lambda x: [exp(linspace(*item)) for item in [(x[0]-1,x[0]+1,17), (x[1]-0.5,x[1]+3.,23), (-14,1.5,81)]],
	'gl3':lambda x: [exp(linspace(*item)) for item in [(x[0]-2,x[0]+2,41), (x[1]-2,x[1]+2,23), (-10,x[0]+2,81)]]}

def guess_para2(funcname, rx, ry, daydis, return_no=1):
	'''Initial MLE search for 2 parameter distributions.
	The number of the distribution log likelihood function is provided with 'funcname'
	rx, ry defines the grid to search.
	daydis is the survivorship data.
	returns the x,y pair that maximizes log likelihood, or return the top 'return_no' best x,y pairs.
	'''
	rp=array(zip(*[item.flatten() for item in meshgrid(rx, ry)]))
	rL=funcname(rp.T, daydis)
	if return_no==1:
		return rp[nanargmax(rL)]
	else:
		srL=rL[isfinite(rL)]
		srL.sort()
		return rp[map(lambda x: int(argwhere(rL==x)), srL[-return_no:])]

def guess_para3(funcname, rx, ry, rz, daydis, return_no=1):
	'''Initial MLE search for 3 parameter distributions.
	Similar to guess_para2().
	'''
	_base_zeros=zeros((len(rx), len(ry), len(rz)))
	rz3D=_base_zeros+rz
	ry3D=_base_zeros+ry[...,newaxis]
	rx3D=_base_zeros+rx[...,newaxis,newaxis]
	rp=array(zip(*[item.flatten() for item in [rx3D, ry3D, rz3D]]))
	rL=funcname(rp.T, daydis)
	if return_no==1:
		return rp[nanargmax(rL)]
	else:
		srL=rL[isfinite(rL)]
		srL.sort()
		return rp[map(lambda x: int(argwhere(rL==x)), srL[-return_no:])]

def guess_para4(funcname, rx, ry, rz, ri, daydis, return_no=1):
	'''Initial MLE search for 4 parameter distributions.
	Similar to guess_para2().
	'''
	_base_zeros=zeros((len(rx), len(ry), len(rz), len(ri)))
	ri3D=_base_zeros+ri
	rz3D=_base_zeros+rz[...,newaxis]
	ry3D=_base_zeros+ry[...,newaxis,newaxis]
	rx3D=_base_zeros+rx[...,newaxis,newaxis,newaxis]
	rp=array(zip(*[item.flatten() for item in [rx3D, ry3D, rz3D, ri3D]]))
	rL=funcname(rp.T, daydis)
	if return_no==1:
		return rp[nanargmax(rL)]
	else:
		srL=rL[isfinite(rL)]
		srL.sort()
		try:
			return rp[map(lambda x: int(argwhere(rL==x)), srL[-return_no:])]
		except:
			pass
def rec_search2(fucname, daydis, p0=array([1,1]), n=5, stepx=1, stepy=1):
	'''Recursive Initial MLE search for 2 parameter distributions.
	Starting with a intial parameter p0:
	1st round: Search a evently spaced grid around p0, of which range defined by n*stepx and n*stepy.
	2st round: Search a reduced grid, range (n-1)*stepx and (n-1)*stepy, around the best parameter found in 1st round.
	......->to nth round.
	Attention, works poorly for functions that changes sharply, such as Gompertz.
	'''
	list_lgrg=range(n)[::-1]
	list_p=[p0,]*(n)
	list_x=[item*stepx for item in list_lgrg]
	list_y=[item*stepy for item in list_lgrg]
	def temp_fuc(n):
		list_p[n+1]=guess_para2(fucname, \
					exp(linspace(log(list_p[n][0])-list_x[n], log(list_p[n][0])+list_x[n],25)),\
					exp(linspace(log(list_p[n][1])-list_y[n], log(list_p[n][1])+list_y[n],25)),daydis)
	map(temp_fuc, range(n-1))
	return list_p[-1]

#Attn: rec_search2 works very fast, but rec_search3 is not.
def rec_search3(fucname, daydis, p0=array([1,1,1]), n=4, stepx=1, stepy=1, stepz=1):
	'''Recursive Initial MLE search for 2 parameter distributions.
	Similar to rec_search2().
	'''
	list_lgrg=range(n)[::-1]
	list_p=[p0,]*(n)
	list_x=[item*stepx for item in list_lgrg]
	list_y=[item*stepy for item in list_lgrg]
	list_z=[item*stepz for item in list_lgrg]
	def temp_fuc(n):
		list_p[n+1]=guess_para3(fucname, \
					exp(linspace(log(list_p[n][0])-list_x[n], log(list_p[n][0])+list_x[n],21)),\
					exp(linspace(log(list_p[n][1])-list_y[n], log(list_p[n][1])+list_y[n],21)),\
					exp(linspace(log(list_p[n][2])-list_z[n], log(list_p[n][2])+list_z[n],21)), daydis)
	map(temp_fuc, range(n-1))
	return list_p[-1]
def guess_ols(in_array):
	'''Initial MLE search, by minimizing OLS of hazard function. 
	'''
	_result={}
	_return={}
	midpoints=mean(in_array, axis=0)
	max_LS=median(midpoints[isfinite(midpoints)])*5
	midpoints=midpoints/max_LS
	_result['evttime']=unique(append([0.], midpoints))
	_result['empcdf']=mean(_result['evttime']>=midpoints[...,newaxis], axis=0)
	_result['emppdf']=append([0.], diff(_result['empcdf'])/diff(_result['evttime']))
	_result['emphzd']=log(append(_result['emppdf'][:-1]/(1-_result['empcdf'][:-1]), [float(nan)]))
	for item in [gp2,lg2,wb2,gp3,lg3,wb3]:
		_return[item.__name__]=olsdicp[item.__name__][nanargmin(\
			map(lambda x: dot(x,x), \
			olsdich[item.__name__][:,(np.round_(_result['evttime'],2)*100).astype(int)[1:-1]]-_result['emphzd'][1:-1]\
			    ))]
		#print _return[item.__name__]
		#print log(item.psolver(_return[item.__name__], in_array/max_LS))
	return _return
#==========<end Initial search for MLE>====================================================================

#==========<derivatives and estimating vcov matrix from Hessian matrix>===================================
def derv_array(p, xdata=0, fucname=0, thes=0.999):#Fully vectorized version by cause problem with fisher information matrix calculation
		_lenp=len(p)
		_adjust=zeros((_lenp, 2*_lenp))
		_preadjp=_adjust+array(p)[...,newaxis]
		_step=array([1.0e-6, ]*_lenp)
		evenidx=range(0, 2*_lenp,2)
		oddidx=range(0, 2*_lenp,2)
		iterno=0
		while iterno<5:
			for i in range(_lenp):
				_adjust[i][i*2], _adjust[i][i*2+1]=(_step[i], -_step[i])
			_postadjp=_adjust+array(p)[...,newaxis]
			result=(fucname(_postadjp, xdata)-fucname(_preadjp, xdata)).reshape((2,_lenp))/vstack((_step, -_step)).T
			change_step=3**(1-((thes<(result[:,0]/result[:,1]))&((result[:,0]/result[:,1])<(1/thes))))
			if (change_step==1).all():
				break
			else:
				_step=_step/change_step
				iterno+=1
		return mean(result, axis=1)

def derv(p, x, fucname, terms=-1):
	'''Partial derivative of parameters (p) of log likelihood function (fucname) given survivorship data x.
	if terms=n, return the partial derivative of p[n], default: return all.
	Note: uses central difference method.
	'''
	result=[]
	p1=copy.deepcopy(p)
	p2=copy.deepcopy(p)
	for i in range(len(p)):
		k=2.22e-13
		if 1>abs(p[i])>2.22e-10:#2.22e-16 is the limit for float number in python
			k=p[i]/1000.0
		elif abs(p[i])<=2.22e-10:
			k=2.22e-13
		elif abs(p[i])>=1:
			k=1.0e-3
		p1[i]+=k/2
		p2[i]-=k/2
		high=fucname(p1,x)
		low=fucname(p2,x)

		if any([math.isnan(high), math.isinf(high), math.isnan(low), math.isinf(low)]):
			result.append(1.0e5-100.0*p[i])
			#this is very empirical, but I found for dMLEg/d(-loga) and dMLEg/d(-logb)
			#when a and b are small, MLE becomes NAN. Set the dMLe's to be positive
			#is to move fsolve iterations out of the xy range that will leads to inf MLE.
		else:
			if abs(low)<1.0e-5:
				for j in range(100):
					sc_fct=1.0e+8
					try:
						ratio=(high*sc_fct-low*sc_fct)/(k*sc_fct)
						result.append(ratio)
						break
					except RuntimeWarning:
						sc_fct*=100
			elif abs(high)>1.0e5:
				for j in range(100):
					sc_fct=1.0e-8
					try:
						ratio=(high*sc_fct-low*sc_fct)/(k*sc_fct)
						result.append(ratio)
						break
					except RuntimeWarning:
						sc_fct/=100.0
			else:
				ratio=(high-low)/k
				result.append(ratio)
		p1[i]-=k/2
		p2[i]+=k/2
	if terms==-1:
		return map(float, result)
	else:
		return float(result[terms])

def hessmatrix(p, x, fucname):
	'''Calculate Fisher information matrix using negative inverse of Hessian matrix.
	'''
	result=[]
	for i in range(len(p)):
		fuc=lambda p, x: derv(p, x, fucname,i)
		result.append(derv(p, x, fuc))
	try:
		result=-matrix(result)
		try:
			return result.I
		except linalg.LinAlgError:
			return array([[float('Nan')]*len(p)]*len(p))
	except ValueError:
		return array([[float('Nan')]*len(p)]*len(p))

def derv2(p, x, fucname):
	'''Alternative version of partial derivative calculation for 3 parameter log likelihood function.
	Using average forward difference method.
	'''
	stp=1e-8
	step=p+stp*array([[0, 0, 0],\
			 [0, 0, 1],\
			 [0, 1, 0],\
			 [0, 1, 1],\
			 [1, 0, 0],\
			 [1, 0, 1],\
			 [1, 1, 0],\
			 [1, 1, 1]])
	daa=[fucname(item, x) for item in step]
	daa=map(array, daa)
	hea=-array([daa[0]-daa[4]+daa[3]-daa[7]+daa[1]-daa[5]+daa[2]-daa[6],\
		   daa[0]-daa[2]+daa[5]-daa[7]+daa[1]-daa[3]+daa[4]-daa[6],\
		   daa[0]-daa[1]+daa[6]-daa[7]+daa[2]-daa[3]+daa[4]-daa[5]])/4
	return (hea/stp).flatten().tolist()

def hessmatrix2(p, x, fucname):
	'''Alternative version of Fisher information matrix calculation for 3 parameter logliklihood function.
	Using average forward difference method.
	'''
	stp=1e-8
	step=p+stp*array([[0, 0, 0],\
			 [0, 0, 1],\
			 [0, 1, 0],\
			 [0, 1, 1],\
			 [1, 0, 0],\
			 [1, 0, 1],\
			 [1, 1, 0],\
			 [1, 1, 1]])
	daa=[derv2(item, x, fucname) for item in step]
	daa=map(array, daa)
	hea=-array([daa[0]-daa[4]+daa[3]-daa[7]+daa[1]-daa[5]+daa[2]-daa[6],\
		   daa[0]-daa[2]+daa[5]-daa[7]+daa[1]-daa[3]+daa[4]-daa[6],\
		   daa[0]-daa[1]+daa[6]-daa[7]+daa[2]-daa[3]+daa[4]-daa[5]])/4
	try:
		try:
			return -matrix(hea).I*stp
		except linalg.LinAlgError:
			return matrix(zeros((3,3))+float('nan'))
	except ValueError:
		return matrix(zeros((3,3))+float('nan'))
#==========<end  derivatives and vcov matrix>=============================================================

#==========<Define mortality Distributions>===============================================================
'''A distribution should be defined as a class, following the rules:
1, inherit from the class mortfunc()
2, latexcdf defines the LATEX string for its cumulative distribution function
3, latexhzd defines the LATEX string for its hazard function
4, lataxp defines the LATEX strings for the names of the parameters
5, cdf(p, x) defines the cumulative distribution function, of parameter p and survivorship data x.
6, hzd(p, x) defines the hazard function.
mortfunc() class:
pdf(p, x), pdf=(1-cdf)*hzd
LL(), iLL(), eLL(): log likelihood function, for interval censored data and exact (time-event) data.
MLEsolver(), iMLEsolver(), eMLEsolver(): MLE finder by solving d(MLE)/dp=0, using Powell's hybrid method. Initial guess of MLE must be provided as argv p.
lgvcov(), ilgvcov(), elgvcov(): variance-covariance matrix of logarithm of parameters. !!!parameter must be pass as agrv p in log scale!!!
'''
def gumbel_hzd(x):
	return exp(-x)
def gumbel_cdf(x):
	return exp(exp(-x))
class mortfunc(object):
	@classmethod
	def pdf(cls, p, x):
		return (1-cls.cdf(p, x))*cls.hzd(p, x)
	@classmethod
	def iLL(cls, p, xpair):
		if len(p.shape)==1:
			p=p.reshape((p.size,1))
		#return sum(log(cls.cdf(p[..., newaxis], xpair[1])-cls.cdf(p[..., newaxis], xpair[0])), axis=1)
		xpair_itv=xpair[:, less(*xpair)]
		xpair_ext=xpair[:, equal(*xpair)] #to enable mixed censoring, left-censor as [0,t], right-censor as [t, +inf) and exact time-event as [t,t]
		return sum(log(cls.cdf(p[..., newaxis], xpair_itv[1])-cls.cdf(p[..., newaxis], xpair_itv[0])), axis=1)+\
		       sum(log(cls.pdf(p[..., newaxis], xpair_ext[0])), axis=1)
	@classmethod
	def eLL(cls, p, x):
		if len(p.shape)==1:
			p=p.reshape((p.size,1))
		return sum(log(cls.pdf(p[..., newaxis], x)), axis=1)
	@classmethod
	def LL(cls, p, x):
		if x.shape[0]==2:
			return cls.iLL(p, x)
		else:
			return cls.eLL(p, x)
	@classmethod
	def iMLEsolver(cls, p, xpair):
		p=exp(p)
		#_tempfuc=lambda x: derv(x, xpair, cls.iLL)
		#return list(optimize.fsolve(_tempfuc, p))
		_tempfuc=lambda x: -cls.iLL(x, xpair)
		return list(optimize.fmin(_tempfuc, p, maxfun=10000, disp=False))
	@classmethod
	def eMLEsolver(cls, p, xsingle):
		p=exp(p)
		#_tempfuc=lambda x: derv(x, xsingle, cls.eLL)
		#return list(optimize.fsolve(_tempfuc, p))
		_tempfuc=lambda x: -cls.eLL(x, xsingle)
		return list(optimize.fmin(_tempfuc, p, maxfun=10000, disp=False))
	@classmethod
	def psolver(cls, p, x):
		if x.shape[0]==2:
			return cls.iMLEsolver(p, x)
		else:
			return cls.eMLEsolver(p, x)
	@classmethod
	def ilgvcov(cls, p, xpair):#newchange
		_tempfuc=lambda p, x: cls.iLL(exp(p), x)
		hess1=hessmatrix(p, xpair, _tempfuc)
		if len(p)==3:
			hess2=hessmatrix2(p, xpair, _tempfuc)
		if (array(hess1).diagonal()>0).all() and isfinite(array(hess1).diagonal()).all():
			return hess1
		elif len(p)==3 and (array(hess2).diagonal()>0).all() and isfinite(array(hess2).diagonal()).all():
			return hess2
		else:
			return matrix(zeros((len(p),len(p)))+float('nan'))

	@classmethod
	def elgvcov(cls, p, xsingle):
		_tempfuc=lambda p, x: cls.eLL(exp(p), x)
		hess1=hessmatrix(p, xsingle, _tempfuc)
		if len(p)==3:
			hess2=hessmatrix2(p, xsingle, _tempfuc)
		if (array(hess1).diagonal()>0).all() and isfinite(array(hess1).diagonal()).all():
			return hess1
		elif len(p)==3 and (array(hess2).diagonal()>0).all() and isfinite(array(hess2).diagonal()).all():
			return hess2
		else:
			return matrix(zeros((len(p),len(p)))+float('nan'))
	@classmethod
	def lgvcov(cls, p, x):
		if x.shape[0]==2:
			return cls.ilgvcov(p, x)
		else:
			return cls.elgvcov(p, x)

class wb2(mortfunc):
	'''Weibull Distribution'''
	latexcdf=r'$1-e^{-(\frac{x}{\lambda})^{\kappa}}$'
	latexhzd=r'$\frac{\kappa}{\lambda}(\frac{x}{\lambda})^{\kappa-1}$'
	latexp=[r'$\lambda$', r'$\kappa$']      
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-exp(-(x/p[0])**p[1])
	@staticmethod
	def hzd(p, x):
		return (p[1]/p[0])*((x/p[0])**(p[1]-1))
	@staticmethod
	def logm(p, not_use=None):
		l=exp(p[0])
		k=exp(p[1])
		return log(l*special.gamma(1.0/k+1))
	@staticmethod
	def logq(p, qt=0.5):
		l=exp(p[0])
		k=exp(p[1])
		return log(l*(-log(qt))**(1.0/k))

class wb3(mortfunc):
	'''Weibull dist. with constant'''
	latexcdf=r'$1-e^{-(\frac{x}{\lambda})^{\kappa}-cx}$'
	latexhzd=r'$\frac{\kappa}{\lambda}(\frac{x}{\lambda})^{\kappa-1}+c$'
	latexp=[r'$\lambda$', r'$\kappa$', r'$c$']      
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-exp(-((x/p[0])**p[1]+p[2]*x))
	@staticmethod
	def hzd(p, x):
		return (p[1]/p[0])*((x/p[0])**(p[1]-1))+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: wb3.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-wb3.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])

class wb3m(mortfunc):
	'''3 parameter Weibull class, negative C'''
	@staticmethod
	def cdf(p, x):
		return 1-exp(-((x/p[0])**p[1]-p[2]*x))
	@staticmethod
	def hzd(p, x):
		return (p[1]/p[0])*((x/p[0])**(p[1]-1))-p[2]

class wl2(mortfunc):
	'''Weibull Distribution, alt expression'''
	@staticmethod
	def cdf(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return 1-exp(-exp(-tx))
	@staticmethod
	def hzd(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return exp(-tx)*p[1]/x #the second part is d(tx)/dx

class wl3(mortfunc):
	'''Weibull-logistic Distribution'''
	@staticmethod
	def cdf(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return 1-exp(-p[2]*tx+p[2]*log(p[2]*exp(tx)-1))
	#TODO: this will not work, must use numerical integration.
	@staticmethod
	def hzd(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return (1./(exp(tx)+(1.0/p[2])))*p[1]/x

class wl4(mortfunc):
	'''Weibull-logistic Distribution'''
	@staticmethod
	def cdf(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return 1-exp(-p[3]*tx+p[3]*log(p[3]*exp(tx)-1)-p[2]*x)
	@staticmethod
	def hzd(p, x):
		tx=(log(p[0])-log(x))*p[1]
		return(1./(exp(tx)+(1.0/p[2])))*p[1]/x+p[2]

class gp2(mortfunc):
	'''Gompertz Distribution'''
	latexcdf=r'$1-e^{-\frac{\alpha}{\beta}(e^{\beta x}-1)}$'
	latexhzd=r'$\alpha e^{\beta x}$'
	latexp=[r'$\alpha$', r'$\beta$']	
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-exp(-(p[0]/p[1])*(exp(p[1]*x)-1))
		#beta=p[1]
		#eta=p[0]
		#return 1-exp(-eta*(exp(beta*x)-1))
	@staticmethod
	def hzd(p, x):
		return p[0]*exp(p[1]*x)
		#return p[0]*p[1]*exp(p[1]*x)
	@staticmethod
	def logm(p, not_use=None):
		beta=exp(p[1])
		eta=exp(p[0]-p[1])
		#eta=exp(p[0])
		return log(-1/beta*exp(eta)*special.expi(-eta))
	@staticmethod
	def logq(p, qt=0.5):
		beta=exp(p[1])
		eta=exp(p[0]-p[1])
		#eta=exp(p[0])	   
		return log((1.0/beta)*log(1-(1.0/eta)*log(qt)))

class gp3(mortfunc):
	'''Gompertz dist. with constant'''
	latexcdf=r'$1-e^{-\frac{\alpha}{\beta}(e^{\beta x}-1)-cx}$'
	latexhzd=r'$\alpha e^{\beta x}+c$'
	latexp=[r'$\alpha$', r'$\beta$', r'$c$']	
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-exp(-p[2]*x-(p[0]/p[1])*(exp(p[1]*x)-1))
		#beta=p[1]
		#eta=p[0]
		#return 1-exp(-p[2]*x-eta*(exp(beta*x)-1))	      
	@staticmethod
	def hzd(p, x):
		return p[0]*exp(p[1]*x)+p[2]
		#return p[0]*p[1]*exp(p[1]*x)+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: gp3.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-gp3.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])

class gl3(mortfunc):
	'''Gompertz-Logistic dist.'''
	latexcdf=r'$1-exp(\frac{a}{\beta d}log(\frac{d+1}{d\,exp(\beta x)+1}))$'
	latexhzd=r'$\frac{\alpha e^{\beta x}}{1+d e^{\beta x}}$'
	latexp=[r'$\alpha$', r'$\beta$', r'$d$']
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		try:
			return 1-exp(p[0]*log(1+1./p[2])/p[1]/p[2]\
				     -p[0]*log(exp(p[1]*x)+(1./p[2]))/p[1]/p[2])
		except ZeroDivisionError:
			return 1-exp(-(p[0]/p[1])*(exp(p[1]*x)-1))
	@staticmethod
	def hzd(p, x):
		return (p[0]*exp(p[1]*x)/(1.+p[2]*exp(p[1]*x)))
		#return p[0]*p[1]*exp(p[1]*x)+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: gl3.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-gl3.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])

class gl4(mortfunc):
	'''Logistic Gompertz with constant'''
	latexcdf=r'$1-exp(\frac{\alpha}{\beta d}log(\frac{d+1}{d\,exp(\beta x)+1}))\cdot \frac{1}{e^{cx}}$'
	latexhzd=r'$\frac{\alpha e^{\beta x}}{1+d e^{\beta x}}+c$'
	latexp=[r'$\alpha$', r'$\beta$', r'$c$', r'$d$']
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		try:
			return 1-exp(p[0]*log(1+1./p[3])/p[1]/p[3]\
				     -p[0]*log(exp(p[1]*x)+(1./p[3]))/p[1]/p[3]-p[2]*x)
		except ZeroDivisionError:
			return 1-exp(-p[2]*x-(p[0]/p[1])*(exp(p[1]*x)-1))
	@staticmethod
	def hzd(p, x):
		return (p[0]*exp(p[1]*x)/(1.+p[3]*exp(p[1]*x)))+p[2]
		#return p[0]*p[1]*exp(p[1]*x)+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: gl4.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-gl4.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])

class gp3m(mortfunc):
	'''3 parameter Gompertz class, negative c'''
	@staticmethod
	def cdf(p, x):
		return 1-exp(p[2]*x-(p[0]/p[1])*(exp(p[1]*x)-1))
	@staticmethod
	def hzd(p, x):
		return p[0]*exp(p[1]*x)-p[2]

class lg2(mortfunc):
	'''Logistic Distribution'''
	latexcdf=r'$1-\frac{e^{\frac{\mu-x}{s}}}{1+e^{\frac{\mu-x}{s}}}$'
	latexhzd=r'$\frac{1}{s(1+e^{\frac{\mu-x}{s}})}$'
	latexp=[r'$\mu$', r'$s$']       
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1/(1+exp((p[0]-x)/p[1]))#u is p[0] and s is p[1]
		#return 0.5+0.5*tanh((x-p[0])/2/p[1])
	@staticmethod
	def hzd(p, x):
		return 1/(1+exp((p[0]-x)/p[1]))/p[1]#note that it is proportional to cdf
		#return (0.5+0.5*tanh((x-p[0])/2/p[1]))/p[1]
	@staticmethod
	def logm(p, not_use=None):
		return p[0]
	@staticmethod
	def logq(p, qt=0.5):
		mu=exp(p[0])
		s=exp(p[1])
		return log(s*log(1.0/qt-1)+mu)

class lg3(mortfunc):
	'''Logistic dist. with constant'''
	latexcdf=r'$1-\frac{e^{\frac{\mu-x}{s}}}{1+e^{\frac{\mu-x}{s}}}\cdot \frac{1}{e^{cx}}$'
	latexhzd=r'$\frac{1}{s(1+e^{\frac{\mu-x}{s}})}+c$'
	latexp=[r'$\mu$', r'$s$', r'$c$']       
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		_temp=exp((p[0]-x)/p[1])
		return 1-_temp/(1+_temp)/exp(p[2]*x)
		#return 1-(0.5-0.5*tanh((x-p[0])/2/p[1]))*exp(-p[2]*x)
	@staticmethod
	def hzd(p, x):
		return 1/(1+exp((p[0]-x)/p[1]))/p[1]+p[2]
		#return (0.5+0.5*tanh((x-p[0])/2/p[1]))/p[1]+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: lg3.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-lg3.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])

class lg3m(mortfunc):
	'''2 parameter logistic class, negative C'''
	@staticmethod
	def cdf(p, x):
		_temp=exp((p[0]-x)/p[1])
		return 1-_temp/(1+_temp)/exp(-p[2]*x)
	@staticmethod
	def hzd(p, x):
		return 1/(1+exp((p[0]-x)/p[1]))/p[1]-p[2]

class llg2(mortfunc):
	'''Log-Logistic dist.'''
	latexcdf=r'$\frac{1}{1+(x/\alpha )^{-\beta}}$'
	latexhzd=r'$\frac{(\beta /\alpha )(x/\alpha )^{\beta -1}}{1+(x/\alpha )^{\beta }}$'
	latexp=[r'$\alpha$', r'$\beta$']
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-1.0/(1.+(x/p[0])**p[1])
		#return 1-(0.5-0.5*tanh((x-p[0])/2/p[1]))*exp(-p[2]*x)
	@staticmethod
	def hzd(p, x):
		return p[1]/p[0]*((x/p[0])**(p[1]-1))/(1.+(x/p[0])**p[1])
		#return (0.5+0.5*tanh((x-p[0])/2/p[1]))/p[1]+p[2]
	@staticmethod
	def logm(p, up):
		return log(p[0]*pi/p[1]/sin(pi/p[1]))
	@staticmethod
	def logq(p, qt=0.5):
		return log(p[0])+(log(1./qt-1)/p[1])

class llg3(mortfunc):
	'''Log-Logistic dist. with constant'''
	latexcdf=r'$1-\frac{1}{1+(x/\alpha )^\beta }\cdot \frac{1}{e^{cx}}$'
	latexhzd=r'$\frac{(\beta /\alpha )(x/\alpha )^{\beta -1}}{1+(x/\alpha )^{\beta }}+c$'
	latexp=[r'$\alpha$', r'$\beta$', r'$c$']
	def __init__(self, latexcdf, latexhzd, latexp):
		self.latexcdf=latexcdf
		self.latexhzd=latexhzd
		self.latexp=latexp
	@staticmethod
	def cdf(p, x):
		return 1-1.0/(1.+(x/p[0])**p[1])/exp(x*p[2])
		#return 1-(0.5-0.5*tanh((x-p[0])/2/p[1]))*exp(-p[2]*x)
	@staticmethod
	def hzd(p, x):
		return p[1]/p[0]*((x/p[0])**(p[1]-1))/(1.+(x/p[0])**p[1])+p[2]
		#return (0.5+0.5*tanh((x-p[0])/2/p[1]))/p[1]+p[2]
	@staticmethod
	def logm(p, up):
		return log(integrate.quad(lambda x: llg3.pdf(exp(p), x)*x, 0, up)[0])
	@staticmethod
	def logq(p, qt=0.5):
		tempfuc=lambda x: (1-llg3.cdf(exp(p), x)-qt)**2
		return log(optimize.fmin(tempfuc, 0, disp=0)[0])
#==========<end Define mortality Distributions>===========================================================

#==========<Scale the data to avoid exponential overflow and find MLE of the rescaled data>===============
def confi_contf2d(result_p, dist, daydis):
	'''Calculate the 2D likelihood confidence limit.
	dist: the class of distribution (gp2, wb2, lg2, etc)
	daydis: time-to-event data.
	result_p: the MLE of dist|daydis
	returns:
	1, result_p: the MLE of dist|daydis, if a better MLE solution is found return the better one.
	2, a list. [rx, ry, V]: rx, ry defines a grid; V is the reduction of log likelihood value from the MLE, on every points of this grid.
	'''
	glc=0.5
	_count=0
	while _count<3:
		if 'gp' in dist.__name__:
			rx=linspace(log(result_p[0])-2, log(result_p[0])+2,51)
			ry=linspace(log(result_p[1])-0.75, log(result_p[1])+0.75, 101)
		else:
			rx=linspace(log(result_p[0])-glc, log(result_p[0])+glc, 51)
			ry=linspace(log(result_p[1])-glc, log(result_p[1])+glc, 51)
		gX, gY=meshgrid(rx, ry)
		gP=exp(array(zip(*[jtem.flatten() for jtem in (gX, gY)])))
		gPT=gP.T
		gLL=dist.LL(gPT, daydis)
		maxLL=nanmax(gLL)
		maxL=dist.LL(array(result_p), daydis)
		if maxLL>maxL and isfinite(maxLL):
			local_p=gP[nanargmax(gLL)]
			result_p=dist.psolver(log(local_p), daydis)
			maxL=dist.LL(array(result_p), daydis)
			#print local_p, result_p, maxL, maxLL
			if maxLL>maxL:
				result_p=local_p
				maxL=maxLL
			_count+=1
		else:
			break
	if _count==4:
		result_p=gP[nanargmax(gLL)]
		maxL=dist.LL(array(result_p), daydis)
	return list(result_p), [rx, ry, float(maxL)-gLL.reshape(gY.shape)]

def confi_contf3d(result_p, dist, daydis):
	'''Calculate the 3D likelihood confidence limit. Similar to confi_contf2d().
	'''
	_count=0
	while _count<4:
		if dist.__name__ =='gp3':
			rx=linspace(log(result_p[0])-2, log(result_p[0])+2,43)
			ry=linspace(log(result_p[1])-0.5, log(result_p[1])+1,43)
			rz=linspace(log(abs(result_p[2]))-2, log(result_p[2])+2,43)
		elif dist.__name__=='gl3':
			rx=linspace(log(result_p[0])-2, log(result_p[0])+2,43)
			ry=linspace(log(result_p[1])-0.5, log(result_p[1])+1,43)
			rz=linspace(log(abs(result_p[2]))-2, log(result_p[2])+2,43)
		else:
			rx=linspace(log(result_p[0])-0.25, log(result_p[0])+0.25,43)
			ry=linspace(log(result_p[1])-0.5, log(result_p[1])+0.5,43)
			rz=linspace(log(abs(result_p[2]))-1, log(abs(result_p[2]))+1,43)
		gX, gY=meshgrid(rx, ry)
		gX, gZ=meshgrid(rx, rz)
		gY2, gZ2=meshgrid(ry, rz)
		zeros1=zeros(gX.shape)
		zeros2=zeros(gY2.shape)
		gP1=exp(array(zip(*[jtem.flatten() for jtem in (gX, gY, zeros1+log(result_p[2]))])).T)
		gP2=exp(array(zip(*[jtem.flatten() for jtem in (gX, zeros1+log(result_p[1]), gZ)])).T)
		gP3=exp(array(zip(*[jtem.flatten() for jtem in (zeros2+log(result_p[0]), gY2, gZ2)])).T)
		gP=hstack((gP1,gP2,gP3))
		gLxy=dist.LL(gP1, daydis)
		gLxz=dist.LL(gP2, daydis)
		gLyz=dist.LL(gP3, daydis)
		gLL=hstack((gLxy, gLxz, gLyz))
		maxLL=nanmax(gLL)
		maxL=dist.LL(array(result_p), daydis)
		if maxLL>maxL and isfinite(maxLL):
			local_p=gP.T[nanargmax(gLL)]
			result_p=dist.psolver(log(local_p), daydis)
			maxL=dist.LL(array(result_p), daydis)
			#print local_p, result_p, maxL, maxLL
			if maxLL>maxL:
				result_p=local_p
				maxL=maxLL
			_count+=1
		else:
			break
	if _count==4:
		result_p=gP.T[nanargmax(gLL)]
		maxL=dist.LL(array(result_p), daydis)
	gLxy=gLxy.reshape(gY.shape)
	gLxz=gLxz.reshape(gZ.shape)
	gLyz=gLyz.reshape(gZ2.shape)
	return list(result_p), [rx, ry, rz, float(maxL)-gLxy, float(maxL)-gLxz, float(maxL)-gLyz]

#Attn: Don't make a class to store the result, otherwise you will get problem with multiprocessing.
def fit_2pdist(daydis, dist, resultdic, sub_daydis=None, fullresult=True):
	'''Calculate the MLE for 2p distributions.
	Arguments:
	1, daydis. Time to event data
	2, dist. A mortfunc() distribution class
	3, resultdic. A dictionary to store the result.
	Returns: (None) A new value is added to the result dictionary containing a list: [result_p, vcov, maxL, confi_ctf]
	1, result_p. MLE, a list
	2, vcov. A matrix, variance-covariance matrix
	3, maxL. Maximum likelihood value, a np.array.
	4, confi_ctf: A list of values for making 2D confidence limit contour plot.
	'''
	if sub_daydis==None:
		guess_p=rec_search2(dist.LL, daydis)
	else:
		guess_p=rec_search2(dist.LL, sub_daydis)
	result_p=dist.psolver(log(guess_p), daydis)
	result_p, confi_ctf=confi_contf2d(result_p, dist, daydis)
	vcov=dist.lgvcov(log(array(result_p)), daydis)
	maxL=dist.LL(array(result_p), daydis)
	resultdic[dist.__name__]=[result_p, vcov, maxL, confi_ctf]

def fit_3pdisto(daydis, dist, resultdic, fullresult=True):
	'''Old MLE calculate for 3p distributions, using brutal force initial parameter vector searching.
	'''
	if dist==gl3:
		ref_dist='gp2'
	else:
		ref_dist=dist.__name__.replace('3', '2')
	search_g=search_dic[dist.__name__](map(log, resultdic[ref_dist][0]))
	guess_pl=guess_para3(dist.LL, search_g[0], search_g[1], search_g[2], daydis, 10)
	resultlist=[]
	for i in range(10):
		guess_p=exp(mean(log(guess_pl[i:]), axis=0))
		result_p=dist.psolver(log(array(guess_p)), daydis)
		maxL=dist.LL(array(result_p), daydis)
		resultlist.append([result_p, maxL])
	resultlist.sort(key=lambda x: x[1])
	result_p, confi_ctf=confi_contf3d(resultlist[-1][0], dist, daydis)
	if any([item<0 for item in result_p]):
		result_p=resultdic[ref_dist][0]+[0,]
	vcov=dist.lgvcov(log(array(result_p)), daydis)
	maxL=dist.LL(array(result_p), daydis)
	resultdic[dist.__name__]=[result_p, vcov, maxL, confi_ctf]

def fit_4pdisto(daydis, dist, resultdic, fullresult=True):
	'''Old MLE calculate for 4p distributions, using brutal force initial parameter vector searching.
	'''
	a1=array(resultdic['gp3'][0])
	a2=array(resultdic['gl3'][0])
	if a1[-1]==0.:
		d3=exp(linspace(-14, 1, 10))
	else:
		d3=exp(linspace(log(a1[-1])-1.1, log(a1[-1])+1.1, 10))
	if a2[-1]==0.:
		d4=exp(linspace(-14, 1, 40))
	else:
		d4=exp(linspace(log(a2[-1])-1.2, log(a2[-1])+1.2, 40))
	a3=log(hstack((sqrt(a1*a2)[:2], a1[-1], a2[-1])))
	subd4=[d4[i*10:i*10+10] for i in range(4)]
	del d4
	guess_pl=vstack([guess_para4(dist.LL, exp(linspace(a3[0]-1.5, a3[0]+1.5, 40)),\
			    exp(linspace(a3[1]-0.5, a3[1]+0.5, 10)),\
			    d3,\
			    item,\
			    daydis,10) for item in subd4])
	guess_L=[-float(dist.LL(item, daydis)) for item in guess_pl]
	guess_pl=[X for Y, X in sorted(zip(guess_L, guess_pl))]
	resultlist=[]
	for i in range(10):
		guess_p=exp(mean(log(guess_pl[i:]), axis=0))
		result_p=dist.psolver(log(array(guess_p)), daydis)
		maxL=dist.LL(array(result_p), daydis)
		vcov=dist.lgvcov(log(array(result_p)), daydis)
		if isfinite(maxL):
			resultlist.append([result_p, vcov, maxL])
	for item in guess_pl:
		maxL=dist.LL(item, daydis)
		vcov=dist.lgvcov(log(array(item)), daydis)
		resultlist.append([item, vcov, maxL])
	for item, jtem in zip([[a1[0], a1[1], a1[2], 0.], [a2[0], a2[1], 0., a2[2]]], ('gp3', 'gl3')):
		maxL=resultdic[jtem][2]
		if jtem=='gp3':
			vcov=hstack((vstack((resultdic[jtem][1], zeros((1,3))+nan)), zeros((4,1))+nan))
		else:
			vcov=hstack((vstack((resultdic[jtem][1], zeros((1,3))+nan))[[0,1,3,2],:], zeros((4,1))+nan))[:,[0,1,3,2]]
		resultlist.append([item, vcov, maxL])
	resultlist.sort(key=lambda x: float(x[2]))
	resultlist=[item for item in resultlist if not any([jtem<0 for jtem in item[0]])]
	#print [item[1] for item in resultlist]
	resultdic[dist.__name__]=resultlist[-1]+[[],]

def neg_fuc(p, x, fucname):
	return -fucname(exp(p), x)

def fit_3pdist(daydis, dist, resultdic, sub_daydis=None, fullresult=True):
	'''.Calculate the MLE for 3p distributions. Similar to fit_2pdist()
	Must be called when the relevant 2p distributions has already been calculated.
	Using fmin() to get the initial parameter guess, which avoids the calculation of derivatives.
	Then: dist.psolver to find the MLE.
	'''
	max_day=nanmax(daydis)
	if dist==gl3:
		ref_dist='gp2'
	else:
		ref_dist=dist.__name__.replace('3', '2')
	rlist=[]
	for cvalue in (-14, -6, 1.5):
		guess_p=log(array(resultdic[ref_dist][0]+[exp(cvalue),]))
		if sub_daydis==None:
			guess_p=exp(optimize.fmin(neg_fuc, guess_p, args=(daydis, dist.LL), disp=0))
		else:
			guess_p=exp(optimize.fmin(neg_fuc, guess_p, args=(sub_daydis, dist.LL), disp=0))
		result_p=dist.psolver(log(array(guess_p)), daydis)
		maxL=dist.LL(array(result_p), daydis)
		#TODO
		#if dist==gl3:
		#	print result_p, maxL, dist.cdf(result_p, max_day*0.9)
		if 0.95>dist.cdf(result_p, max_day*0.9)>0.1:
			rlist.append((result_p, maxL))
	result_p=array(resultdic[ref_dist][0]+[0.,])
	rlist.append((result_p, resultdic[ref_dist][2]))
	rlist.sort(key=lambda x: x[-1])
	result_p, confi_ctf=confi_contf3d(rlist[-1][0], dist, daydis)
	maxL=dist.LL(array(result_p), daydis)
	vcov=dist.lgvcov(log(array(result_p)), daydis)
	if any([item<0 for item in result_p]) or maxL<resultdic[ref_dist][2] or result_p[-1]==0.:
		result_p=resultdic[ref_dist][0]+[0.,]
		junk, confi_ctf=confi_contf3d(resultdic[ref_dist][0]+[exp(-20.),], dist, daydis)
		vcov=hstack((vstack((resultdic[ref_dist][1], zeros((1,2))+nan)), zeros((3,1))+nan))
		maxL=resultdic[ref_dist][2]
	resultdic[dist.__name__]=[result_p, vcov, maxL, confi_ctf]
	#print result_p, vcov, maxL, resultdic[dist.__name__][2]

def fit_4pdist(daydis, dist, resultdic, sub_daydis=None, fullresult=True):
	'''Calculate the MLE for 4p distributions. Similar to fit_2pdist()
	Must be called when the relevant 4p distributions has already been calculated. 
	'''
	a1=array(resultdic['gp3'][0])
	a2=array(resultdic['gl3'][0])
	for item in (a1, a2):
		if item[-1]==0.:
			item[-1]=exp(-14.)
	guess_p=log(hstack((sqrt(a1*a2)[:2], a1[-1], a2[-1])))
	resultlist=[]
	if sub_daydis==None:
		guess_p=exp(optimize.fmin(neg_fuc, guess_p, args=(daydis, dist.LL), disp=0))
	else:
		guess_p=exp(optimize.fmin(neg_fuc, guess_p, args=(sub_daydis, dist.LL), disp=0))
	result_p=dist.psolver(log(array(guess_p)), daydis)
	maxL=dist.LL(array(result_p), daydis)
	vcov=dist.lgvcov(log(array(result_p)), daydis)
	if isfinite(maxL):
		resultlist.append([result_p, vcov, maxL])
	for item, jtem in zip([[a1[0], a1[1], a1[2], 0.], [a2[0], a2[1], 0., a2[2]]], ('gp3', 'gl3')):
		maxL=resultdic[jtem][2]
		if jtem=='gp3':
			vcov=hstack((vstack((resultdic[jtem][1], zeros((1,3))+nan)), zeros((4,1))+nan))
		else:
			vcov=hstack((vstack((resultdic[jtem][1], zeros((1,3))+nan))[[0,1,3,2],:], zeros((4,1))+nan))[:,[0,1,3,2]]
		resultlist.append([item, vcov, maxL])
	#if resultdic['gp3'][0][-1]==resultdic['gl3'][0][-1]:
	#	resultlist[-1][0][-1]=0.
	#	resultlist[-1][0][-2]=0.
	#else:
	resultlist.sort(key=lambda x: float(x[2]))
	#print [item[1] for item in resultlist]
	resultdic[dist.__name__]=resultlist[-1]+[[],]
	#print resultlist[-1], resultdic[dist.__name__][2]

def AICBIC(daydis, resultdic):
	'''Calculate AIC, AICc and BIC.
	'''
	for item in resultdic.keys():
		L=2*float(resultdic[item][2])
		n=daydis.shape[1]
		k=int(item[-1])
		AIC=2*k-L
		AICc=(2.*k*n/(n-k-1))-L
		BIC=k*log(n)-L
		resultdic[item].append(AIC)
		resultdic[item].append(AICc)
		resultdic[item].append(BIC)

def adjust_Max_ls(in_array):
	'''     Rescale the data to a HYPOTHETICAL maximum lifespan, solve for maximum likelihood and scale the result back to original scale.
	the HYPOTHETICAL maximum lifespan is defined as 3 times the median.
	'''
	_result={}
	midpoints=mean(in_array, axis=0)
	midpoints=midpoints[isfinite(midpoints)]
	max_LS=nanmax(in_array) #nanmax cause problems?
	_result['evttime']=unique(append([0.], midpoints))
	_result['midpoints']=midpoints
	_result['empcdf']=mean(_result['evttime']>=midpoints[...,newaxis], axis=0)
	_result['emppdf']=append([0.], diff(_result['empcdf'])/diff(_result['evttime']))
	_result['emphzd']=append(_result['emppdf'][:-1]/(1-_result['empcdf'][:-1]), [float(nan)])
	#_result['above_T']=lambda x: sum(max_LS*x<midpoints[...,newaxis])
	#_result['above_eq_T']=lambda x: sum(max_LS*x<=midpoints[...,newaxis])
	max_LS=median(midpoints[isfinite(midpoints)])*5
	#print max_LS, median(in_array)*5
	_result['fitresult']=fit_dists(in_array/5/median(midpoints[isfinite(midpoints)]))
	#for some unkown reason, in python 2.7.5, scipy0.12.0, numpy1.7.1 you can't 'in_array/max_LS' any more
	#in python 2.7.2, scipy0.10.1, numpy 1.7.1, it still works. (maccay data female II) 
	#print _result['fitresult']['gl3'][:3]
	_result['var_log(qtm)']=[]
	fitted_distributions=[lg2, gp2, wb2, lg3, gp3, wb3, llg2, llg3, gl3, gl4]
	fiited_distrinames=[item.__name__ for item in fitted_distributions]
	for item, jtem in zip(fiited_distrinames, fitted_distributions):
		if '3' in item and 'gl' not in item:
			_result['fitresult'][item][0][2]/=max_LS
			_result['fitresult'][item][3][2]-=log(max_LS)
		if 'wb' in item or 'llg' in item:
			_result['fitresult'][item][0][0]*=max_LS
			_result['fitresult'][item][3][0]+=log(max_LS)
		if 'lg' in item and 'llg' not in item:
			_result['fitresult'][item][0][0]*=max_LS
			_result['fitresult'][item][0][1]*=max_LS
			_result['fitresult'][item][3][0]+=log(max_LS)
			_result['fitresult'][item][3][1]+=log(max_LS)
		if 'gp' in item or 'gl3' in item:
			_result['fitresult'][item][0][0]/=max_LS
			_result['fitresult'][item][0][1]/=max_LS
			_result['fitresult'][item][3][0]-=log(max_LS)
			_result['fitresult'][item][3][1]-=log(max_LS)
		if 'gl4' == item:
			_result['fitresult'][item][0][0]/=max_LS
			_result['fitresult'][item][0][1]/=max_LS
			_result['fitresult'][item][0][2]/=max_LS
		#_result['fitresult'][item][1]=jtem.lgvcov(log(array(_result['fitresult'][item][0])), in_array)
		#this is redundent, variance-covariance will not change.
	for jtem, ktem in zip([lg2, gp2, wb2, llg2, gl3, lg3, lg3, gp3, gp3, wb3, wb3, llg3, llg3, gl4, gl4],
			      [False,False,False,False,False,False,True,False,True,False,True,False,True,False,True]):
		temp=[mean_qt_var(_result['fitresult'][jtem.__name__],jtem.logq, item, ktem) for item in hstack((0.01, linspace(0.05,0.95,19), 0.99))]
		temp.append(mean_qt_var(_result['fitresult'][jtem.__name__],jtem.logm, max_LS, ktem))
		_result['var_log(qtm)'].append((jtem.__name__+'\t'+str(ktem)).replace('2\tFalse','2\t'))
		_result['var_log(qtm)'].append(array([item[0] for item in temp]))
		_result['var_log(qtm)'].append(array([item[1] for item in temp]))
	if in_array.shape[0]==2 and not any(equal(*in_array)):
		pass
	else:
		for item in [lg2, gp2, wb2, lg3, gp3, wb3, llg2, llg3, gl3, gl4]:
			_result['fitresult'][item.__name__][2]=item.LL(array(_result['fitresult'][item.__name__][0]),in_array)
	for item, jtem in zip(['wb2', 'lg2', 'gp2', 'llg2', 'gp2', 'gp3'], ['wb3', 'lg3', 'gp3', 'llg3', 'gl3', 'gl4']):
		if _result['fitresult'][jtem][0][-1]==0.:
			_result['fitresult'][jtem][2]=_result['fitresult'][item][2]
	if _result['fitresult']['gl4'][0][-2]==0.:
		_result['fitresult']['gl4'][2]=_result['fitresult']['gl3'][2]
		#print _result['fitresult']['gl4'][2]
	return _result

def fit_dists(daydis, full_result=True):
	resultdic={}
	if daydis.shape[1]>300:
		sub_daydis=daydis[:, sorted(random.randint(0, daydis.shape[1]-1, size=(300)))]
	else:
		sub_daydis=None
	fit_2pdist(daydis, wb2, resultdic, sub_daydis)
	fit_2pdist(daydis, lg2, resultdic, sub_daydis)
	fit_2pdist(daydis, gp2, resultdic, sub_daydis)
	fit_3pdist(daydis, wb3, resultdic, sub_daydis)
	fit_3pdist(daydis, lg3, resultdic, sub_daydis)
	fit_3pdist(daydis, gp3, resultdic, sub_daydis)
	if full_result:
		fit_2pdist(daydis, llg2, resultdic, sub_daydis)
		fit_3pdist(daydis, llg3, resultdic, sub_daydis)
		fit_3pdist(daydis, gl3,  resultdic, sub_daydis)
		fit_4pdist(daydis, gl4,  resultdic, sub_daydis)
	try:
		for item, jtem in zip(['wb2', 'lg2', 'gp2', 'llg2', 'gp2', 'gp3'], ['wb3', 'lg3', 'gp3', 'llg3', 'gl3', 'gl4']):
			if jtem=='gl4':
				d=3
			else:
				d=2
			if float(resultdic[item][2])>float(resultdic[jtem][2]) or any([ktem<0 for ktem in resultdic[jtem][0]]):
				resultdic[jtem][1]=hstack((vstack((resultdic[item][1], zeros((1,d))+nan)), zeros((d+1,1))+nan))
				resultdic[jtem][0]=resultdic[item][0]+[0,]
				resultdic[jtem][2]=resultdic[item][2]

	except KeyError:
		pass
	AICBIC(daydis, resultdic)
	return resultdic

def mean_qt_var(fitresult, funcname, qt=0.5, drop_c=False):
	'''     Calculate the variance for log(qauntiles) or log(mean). if drop_c is True, constant mortality in mixed models are ignored.
	'''
	eqvdic={wb3.logm:wb2.logm, lg3.logm:lg2.logm, gp3.logm:gp2.logm,
		wb3.logq:wb2.logq, lg3.logq:lg2.logq, gp3.logq:gp2.logq,
		llg3.logq:llg2.logq, llg3.logm:llg2.logm,
		gl4.logq:gl3.logq, gl4.logm:gl3.logm}
	if drop_c:
		if funcname in [gl4.logm, gl4.logq]:
			partial_d=derv(log(array(fitresult[0])[[0,1,3]]), qt, eqvdic[funcname])
			mpartial_d=matrix(partial_d)
			return eqvdic[funcname](log(array(fitresult[0])[[0,1,3]]), qt), sum(array(mpartial_d.T*mpartial_d)*array(fitresult[1][[0,1,3]][:,[0,1,3]]))
		else:
			partial_d=derv(log(array(fitresult[0][:2])), qt, eqvdic[funcname])
			mpartial_d=matrix(partial_d)
			return eqvdic[funcname](log(array(fitresult[0][:2])), qt), sum(array(mpartial_d.T*mpartial_d)*array(fitresult[1][:2, :2]))
	else:
		partial_d=derv(log(array(fitresult[0])), qt, funcname)
		mpartial_d=matrix(partial_d)
		return funcname(log(array(fitresult[0])), qt), sum(array(mpartial_d.T*mpartial_d)*array(fitresult[1]))
#==========<End scale the data and estimate MLE>==========================================================

#==========<Partially fixing one or more parameters>======================================================
def guessp_fixed(plist, fixposls, fucfix=median): #TODO
	'''1, plist: a list of parameter lists
	2, Fixposls: a list of position where parameter is fixed for all classes/levels
	Returns: A intitial guess of parameters, with fixed parameters in the end.'''
	rglen_cell=range(len(plist[0]))
	p_varie=[[item[i] for item in plist] for i in rglen_cell if i not in fixposls]
	p_fixed=map(fucfix, [[item[i] for item in plist] for i in rglen_cell if i in fixposls])
	p_varie.append(p_fixed)
	return list(chain(*p_varie))
def partial_fixed(plist, dlist, dist, fixposls, neg=False):
	'''1, plist: a list of parameter values.
	2, dlist: a list of event-time data, one for each classes/levels
	3, dist: a mortfunc class
	4, fixposls: a list of position where the parameters are fixed for all classes/levels
	Returns: The log-likelihood value.
	'''
	nu_clss=len(dlist)             #number of classes/fixed levels
	p_varie=[plist[nu_clss*i:nu_clss*(i+1)] for i in range(len(plist)/nu_clss)]
	map(p_varie.insert, fixposls, [[plist[i],]*nu_clss for i in range(-len(fixposls), 0)])
	if neg:
		p_varie=[map(exp, item) for item in p_varie]
		return -sum(map(dist.LL, map(array, zip(*p_varie)), dlist))
	else:
		return sum(map(dist.LL, map(array,zip(*p_varie)), dlist))
def fit_partial_fixed(resultdic, dlist, dist, fixposls):
	'''
	1, resultdic: a diction of paramters. The parameters of mortfunc dist are stored as a key 'dist':
		Can be obtained from read_fitresult()
	2, dlist: datalist, from which resultdic is obtained.
	3, dist: a mortfunc class
	4, fixposls: the indices of parameters that are to be fixed.
	return:
	1, MLE parameter fit
	2, variance covariance matrix of log(MLE parameters)
	3, L: loglikelihood
	4, AIC
	5, AICc
	6, BIC
	7, k: number of parameter in the model
	'''
	if len(fixposls)>0:
		#todo
		#begin of dirty talk
		guess_p0=guessp_fixed(resultdic[dist.__name__], fixposls, median)
		guess_p1=guessp_fixed(resultdic[dist.__name__], fixposls, max)
		guess_p2=guessp_fixed(resultdic[dist.__name__], fixposls, min)
		L0=partial_fixed(map(log, guess_p0), dlist, dist, fixposls, True)
		L1=partial_fixed(map(log, guess_p1), dlist, dist, fixposls, True)
		L2=partial_fixed(map(log, guess_p2), dlist, dist, fixposls, True)
		guess_L=[[guess_p0, L0],[guess_p1, L1],[guess_p2, L2]]
		filter_gL=[(p,l) for p,l in ([guess_p0, L0],[guess_p1, L1],[guess_p2, L2]) if isfinite(l)]
		if len(filter_gL)>0:
			guess_L=filter_gL
		guess_L.sort(key=lambda x: x[1])
		guess_p=array(guess_L[0][0])
		guess_p=(guess_p+(guess_p==0.)*1e-19).tolist()
		#end of dirty talk
		#l_gp=log(array(guess_p))
		PFF=lambda p, x: -array(partial_fixed(exp(p),x,dist,fixposls))
		PPF=lambda p, x: array(derv(p, x, PFF))#;print 'start:',guess_p,"\n"
		if (1e-19 not in guess_p) or ('gl' in dist.__name__): #use SLSQP's constants/f_eqcon leads to singularity problem too often
			guess_p=optimize.fmin_slsqp(PFF, x0=map(log, guess_p), fprime=PPF, args=(dlist,), disp=0, full_output=1)#;print 'SLSQP:', exp(guess_p[0]), guess_p[1], "\n"; 
			guess_p=exp(guess_p[0])
		else:
			guess_p=exp(optimize.fmin(partial_fixed, x0=map(log, guess_p), 
                                  args=(dlist, dist, fixposls, True), disp=0))#;print 'Nelder-Mead:', guess_p
		#guess_p=optimize.fsolve(derv, guess_p, args=(dlist, lambda X, Y: partial_fixed(X, Y, dist, fixposls)))
	else:
		guess_p=list(itertools.chain(*zip(*resultdic[dist.__name__])))
	VCOV=hessmatrix(log(guess_p), dlist, lambda X, Y: -partial_fixed(X, Y, dist, fixposls, True))
	L=partial_fixed(guess_p, dlist, dist, fixposls)
	n=sum([item.shape[1] for item in dlist])*1.
	k=len(guess_p)*1.
	AIC=2*k-L
	AICc=(2.*k*n/(n-k-1))-L
	BIC=k*log(n)-L#;print L
	return guess_p, matrix(VCOV), L, AIC, AICc, BIC, int(k)

def head_first(head, col_name):
	result=[['<TH COLSPAN="%s"><BR><H3>%s</H3></TH>'%(len(col_name), head)]]
	result.append(['<TH>%s</TH>'%item for item in col_name])
	return result

def data_table(data_list):
	max_len=max(map(len, data_list))
	result=[['<TD>%s</TD>'%x[i] for x in data_list if (len(x)>1) or (i==0)]\
		for i in range(max_len)]
	if 1 in map(len, data_list):
		result[0]=['<TD>%s</TD>'%x[0] if len(x)>1 \
			   else '<TD ROWSPAN="%s">%s</TD>'%(max_len, x[0])\
			   for x in data_list]
	return result
	
def make_table(head, col_name, data_list):
	result=['<TABLE style="float: left;" BORDER="5" CELLPADDING="1" CELLSPACING="2">']
	for item in itertools.chain(*(head_first(head, col_name), data_table(data_list))):
		result.append('    <TR ALIGN="CENTER">')
		for jtem in item:
			result.append('        '+jtem)
		result.append('    </TR>')
	result.append('</TABLE>')
	return result

def convert_pf_result(in_list, Ns, Np, fixls):
	nv=(len(in_list)-len(fixls))/Ns
	result=[in_list[i*Ns:i*Ns+Ns] for i in range(nv)]
	result.extend([[item] for item in in_list[-len(fixls):]])
	idx=[item for item in range(Np) if item not in fixls]+fixls
	_return=range(Np)
	for item, jtem in zip(idx, result):
		_return[item]=jtem
	return _return

def result2html(PF_result, dist_name, Ns, Np, fixls, user_title='User Supplied Title'):
	head_dic={'gp2':('Gompertz',['&#945', '&#946']),
	          'gp3':('Gompertz-Makeham',['&#945', '&#946', 'C']),
	          'gl3':('Gompertz-Logistic',['&#945', '&#946', 'd']),
	          'gl4':('Gompertz-Logistic-Makeham', ['&#945', '&#946', 'C', 'd']),
	          'lg2':('Logistic',['&#956','s']),
	          'lg3':('Logistic-Makeham',['&#956','s', 'C']),
	          'wb2':('Weibull',['&#955','k']),
	          'wb3':('Weibull-Makeham',['&#955','k', 'C']),
	          'llg2':('log-Logistic',['&#945', '&#946',]),
	          'llg3':('log-Logistic-Makeham',['&#945', '&#946', 'C'])}
	st=['<h2>%s</h2>'%(user_title)] #'<!DOCTYPE html>', '<html>', '<body>',
	sb=itertools.chain(st,
	                   ['<p>log-Likelihood = %s;</br> AIC = %s; AICc = %s;</br> BIC = %s;</br> Degrees of Freedom = %s</p>'%PF_result[-5:]],
	                  ['<TABLE border="0">','<TD>'],
	                  make_table(head_dic[dist_name][0], head_dic[dist_name][1], convert_pf_result(PF_result[0], Ns, Np, fixls)),
	                  ['</TD>','<TD>'],
	                  make_table("&#963 of log(parameters)", head_dic[dist_name][1], convert_pf_result(np.sqrt(PF_result[1].diagonal()).tolist()[0], Ns, Np, fixls)),
	                  ['</TD>', '</TABLE>'])# ,'</body>','</html>'])
	return list(sb)
def all_partial_fixed_model(pdict, datals, dist, usr_title='User Supplied Title'):
	head_dic={'gp2':('Gompertz',['&#945', '&#946']),
	          'gp3':('Gompertz-Makeham',['&#945', '&#946', 'C']),
	          'gl3':('Gompertz-Logistic',['&#945', '&#946', 'd']),
	          'gl4':('Gompertz-Logistic-Makeham', ['&#945', '&#946', 'C', 'd']),
	          'lg2':('Logistic',['&#956','s']),
	          'lg3':('Logistic-Makeham',['&#956','s', 'C']),
	          'wb2':('Weibull',['&#955','k']),
	          'wb3':('Weibull-Makeham',['&#955','k', 'C']),
	          'llg2':('log-Logistic',['&#945', '&#946',]),
	          'llg3':('log-Logistic-Makeham',['&#945', '&#946', 'C'])}
	n_p=int(dist.__name__[-1])
	n_s=len(datals)
	a=list(itertools.product(*[[0,1]]*n_p))[::-1]
	a.sort(key=sum)
	b=[[i for i, val in enumerate(item) if val==1] for item in a]
	a=[['-' if val==0 else 'Fixed' for val in item] for item in a]
	html_list=[]
	L=[]
	AIC=[]
	AICc=[]
	BIC=[]
	DF=[]
	LK=[]
	for i,item in enumerate(b):
		r=fit_partial_fixed(pdict, datals, dist, item)
		html_list.append(result2html(r,dist.__name__,n_s,n_p,item, 
		                             'Model #%s<a id="Model_%s"></a><a href="#Summary" title="Back to Summary">&#x21a9</a>'%(i,i)))
		L.append(r[-5])
		AIC.append(r[-4])
		AICc.append(r[-3])
		BIC.append(r[-2])
		DF.append(r[-1])
		LK.append("<a href='#Model_%s'>Model #%s</a>"%(i,i))
		#print item
	s=make_table(head_dic[dist.__name__][0], head_dic[dist.__name__][1]+['L','AIC','AICc','BIC','DF','LK'], 
	             list(itertools.chain(zip(*a),[L,AIC,AICc,BIC,DF,LK])))
	s.insert(0, '<a id="Summary"></a><h1>%s</h1>'%usr_title)
	s[1]=s[1].replace('style="float: left;"', '')
	return itertools.chain(*itertools.chain([s],html_list))

class partial_fixed_HTML(object):
	def __init__(self, result_txt, data_txt, ndata=2, limit=0):
		with open(result_txt) as f:
			data=f.readlines()
		data=data[1:]
		data=[item.split('\t') for item in data]
		pdict={}
		for dist in ['lg2','gp2','wb2','lg3','gp3','wb3','llg3','llg4','gl3','gl4']:
			if limit>0:
				pdict[dist]=[[float(item[i]) for i in [8,10,12,14][:int(dist[-1])]] for item in data if dist in item][:limit]
			else:
				pdict[dist]=[[float(item[i]) for i in [8,10,12,14][:int(dist[-1])]] for item in data if dist in item]
		self.pdict=pdict
		if limit>0:
			self.dlist=input_data(data_txt, ndata)[0][:limit]
		else:
			self.dlist=input_data(data_txt, ndata)[0]
		self.have_html=('HTML' in dir())
	def summary(self, dist, title='User Supplied Title'):
		s=all_partial_fixed_model(self.pdict, self.dlist, dist, usr_title=title)
		s='\n'.join(s)
		return s

def read_fitresult(fit_fn):
	f1=open(fit_fn)
	head=f1.readline()[:-1].split('\t')
	didx=head.index('Dist.')
	pidx=[i for i, val in enumerate(head) if 'P_' in val and 'v' not in val]
	resultdic={}
	while True:
		eachline=f1.readline()
		if not eachline:
			break
		else:
			eachline=eachline[:-1].split('\t')
			if eachline[didx] not in resultdic:
				resultdic[eachline[didx]]=[]
			resultdic[eachline[didx]].append([float(eachline[i]) for i in pidx if i<len(eachline)])
	f1.close()
	return resultdic

#to try: do these:
#atest21=atest2[:, random.randint(0, 162, size=(100,))]
#atest22=atest2[:, random.randint(0, 162, size=(100,))]
#atest23=atest2[:, random.randint(0, 162, size=(100,))]
#rlist=[{},{},{}]
#fit_2pdist(atest21, gp2, rlist[0])
#fit_2pdist(atest22, gp2, rlist[1])
#fit_2pdist(atest23, gp2, rlist[2])
#rrdict={'gp2': [item['gp2'][0] for item in rlist]}
#fit_partial_fixed(rrdict, dlist, gp2, [0,1])
#optimize.fsolve(derv, pguess, args=(dguess, lambda X, Y: partial_fixed(X, Y, gp2, [0,])), full_output=1)
#optimize.fmin(partial_fixed, pguess, args=(dguess, gp2, [0,]))
#==========<End Partially fixing one or more parameters>==================================================

#==========<Plotting result and write summary statistics>=================================================
clh=lambda v, m: exp(log(m)+1.9599639845400543*sqrt(v)/m)
cll=lambda v, m: exp(log(m)-1.9599639845400543*sqrt(v)/m)
chisqlevels2=append([-0.001,],0.5*stats.chi2.ppf(hstack((linspace(0.1, 0.9, 9), array((0.95, 0.99)))),2))
chisqlevels3=append([-0.001,],0.5*stats.chi2.ppf(hstack((linspace(0.1, 0.9, 9), array((0.95, 0.99)))),3))

def cov2cor(mcov):
	'''Convert variance-covariance matrix to correlation matrix
	'''
	return mcov/sqrt(mcov.diagonal().T*mcov.diagonal())

def convert_str(x): #Retained from v0.6.
	'''A customized float to sting convertor, used in drawing fitting distribution plot.
	'''
	if math.isnan(x):
		result= 'Nan'
	else:
		#if 1.0e-14<abs(x)<0.01 or abs(x)>=1000.:
		if 0.<abs(x)<0.01 or abs(x)>=1000.:
			result= "%.2e" %abs(x)
			result=result.replace('e', r'\cdot 10^{')
			result=result.replace('-0', '-')
			result=result.replace('+0', '')
			result=result+'}'
		#elif abs(x)<1.0e-14:
		#	result='0'
		elif x==0.:
			result='0'
		else:
			result= "%3.3f" %abs(x)
		if abs(x)!=x:
			result='-'+result
	result=r'$'+result+'$'
	return result

class plotting_result(object):
	_colors={'gp2':'g',
			  'gp3':'g',
			  'gl3':'c',
			  'gl4':'c',
			  'lg2':'r',
			  'lg3':'r',
			  'wb2':'b',
			  'wb3':'b',
			  'llg2':'m',
			  'llg3':'m'}
	_dashes=[[5,2,2,2,2,2,5,0.1], [5,2,2,2,2,2,2,2,5,0.1]]
	_labels={'gp2':'Gompertz',
			  'gp3':'Gompertz-3P',
			  'gl3':'Gompertz-Logistic',
			  'gl4':'Gompertz-Logistic-4P',
			  'lg2':'Logistic',
			  'lg3':'Logistic-3P',
			  'wb2':'Weibull',
			  'wb3':'Weibull-3P',
			  'llg2':'log-Logistic',
			  'llg3':'log-Logistic-3P'}
	def __init__(self, resultdic):
		self.resultdic=resultdic
		self.maxl=max(resultdic['evttime'])
		self.xval=linspace(0,1.005,100)*max(resultdic['evttime'])
		self.xhzd=linspace(resultdic['evttime'][1]*0.95,1.005*max(resultdic['evttime']),100)
	def plot_cdf(self, ax, distLS, show_label=False):
		_maxl=self.maxl
		_xval=self.xval
		if not isinstance(distLS, list):
			distLS=[distLS]
		#first: plot emperical data
		ax.plot(self.resultdic['evttime'], 1-self.resultdic['empcdf'], 'k+', label='0bserved')
		#then: plot the fitted data
		for dist in distLS:
			l,=ax.plot(_xval,
					   1-dist.cdf(self.resultdic['fitresult'][dist.__name__][0], _xval),
					   color=self._colors[dist.__name__], alpha=0.7, label=self._labels[dist.__name__])
			if dist.__name__ in ['gp2', 'wb2', 'lg2', 'llg2', 'gl3']:
				l.set_dashes(self._dashes[0])
			else:
				l.set_dashes(self._dashes[1])
		ax.set_ylim((-0.01, 1.01))
		ax.set_xlim((-0.01*_maxl,1.01*_maxl))
		ax.set_title('Survivorship (1-Cumulative Probability)')
		if show_label:
			ax.legend(loc='lower left', frameon=False)
	def plot_pdf(self, ax, distLS, show_label=False):
		_maxl=self.maxl
		_xval=self.xval
		if not isinstance(distLS, list):
			distLS=[distLS]
		#first: plot emperical data
		ax.hist(self.resultdic['midpoints'], 
				bins=max([int(sqrt(len(self.resultdic['midpoints'])))+1,10]), normed=True, alpha=0.3)
		#then: plot the fitted data
		for dist in distLS:
			l,=ax.plot(_xval,
					   dist.pdf(self.resultdic['fitresult'][dist.__name__][0], _xval),
					   color=self._colors[dist.__name__], alpha=0.7, label=self._labels[dist.__name__])
			if dist.__name__ in ['gp2', 'wb2', 'lg2', 'llg2', 'gl3']:
				l.set_dashes(self._dashes[0])
			else:
				l.set_dashes(self._dashes[1])
		ax.set_ylim(ymin=-0.01*ax.get_ylim()[1])
		ax.set_xlim((-0.01*_maxl,1.01*_maxl))
		ax.set_title('Incidence of death (Probability Density)')
		if show_label:
			ax.legend(loc='upper left', frameon=False)
	def plot_hzd(self, ax, distLS, show_label=False):
		_maxl=self.maxl
		_xval=self.xhzd
		if not isinstance(distLS, list):
			distLS=[distLS]
		#first: plot emperical data
		lhzd=log(self.resultdic['emphzd'])
		ax.plot(self.resultdic['evttime'], lhzd, 'k+')
		hzdlim_temp=(floor(min(lhzd[isfinite(lhzd)])), ceil(max(lhzd[isfinite(lhzd)])))
		#then: plot the fitted data
		for dist in distLS:
			l,=ax.plot(_xval,
					   log(dist.hzd(self.resultdic['fitresult'][dist.__name__][0], _xval)),
					   color=self._colors[dist.__name__], alpha=0.7, label=self._labels[dist.__name__])
			if dist.__name__ in ['gp2', 'wb2', 'lg2', 'llg2', 'gl3']:
				l.set_dashes(self._dashes[0])
			else:
				l.set_dashes(self._dashes[1])
		ax.set_ylim(hzdlim_temp)
		ax.set_xlim((-0.01*_maxl,1.01*_maxl))
		ax.set_title(r'$log($'+'Mortality Rate'+r'$)$'+r' ($log($'+r'Hazard$)$)')
		if show_label:
			ax.legend(loc='upper left', frameon=False)
	def plot_numbers(self, ax, distLS, show_label=False):
		result_dic=self.resultdic
		if not isinstance(distLS, list):
			distLS=[distLS]
		linespace=1./3/(len(distLS))
		for i, dist in enumerate(distLS):
			ax.text(-0.1, linespace*(3+4*i), dist.__doc__,va='center',color=self._colors[dist.__name__])
			ax.text(0.6, linespace*(3+4*i), r'$log(\mathcal{L})=$'+convert_str(result_dic['fitresult'][dist.__name__][2]),va='center')
			ax.text(-0.08, linespace*(1+4*i), r'$h(x)=$'+dist.latexhzd,va='center')
			if dist in (gl4, gl3):
				xposlist=[0.26, 0.444, 0.637, 0.82]
			else:
				xposlist=[0.26, 0.54, 0.82]
			for j in range(int(dist.__name__[-1])):
				ax.text(xposlist[j], linespace*(1+4*i), dist.latexp[j]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][j])+
					 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][j,j]),va='center')
		ax.axis('off')
	def plot_confi(self, ax, dist, title='User Supplied Title'):
		result_dic=self.resultdic
		if 'gp' in dist.__name__:
			distLS=[gp2, gp3]
		elif 'wb' in dist.__name__:
			distLS=[wb2, wb3]
		elif 'gl' in dist.__name__:
			distLS=[gp2, gl3]
		elif 'llg' in dist.__name__:
			distLS=[llg2, llg3]
		elif 'lg' in dist.__name__:
			distLS=[lg2, lg3]
		else:
			raise NotImplementedError('Distribution not supported')
		def add_sub_axes(axis, rect): #@SO, unutbu
			def axis_to_fig(axis): #@SO, unutbu
				def transform(coord):
					return fig.transFigure.inverted().transform(axis.transAxes.transform(coord))
				return transform
			fig = axis.figure
			left, bottom, width, height = rect
			trans = axis_to_fig(axis)
			figleft, figbottom = trans((left, bottom))
			figwidth, figheight = trans([width,height]) - trans([0,0])
			return fig.add_axes([figleft, figbottom, figwidth, figheight])
		#Confidence contourf for the 2P model
		sax1=add_sub_axes(ax, (0.0 ,0.45,0.4,0.4))
		ctf=sax1.contourf(result_dic['fitresult'][distLS[0].__name__][3][0], 
				  result_dic['fitresult'][distLS[0].__name__][3][1], 
				  result_dic['fitresult'][distLS[0].__name__][3][2], levels=chisqlevels2)
		sax1.set_title(distLS[0].__name__)
		sax1.set_xlim(mean(sax1.get_xlim())+array([-1,1])*max(diff(sax1.get_xlim()), diff(sax1.get_ylim()))/2)
		sax1.set_ylim(mean(sax1.get_ylim())+array([-1,1])*max(diff(sax1.get_xlim()), diff(sax1.get_ylim()))/2)
		sax1.set_xlabel(r'$log($'+distLS[0].latexp[0]+r'$)$')
		sax1.set_ylabel(r'$log($'+distLS[0].latexp[1]+r'$)$')
		sax1.set_title(distLS[0].__doc__)
		sax1.set_aspect('equal')
		
		#Confidence Contourf for the 3p model
		saxL=[add_sub_axes(ax , item) for item in [(0.5 ,0.45,0.4,0.4), 
                                                   (0.0 ,0.0 ,0.4,0.4), 
                                                   (0.5 ,0.0 ,0.4,0.4)]]
		_3D_ctf_Xidx=[0,0,1]
		_3D_ctf_Yidx=[1,2,2]
		_3D_ctf_Vidx=[3,4,5]
		_3D_ctf_axis=[r'$log($'+item+r'$)$' for item in distLS[1].latexp]
		for i in range(3):
			saxL[i].contourf(result_dic['fitresult'][distLS[1].__name__][3][_3D_ctf_Xidx[i]], 
					 result_dic['fitresult'][distLS[1].__name__][3][_3D_ctf_Yidx[i]], 
					 result_dic['fitresult'][distLS[1].__name__][3][_3D_ctf_Vidx[i]], levels=chisqlevels3)
			saxL[i].set_title(distLS[1].__doc__+'\n ('+_3D_ctf_axis[_3D_ctf_Xidx[i]]+' v.s. '+_3D_ctf_axis[_3D_ctf_Yidx[i]]+') ')
			saxL[i].set_xlim(mean(saxL[i].get_xlim())+array([-1,1])*max(diff(saxL[i].get_xlim()), diff(saxL[i].get_ylim()))/2)
			saxL[i].set_ylim(mean(saxL[i].get_ylim())+array([-1,1])*max(diff(saxL[i].get_xlim()), diff(saxL[i].get_ylim()))/2)
			saxL[i].set_xlabel(_3D_ctf_axis[_3D_ctf_Xidx[i]])
			saxL[i].set_ylabel(_3D_ctf_axis[_3D_ctf_Yidx[i]])
			saxL[i].set_aspect('equal')
		saxc=add_sub_axes(ax, (0.93, 0.025 ,0.03 ,0.8))
		cb=plt.colorbar(ctf, cax=saxc, ticks=chisqlevels2)
		cb.ax.set_yticklabels(map(str, hstack((linspace(0.1, 0.9, 9), array((0.95, 0.99)))).tolist()))
		
		#The numbers
		saxt=add_sub_axes(ax, (0.0 ,0.9 ,1. ,0.1))
		saxt.text(0.0 ,0.9,u'Confidence regions for of parameter estimates for '+distLS[0].__doc__, fontsize=12, va='center')
		saxt.text(0.55,0.6,u'without constant', fontsize=12, va='center', ha='right')
		saxt.text(0.95,0.6,u'with constant', fontsize=12, va='center', ha='right')
		saxt.text(0.0 ,0.3,'Cumulative dist. function (CDF)', fontsize=12, va='center')
		saxt.text(0.55,0.3, distLS[0].latexcdf, fontsize=16, va='center', ha='right')
		saxt.text(0.95,0.3, distLS[1].latexcdf, fontsize=16, va='center', ha='right')
		saxt.text(0.0 ,0.0,'Hazard function', fontsize=12, va='center')
		saxt.text(0.55,0.0, distLS[0].latexhzd, fontsize=16, va='center', ha='right')
		saxt.text(0.95,0.0, distLS[1].latexhzd, fontsize=16, va='center', ha='right')
		saxt.axis('off')
		ax.axis('off')
		ax.set_title(title, position=(0.5, 1.02), size='xx-large')
	def fig_fit(self, title='User Supplied Title', figsize=(14,16), dpi=100):
		f=plt.figure(figsize=figsize, dpi=dpi)
		ax1=f.add_subplot(4,2,1)
		ax2=f.add_subplot(4,2,3)
		ax3=f.add_subplot(4,2,5)
		ax4=f.add_subplot(4,2,7)
		self.plot_numbers(ax1, [gp2, wb2, lg2, llg2, gl3], True)
		self.plot_cdf(ax2, [gp2, wb2, lg2, llg2, gl3], True)
		self.plot_pdf(ax3, [gp2, wb2, lg2, llg2, gl3] )
		self.plot_hzd(ax4, [gp2, wb2, lg2, llg2, gl3] )
		ax1=f.add_subplot(4,2,2)
		ax2=f.add_subplot(4,2,4)
		ax3=f.add_subplot(4,2,6)
		ax4=f.add_subplot(4,2,8)
		self.plot_numbers(ax1, [gp3, wb3, lg3, llg3, gl4], True)
		self.plot_cdf(ax2, [gp3, wb3, lg3, llg3, gl4], True)
		self.plot_pdf(ax3, [gp3, wb3, lg3, llg3, gl4] )
		self.plot_hzd(ax4, [gp3, wb3, lg3, llg3, gl4] )
		f.suptitle(title, size='xx-large')
		return f
	def fig_confi(self, dist, title='User Supplied Title', figsize=(14,16), dpi=100):
		f=plt.figure(figsize=figsize, dpi=dpi)
		ax=f.add_subplot(1,1,1)
		self.plot_confi(ax, dist, title=title)
		return f

#======The above version causes problem with multiprossessing, (numpy 0.7.1, scipy 1.1.1, matplotlib 0.12.0)
#======Until it can be fixed, mutltiprossessing has to use the following ugly code.

def plot_confi_region(result_dic, super_title='Title', folder=''):
	'''Plot multidimensional confidence limit.
	'''
	if isinstance(result_dic['fitresult']['lg2'][3][0], (float, int)):
		pass
	else:
		try:
			for item, item3 in zip([lg2, gp2, wb2, gp2, llg2], [lg3, gp3, wb3, gl3, llg3]):
				gs = gridspec.GridSpec(3, 3,
								   width_ratios=[5,5,0.25],
								   height_ratios=[1,5,5]
								   )
				fig=plt.figure(figsize=(10,11.75))
				plt.subplot(gs[1,0], aspect='equal')
				#X, Y=meshgrid(result_dic['fitresult'][item.__name__][3][0], result_dic['fitresult'][item.__name__][3][1])
				ctf=plt.contourf(result_dic['fitresult'][item.__name__][3][0], result_dic['fitresult'][item.__name__][3][1], result_dic['fitresult'][item.__name__][3][2], levels=chisqlevels2)
				plt.title(item.__name__)
				plt.xlim(mean(plt.xlim())+array([-1,1])*max(diff(plt.xlim()), diff(plt.ylim()))/2)
				plt.ylim(mean(plt.ylim())+array([-1,1])*max(diff(plt.xlim()), diff(plt.ylim()))/2)
				plt.xlabel(r'$log($'+item.latexp[0]+r'$)$')
				plt.ylabel(r'$log($'+item.latexp[1]+r'$)$')
				plt.title(item.__doc__)

				grid_to_plot=(gs[1,1],gs[2,0],gs[2,1])
				_3D_ctf_Xidx=[0,0,1]
				_3D_ctf_Yidx=[1,2,2]
				_3D_ctf_Vidx=[3,4,5]
				_3D_ctf_axis=[r'$log($'+iitem+r'$)$' for iitem in item3.latexp]
				for i in range(3):
					plt.subplot(grid_to_plot[i], aspect='equal')
					#X, Y=meshgrid(result_dic['fitresult'][item3.__name__][3][_3D_ctf_Xidx[i]], result_dic['fitresult'][item3.__name__][3][_3D_ctf_Yidx[i]])
					plt.contourf(result_dic['fitresult'][item3.__name__][3][_3D_ctf_Xidx[i]], result_dic['fitresult'][item3.__name__][3][_3D_ctf_Yidx[i]], result_dic['fitresult'][item3.__name__][3][_3D_ctf_Vidx[i]], levels=chisqlevels3)
					plt.title(item3.__doc__+'\n ('+_3D_ctf_axis[_3D_ctf_Xidx[i]]+' v.s. '+_3D_ctf_axis[_3D_ctf_Yidx[i]]+') ')
					plt.xlim(mean(plt.xlim())+array([-1,1])*max(diff(plt.xlim()), diff(plt.ylim()))/2)
					plt.ylim(mean(plt.ylim())+array([-1,1])*max(diff(plt.xlim()), diff(plt.ylim()))/2)
					plt.xlabel(_3D_ctf_axis[_3D_ctf_Xidx[i]])
					plt.ylabel(_3D_ctf_axis[_3D_ctf_Yidx[i]])
				ax=plt.subplot(gs[1:,2])
				cb=plt.colorbar(ctf, cax=ax, ticks=chisqlevels2)
				cb.ax.set_yticklabels(map(str, hstack((linspace(0.1, 0.9, 9), array((0.95, 0.99)))).tolist()))

				ax=plt.subplot(gs[0,:])
				plt.text(0,0.8,u'Confidence regions for of parameter estimates for '+item.__doc__, fontsize=12, va='center')
				plt.text(0.55,0.55,u'without constant', fontsize=12, va='center', ha='right')
				plt.text(0.95,0.55,u'with constant', fontsize=12, va='center', ha='right')
				plt.text(0,0.15,'Cumulative dist. function (CDF)', fontsize=12, va='center')
				plt.text(0.55, 0.15, item.latexcdf, fontsize=16, va='center', ha='right')
				plt.text(0.95, 0.15, item3.latexcdf, fontsize=16, va='center', ha='right')
				plt.text(0,-0.375,'Hazard function', fontsize=12, va='center')
				plt.text(0.55, -0.375, item.latexhzd, fontsize=16, va='center', ha='right')
				plt.text(0.95, -0.375, item3.latexhzd, fontsize=16, va='center', ha='right')
				plt.text(0,-0.65, ' ', fontsize=16)
				plt.title(super_title, fontsize=20)
				ax.axis('off')
				fig.tight_layout()
				plt.savefig(folder+super_title+'_'+item3.__name__[:-1]+'.pdf')
				plt.savefig(folder+super_title+'_'+item3.__name__[:-1]+'.png')
				plt.close()
		except KeyError:
			pass

def plot_fittes(result_dic, super_title='Title', folder=''):
	'''Plot fitted distributions, parameter estimates and variance, for Weibull, Logistic and Gompertz. 
	'''
	_colors=['r','g','b']
	_dashes=[[5,2,2,2,2,2,5,0.1], [5,2,2,2,2,2,2,2,5,0.1]]
	_title=['Survivorship (1-Cumulative Probability)', 'Incidence of death (Probability Density)', r'$log($'+'Mortality Rate'+r'$)$'+r' ($log($'+r'Hazard$)$)']
	_labels=['Logistic', 'Gompertz', 'Weibull']
	fig=plt.figure(figsize=(10,11.75))
	gs = gridspec.GridSpec(5,2,
			       width_ratios=[1,1],
			       height_ratios=[0.25,4.5,5,5,5])
	plt.subplot(gs[2,1])
	for i, dist in enumerate([lg3, gp3, wb3]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    1-dist.cdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7, label=_labels[i])
		l.set_dashes(_dashes[1])
	plt.plot(result_dic['evttime'], 1-result_dic['empcdf'], 'k+', label='Observed')
	lgd=plt.legend(loc='lower left', frameon=False)
	plt.ylim((-0.01,1.01))
	xlim_temp=plt.xlim()
	plt.title(_title[0])
	plt.subplot(gs[3,1])
	for i, dist in enumerate([lg3, gp3, wb3]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    dist.pdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[1])
	#plt.plot(result_dic['evttime'], result_dic['emppdf'], 'k+')
	plt.hist(result_dic['midpoints'], bins=max([int(sqrt(len(result_dic['midpoints'])))+1,10]), normed=True, alpha=0.3)
	plt.ylim(ymin=-0.01*plt.ylim()[1])
	plt.xlim(xlim_temp)
	plt.title(_title[1])
	plt.subplot(gs[4,1])
	for i, dist in enumerate([lg3, gp3, wb3]):
		l,=plt.plot(linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100),
			    log(dist.hzd(result_dic['fitresult'][dist.__name__][0], linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100))),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[1])
	lhzd=log(result_dic['emphzd'])
	plt.plot(result_dic['evttime'], lhzd, 'k+')
	hzdlim_temp=(floor(min(lhzd[isfinite(lhzd)])), ceil(max(lhzd[isfinite(lhzd)])))
	plt.ylim(hzdlim_temp)
	plt.xlim(xlim_temp)
	plt.title(_title[2])
	lgd.get_texts()[0].set_fontsize(10)
	plt.xlabel('Time')
	ax=plt.subplot(gs[1,1])
	for i, dist in enumerate([lg3, gp3, wb3]):
		plt.text(-0.1, 0.2+0.4*i, dist.__doc__,va='center')
		plt.text(0.6, 0.2+0.4*i, r'$log(\mathcal{L})=$'+convert_str(result_dic['fitresult'][dist.__name__][2]),va='center')
		plt.text(-0.08, 0.4*i, r'$h(x)=$'+dist.latexhzd,va='center')
		plt.text(0.26, 0.4*i, dist.latexp[0]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][0])+
			 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][0,0]),va='center')
		plt.text(0.54, 0.4*i, dist.latexp[1]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][1])+
			 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][1,1]),va='center')
		plt.text(0.82, 0.4*i, dist.latexp[2]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][2])+
			 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][2,2]),va='center')
	ax.axis('off')

	plt.subplot(gs[2,0])
	for i, dist in enumerate([lg2, gp2, wb2]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    1-dist.cdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7, label=_labels[i])
		l.set_dashes(_dashes[0])
	plt.plot(result_dic['evttime'], 1-result_dic['empcdf'], 'k+', label='Observed')
	lgd=plt.legend(loc='lower left', frameon=False)
	plt.ylim((-0.01,1.01))
	plt.xlim(xlim_temp)
	plt.title(_title[0])
	plt.subplot(gs[3,0])
	for i, dist in enumerate([lg2, gp2, wb2]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    dist.pdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[0])
	#plt.plot(result_dic['evttime'], result_dic['emppdf'], 'k+')
	plt.hist(result_dic['midpoints'], bins=max([int(sqrt(len(result_dic['midpoints'])))+1,10]), normed=True, alpha=0.3)
	plt.ylim(ymin=-0.01*plt.ylim()[1])
	plt.xlim(xlim_temp)
	plt.title(_title[1])
	plt.subplot(gs[4,0])
	for i, dist in enumerate([lg2, gp2, wb2]):
		l,=plt.plot(linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100),
			    log(dist.hzd(result_dic['fitresult'][dist.__name__][0], linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100))),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[0])
	plt.plot(result_dic['evttime'], lhzd, 'k+')
	plt.ylim(hzdlim_temp)
	plt.title(_title[2])
	lgd.get_texts()[0].set_fontsize(10)
	plt.xlim(xlim_temp)
	plt.xlabel('Time')
	ax=plt.subplot(gs[1,0])
	for i, dist in enumerate([lg2, gp2, wb2]):
		plt.text(0, 0.2+0.4*i, dist.__doc__,va='center')
		plt.text(0.6, 0.2+0.4*i, r'$log(\mathcal{L})=$'+convert_str(result_dic['fitresult'][dist.__name__][2]),va='center')
		plt.text(0.05, 0.4*i, r'$h(x)=$'+dist.latexhzd,va='center')
		plt.text(0.4, 0.4*i, dist.latexp[0]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][0])+
			 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][0,0]),va='center')
		plt.text(0.7, 0.4*i, dist.latexp[1]+u'='+convert_str(result_dic['fitresult'][dist.__name__][0][1])+
			 '\n'+r'$\sigma_{log}^{2}$'+'='+convert_str(result_dic['fitresult'][dist.__name__][1][1,1]),va='center')
	ax.axis('off')
	ax=plt.subplot(gs[0,:])
	plt.text(0.5, -0.5, 'Maximum Likelihood Parameter Estimates', fontsize=16, va='center', ha='center')
	plt.title(super_title, fontsize=20)
	ax.axis('off')
	fig.tight_layout()
	plt.savefig(folder+super_title+'_fit.pdf')
	plt.savefig(folder+super_title+'_fit.png')
	plt.close()

def plot_fittes2(result_dic, super_title='Title', folder=''):
	'''Plot fitted distributions, parameter estimates and variance, for Gompertz-logistic and log-logistic.
	'''
	_colors=['m','c','y']
	_dashes=[[5,2,2,2,2,2,5,0.1], [5,2,2,2,2,2,2,2,5,0.1]]
	_title=['Survivorship (1-Cumulative Probability)', 'Incidence of death (Probability Density)', r'$log($'+'Mortality Rate'+r'$)$'+r' ($log($'+r'Hazard$)$)']
	_labels=['Log-Logistic', 'Logistic-Gompertz']
	fig=plt.figure(figsize=(10,11.75))
	gs = gridspec.GridSpec(5,2,
			       width_ratios=[1,1],
			       height_ratios=[0.25,4.5,5,5,5])
	plt.subplot(gs[2,1])
	for i, dist in enumerate([llg3, gl4]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    1-dist.cdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7, label=_labels[i])
		l.set_dashes(_dashes[1])
	plt.plot(result_dic['evttime'], 1-result_dic['empcdf'], 'k+', label='Observed')
	lgd=plt.legend(loc='lower left', frameon=False)
	plt.ylim((-0.01,1.01))
	xlim_temp=plt.xlim()
	plt.title(_title[0])
	plt.subplot(gs[3,1])
	for i, dist in enumerate([llg3, gl4]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    dist.pdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[1])
	#plt.plot(result_dic['evttime'], result_dic['emppdf'], 'k+')
	plt.hist(result_dic['midpoints'], bins=max([int(sqrt(len(result_dic['midpoints'])))+1,10]), normed=True, alpha=0.3)
	plt.ylim(ymin=-0.01*plt.ylim()[1])
	plt.xlim(xlim_temp)
	plt.title(_title[1])
	plt.subplot(gs[4,1])
	for i, dist in enumerate([llg3, gl4]):
		l,=plt.plot(linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100),
			    log(dist.hzd(result_dic['fitresult'][dist.__name__][0], linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100))),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[1])
	lhzd=log(result_dic['emphzd'])
	plt.plot(result_dic['evttime'], lhzd, 'k+')
	hzdlim_temp=(floor(min(lhzd[isfinite(lhzd)])), ceil(max(lhzd[isfinite(lhzd)])))
	plt.ylim(hzdlim_temp)
	plt.xlim(xlim_temp)
	plt.title(_title[2])
	lgd.get_texts()[0].set_fontsize(10)
	plt.xlabel('Time')
	ax=plt.subplot(gs[1,1])
	for i, dist in enumerate([llg3, gl4]):
		plt.text(-0.1, 1.1-0.65*i, dist.__doc__,va='center')
		plt.text(0.6, 1.1-0.65*i, r'$log(\mathcal{L})=$'+convert_str(result_dic['fitresult'][dist.__name__][2]),va='center')
		plt.text(-0.08, 0.9-0.65*i, r'$h(x)=$'+dist.latexhzd, fontsize=15, va='center')
		if dist==gl4:
			xposlist=[-0.02, 0.26, 0.54, 0.82]
		else:
			xposlist=[0.26, 0.54, 0.82]
		for pos, xpos in enumerate(xposlist):
			plt.text(xpos, 0.65-0.65*i, dist.latexp[pos]+u'='+\
						 convert_str(result_dic['fitresult'][dist.__name__][0][pos])+\
						 '\n'+r'$\sigma_{log}^{2}$'+'='+\
						 convert_str(result_dic['fitresult'][dist.__name__][1][pos,pos]),va='center')
	ax.axis('off')

	plt.subplot(gs[2,0])
	for i, dist in enumerate([llg2, gl3]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    1-dist.cdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7, label=_labels[i])
		l.set_dashes(_dashes[0])
	plt.plot(result_dic['evttime'], 1-result_dic['empcdf'], 'k+', label='Observed')
	lgd=plt.legend(loc='lower left', frameon=False)
	plt.ylim((-0.01,1.01))
	plt.xlim(xlim_temp)
	plt.title(_title[0])
	plt.subplot(gs[3,0])
	for i, dist in enumerate([llg2, gl3]):
		l,=plt.plot(linspace(0,1.005,100)*max(result_dic['evttime']),
			    dist.pdf(result_dic['fitresult'][dist.__name__][0], linspace(0,1.005,100)*max(result_dic['evttime'])),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[0])
	#plt.plot(result_dic['evttime'], result_dic['emppdf'], 'k+')
	plt.hist(result_dic['midpoints'], bins=max([int(sqrt(len(result_dic['midpoints'])))+1,10]), normed=True, alpha=0.3)
	plt.ylim(ymin=-0.01*plt.ylim()[1])
	plt.xlim(xlim_temp)
	plt.title(_title[1])
	plt.subplot(gs[4,0])
	for i, dist in enumerate([llg2, gl3]):
		l,=plt.plot(linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100),
			    log(dist.hzd(result_dic['fitresult'][dist.__name__][0], linspace(result_dic['evttime'][1]*0.95,1.005*max(result_dic['evttime']),100))),
			    color=_colors[i], alpha=0.7)
		l.set_dashes(_dashes[0])
	plt.plot(result_dic['evttime'], lhzd, 'k+')
	plt.ylim(hzdlim_temp)
	plt.title(_title[2])
	lgd.get_texts()[0].set_fontsize(10)
	plt.xlim(xlim_temp)
	plt.xlabel('Time')
	ax=plt.subplot(gs[1,0])
	for i, dist in enumerate([llg2, gl3]):
		plt.text(0, 1.1-0.65*i, dist.__doc__,va='center')
		plt.text(0.6, 1.1-0.65*i, r'$log(\mathcal{L})=$'+convert_str(result_dic['fitresult'][dist.__name__][2]),va='center')
		plt.text(0.05, 0.9-0.65*i, r'$h(x)=$'+dist.latexhzd, fontsize=15, va='center')
		if dist==gl3:
			xposlist=[0.1, 0.4, 0.7]
		else:
			xposlist=[0.4, 0.7]
		for pos, xpos in enumerate(xposlist):
			plt.text(xpos, 0.65-0.65*i, dist.latexp[pos]+u'='+\
						 convert_str(result_dic['fitresult'][dist.__name__][0][pos])+\
						 '\n'+r'$\sigma_{log}^{2}$'+'='+\
						 convert_str(result_dic['fitresult'][dist.__name__][1][pos,pos]),va='center')
	ax.axis('off')
	ax=plt.subplot(gs[0,:])
	plt.text(0.5, -0.5, 'Maximum Likelihood Parameter Estimates', fontsize=16, va='center', ha='center')
	plt.title(super_title, fontsize=20)
	ax.axis('off')
	fig.tight_layout()
	plt.savefig(folder+super_title+'_fit2.pdf')
	plt.savefig(folder+super_title+'_fit2.png')
	plt.close()

def summary_stats(result_dic):
	'''Summarize log likelihood, parameter estimates and variance into a list, in with each item is the result for one distribution.
	'''
	_result=[]
	for dist in [lg2, gp2, wb2, lg3, gp3, wb3, llg2, llg3, gl3, gl4]:
		_estimate_var=list(chain(*zip(result_dic['fitresult'][dist.__name__][0], array(result_dic['fitresult'][dist.__name__][1]).diagonal())))
		_result.append(list(chain([dist.__name__,],result_dic['fitresult'][dist.__name__][2], \
					result_dic['fitresult'][dist.__name__][-3:], _estimate_var)))
	return _result

def container(input_arwg):
	'''Draw the confidence limit plot, the fitted distribution plot and write a summary text file containing the numerical result.
	Arguments: [array_in, treat_head, super_title, filename_out, folder]
	1, array_in: A time-to-event data array
	2, treat_head: A string, describing each item in array_in, such as 'Duck Male' or 'Chick Female'
	3, super_title: A string, describing the meaning of treat_head, such as 'Species Sex'
	4, filename_out: The filename.
	5, folder: folder path to store result
	'''
	array_in, treat_head, super_title, filename_Out, folder=input_arwg
	_result=adjust_Max_ls(array_in)
	with open((folder+filename_Out).replace('.txt','.cpl'), 'w')as f:
		cPickle.dump(_result, f)
	plot_confi_region(_result, treat_head, folder)
	plot_fittes(_result, treat_head, folder)
	plot_fittes2(_result, treat_head, folder)
	treat_head=treat_head.replace(' ','\t')
	if folder=='':
		folder=os.getcwd()+'/'
	if filename_Out not in os.listdir(os.path.split(folder+filename_Out)[0]):
		to_print=super_title+'\t'.join(['', 'Dist.', 'logL', 'AIC', 'AICc', 'BIC', \
							'P_1', 'vP_1', 'P_2', 'vP_2', 'P_3', 'vP_3', 'P_4', 'vP_4'])
		print>>open(folder+filename_Out,'a'), to_print
		to_print=super_title+'\t'.join(['','','Model','Constant_ignored','0.01','0.05','0.1','0.15','0.2',\
		'0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95',\
		'0.99','Mean'])
		print>>open(folder+'var_'+filename_Out,'a'), to_print
	for item in summary_stats(_result):
		print>>open(folder+filename_Out,'a'), treat_head+'\t'+'\t'.join(map(str, item)) 
	for item1, item2, item3 in zip(_result['var_log(qtm)'][::3], _result['var_log(qtm)'][1::3], _result['var_log(qtm)'][2::3]):
		print>>open(folder+'var_'+filename_Out,'a'), treat_head+'\tlog(time)\t'+item1+'\t'+'\t'.join(map(str, item2))
		print>>open(folder+'var_'+filename_Out,'a'), treat_head+'\tvariance\t'+item1+'\t'+'\t'.join(map(str, item3))
		
def z_test(filename, filenameout):
	'''Conduct 2 sample Z test for every possible comparisions.
	The input file: 'filename', must be a concatenated file of parameters and variances for all treatment combinations. Such as the return from combin_txt()
	returns: 
	1, A list of strings pairs, for each pair of treatments compared
	2, p-value matrix (upper-triangle.)
	'''
	dist_paras=['lg2_u', 'lg2_s', 'gp2_a', 'gp2_b', 'wb2_k', 'wb2_l', \
		    'lg3_u', 'lg3_s', 'lg3_c', 'gp3_a', 'gp3_b', 'gp3_c', 'wb3_k', 'wb3_l', 'wb3_c', \
		    'llg2_a', 'llg2_b', 'llg3_a', 'llg3_b', 'llg2_c', \
		    'gl3_a', 'gl3_b', 'gl3_d', 'gl4_a', 'gl4_b', 'gl4_c', 'gl4_d']
	no_paras=10 #6 if no llg2/3 and gl3/4
	f1=open(filename)
	lines=f1.readlines()
	f1.close()
	f1=open(filenameout, 'w')
	start_idxls=[[i for i, val in enumerate(item) if val=='\t'] for item in lines]
	start_idxidx=start_idxls[0].index(lines[0].index('\tP_1\t'))
	header=[item[:jtem[start_idxidx-5]] for item, jtem in zip(lines[:1], start_idxls[:1])]*2
	header=[item.replace('\t', str(i)+'\t')+str(i) for i, item in enumerate(header)]
	print>>f1, '\t'.join(header+list(chain(*[(item+'.Z', item+'.P') for item in dist_paras])))
	treatstr=[item[:jtem[start_idxidx-5]] for item, jtem in zip(lines[1:], start_idxls[1:])]
	est_SEM=[map(float, each.split('\t')) for each in [item[jtem[start_idxidx]+1:] for item, jtem in zip(lines, start_idxls)][1:]]
	compare=list(combinations(range(len(treatstr))[::no_paras],2))
	resultlist=[]
	resultlist_no=[]
	for item1, item2 in compare:
		_temp=[treatstr[item1]+'\t'+treatstr[item2],]
		_tempN=[]
		for j in range(no_paras):
			Z=[]
			for i in range(len(est_SEM[item1+j]))[::2]:
				try:
					Z.append(log(est_SEM[item1+j][i]/est_SEM[item2+j][i])/sqrt(est_SEM[item1+j][i+1]+est_SEM[item2+j][i+1]))
				except ZeroDivisionError:
					Z.append(float('nan'))
			Z=array(Z)
			P=2*norm.sf(abs(Z))
			_temp+=array([Z, P]).T.flatten().tolist()
			_tempN+=P.flatten().tolist()
		resultlist.append('\t'.join(map(str, _temp)))
		resultlist_no.append(_tempN)
	result_dic={}
	for i, item in enumerate(dist_paras):
		resultarray=zeros((len(treatstr)/no_paras, len(treatstr)/no_paras))+float('nan')
		for ix in range(len(treatstr)/no_paras):
			for iy in range(len(treatstr)/no_paras):
				if (ix*no_paras, iy*no_paras) in compare:
					resultarray[ix,iy]=resultlist_no[compare.index((ix*no_paras, iy*no_paras))][i]
				else:
					pass
		result_dic[item]=resultarray
	for item in resultlist:
		print>>f1, item
	f1.close()
	return treatstr[::no_paras], result_dic

def Ztest_heat_plt(zresult, lables, titlekw, name):
	'''Plot pairwise Z test heatmap. For one parameter.
	'''
	x=range(len(zresult))
	y=range(len(zresult))
	x=[item+0.5 for item in x]
	y=[item+0.5 for item in y]
	X, Y =meshgrid(x, y)
	log10result=fliplr(log(zresult))
	mmsk=[-10.0]*len(zresult)
	putmask(log10result, log10result<=mmsk, -10)
	plt.figure(figsize=(10,8))
	junk=plt.pcolormesh(ma.masked_invalid(log10result), lw=0.0, vmin=-10, vmax=0)
	cb=plt.colorbar(ticks=[-item for item in range(11)])
	cb.ax.set_yticklabels([r'$1.0$',r'$0.1$',r'$0.01$',r'$0.001$',r'$10^{-4}$',r'$10^{-5}$',r'$10^{-6}$',r'$10^{-7}$',r'$10^{-8}$',r'$10^{-9}$', r'$10^{-10}$'])
	junk=plt.xlim((-0.5, len(zresult)+0.5))
	junk=plt.ylim((-0.5, len(zresult)+0.5))
	junk=plt.axes().set_aspect(1)
	tpos=[item-0.5 for item in range(1,len(zresult)+1)]
	junk=plt.yticks(tpos)
	junk=plt.xticks(tpos)
	junk=plt.axes().set_yticklabels(lables, ha='right', va='center', fontsize=7)
	junk=plt.axes().set_xticklabels(lables[::-1], ha='center', va='top', rotation=90, fontsize=8)
	#plt.axes().xaxis.set_ticks_position("none")
	#plt.axes().yaxis.set_ticks_position("none")
	junk=plt.title('Pair-wise Z-test for '+titlekw+' ($p$ value in $log_{10}$ scale)')
	plt.savefig(name+'.png')
	plt.savefig(name+'.pdf')
	plt.close()

def mZtest_heat_plt(indata):
	'''Wrapper for mutliprocessing.
	'''
	Ztest_heat_plt(*indata)

def Z_test_plots(treatstr, p_values, namepdf):
	'''Make a multipage pdf of Z test heatmap, each page is for one parameter.  
	'''
	titles=[r'logistic (2p) $\mu$', r'logistic (2p) $s$',\
		r'Gompertz (2p) $\alpha$', r'Gompertz (2p) $\beta$',\
		r'Weibull (2p) $k$', r'Weibull (2p) $\lambda$',\
		r'logistic (3p) $\mu$', r'logistic (3p) $s$', r'logistic (3p) $c$',\
		r'Gompertz (3p) $\alpha$', r'Gompertz (3p) $\beta$', r'Gompertz (3p) $c$', \
		r'Weibull (3p) $k$', r'Weibull (3p) $\lambda$', r'Weibull (3p) $c$']
	treatstr=[item.replace('\t', ' ') for item in treatstr]
	pv_keys=['lg2_u', 'lg2_s', 'gp2_a', 'gp2_b', 'wb2_k', 'wb2_l', 'lg3_u', 'lg3_s', 'lg3_c', 'gp3_a', 'gp3_b', 'gp3_c', 'wb3_k', 'wb3_l', 'wb3_c']
	for i, val1, val2 in zip(range(len(pv_keys)), pv_keys, titles):
		Ztest_heat_plt(p_values[val1], treatstr, val2, 'Z_test_'+val1)
#       def temp_fuc(indata):
#	       Ztest_heat_plt(p_values[indata[0]], treatstr, indata[1], 'Z_test_'+indata[0])
#       rezt=ppool.map_async(temp_fuc, zip(pv_keys, titles))
#       rezt.get()
	combin_pdf(['Z_test_'+item+'.pdf' for item in pv_keys], namepdf)

def combin_pdf(pdflist, out_file, allpage=False):
	'''Combine a list of pdf files into one big pdf mutlipage file.
	'''
	output=PdfFileWriter()
	open_file_list=[open(item, 'rb') for item in pdflist]
	for item in open_file_list:
		pdfinput=PdfFileReader(item)
		if allpage:
			for i in range(pdfinput.getNumPages()):
				output.addPage(pdfinput.getPage(i))
		else:
			output.addPage(pdfinput.getPage(0))
	f1=open(out_file, 'wb')
	output.write(f1)
	f1.close()
	for item in open_file_list:
		item.close()

def mcombin_pdf(indata):
	'''Warper for multiprocessing.
	'''
	combin_pdf(*indata)

def combin_txt(txtlist, out_file):
	'''Combine may text file, each summarizing the parameter estimates for one treatment combination, into a big one.
	Attention, the output can be used as input for Z test. 
	'''
	out_str=[]
	for item in txtlist:
		f1=open(item, 'r')
		out_str.append(f1.readlines())
		f1.close()
	out_str=list(chain(*out_str))
	out_str=[out_str[0],]+[item for item in out_str[1:] if ('vP_1' not in item and 'Constant_ignored' not in item)]
	f1=open(out_file, 'w')
	for item in out_str:
		print>>f1, item[:-1]
	f1.close()

def AICBIC_whole(datals, fitresult_txt, outfile=None):
	'''Generate the model summary file, calculate the AIC, AICc and BIC for the models fitted to the whole data.
	'''
	n=sum([item.shape[1] for item in datals])
	f1=open(fitresult_txt)
	lines=f1.readlines()
	f1.close()
	start_idxls=[[i for i, val in enumerate(item) if val=='\t'] for item in lines]
	start_idx=start_idxls[0].index(lines[0].index('\tDist.\t'))
	dist_name=[]
	Llist=[]
	for item, jtem in zip(start_idxls[1:], lines[1:]):
		dist_name.append(jtem[item[start_idx]+1:item[start_idx+1]])
		Llist.append(float(jtem[item[start_idx+1]+1:item[start_idx+2]]))
	dist_list=list(set(dist_name))
	AIC=[]
	AICc=[]
	BIC=[]
	for item in dist_list:
		Llist_temp=[jtem for jtem, ktem in zip(Llist, dist_name) if ktem==item]
		L=2*sum(Llist_temp)
		k=len(Llist_temp)*float(item[-1])
		AIC.append(2*k-L)
		AICc.append((2.*k*n/(n-k-1))-L)
		BIC.append(k*log(n)-L)
	Model_rank={}
	for item, jtem in zip([AIC, AICc, BIC], ['By_AIC', 'By_AICc', 'By_BIC']):
		Model_rank[jtem]=[x for (y,x) in sorted(zip(item, dist_list))]
	if outfile:
		f1=open(outfile, 'w')
		print>>f1, 'Models:'
		print>>f1, 'wb2:  Weibull  distribution;     wb3:  Weibull  distribution with Makeham constant'
		print>>f1, 'lg2:  Logistic distribution;     lg3:  Logistic distribution with Makeham constant'
		print>>f1, 'llg2: Log-logistic distr.  ;     llg3: Log-logistic  distri. with Makeham constant'
		print>>f1, 'gp2:  Gompertz distribution;     gp3:  Gompertz-Makeham distribution'
		print>>f1, 'gl3:  Gompertz-Logistic distr.;  gl4:  Gompertz-Logistic-Makeham distribution'
		print>>f1, '=================================================================================='
		print>>f1, ''
		print>>f1, 'Model comparisons:'
		print>>f1, 'Models:\t'+'\t'.join(dist_list)
		print>>f1, 'AIC:\t'+'\t'.join(map(str, AIC))
		print>>f1, 'AICc:\t'+'\t'.join(map(str, AICc))
		print>>f1, 'BIC:\t'+'\t'.join(map(str, BIC))
		print>>f1, 'Total obs.:\t'+str(n)
		print>>f1, '=================================================================================='
		print>>f1, ''
		print>>f1, 'Model rank: (Left to right: small AIC, AICc or BIC value to large value, the smaller the better)'
		for item in ['By_AIC', 'By_AICc', 'By_BIC']:
			print>>f1, item[3:]+':\t'+'\t'.join(Model_rank[item])
		f1.close()
	else:
		return Model_rank, dist_list, AIC, AICc, BIC
#==========<End Plotting result and write summary statistics>=============================================
	
def test_speed():
	import time
	try:
		f1=open('speed_test.txt')
	except IOError:
		timelist=[]
		for i, item in enumerate([atest1, atest2, atest3]):
			timelist.append(time.time())
			container((item, 'speedtest_'+str(i), 'Speedtest', 'speedtest_'+str(i), ''))
		timelist.append(time.time())
		timelist=array(timelist)
		f1=open('speed_test.txt', 'w')
		result= float(mean(diff(timelist)))
		print>>f1,result
		f1.close()
	return result

if __name__ == '__main__':
	freeze_support()
	#print k
	if cpu_count()==1:
		cpuC=1
	else:
		cpuC=cpu_count() #2 #to limit it to 2
	pool=Pool(cpuC)
	for job in ['inital setups']:
		_versionName='ResultCv8_'
		_file_name=r'testdata'
		#define default action and dataset
		if len(sys.argv)==1:
			try:
				arvg=open('MLE.ini').readlines()[0][:-1].split(' ')
			except IOError:
				arvg=[_file_name, 2]
		else:
			arvg=sys.argv[1:]
		#ask for the scaling factor, can be omitted with calculation core v.073
		if len(arvg)==2:
			arvg.append(1.0)
			arvg[1]=int(arvg[1])
		else:
			arvg[2]=float(arvg[2])
			arvg[1]=int(arvg[1])
		#modifies the path to data input file:
		if os.path.split(arvg[0])[0]=='':
			arvg[0]=os.path.join(os.getcwd(),arvg[0])
		#make a new folder to store result, not used after v.06
		if os.getcwd()[-4:]=='/bin':
			storage_path=os.getcwd()[:-4]+'/data'
		else:
			storage_path=os.getcwd()+'/data'
		new_folder=(os.path.split(arvg[0]))[0]
		#read dataset
		_datals, _treatls, _head=input_data(arvg[0], arvg[1])
		_lendata=len(_datals)
		to_do_list=zip(_datals, _treatls, [_head,]*_lendata, [item+'.txt' for item in _treatls], [os.path.join(new_folder, _versionName), ]*_lendata)
		#make a sorted list of treatment strings
		txt_file_ls=[item.split(' ') for item in _treatls]
		for i in range(len(txt_file_ls[0]))[::-1]:
			txt_file_ls.sort(key=lambda x: x[i])
		txt_file_ls=[' '.join(item) for item in txt_file_ls]
		txt_file_ls=[_versionName+item+'.txt' for item in txt_file_ls]
		#make a pdf list of treatment strings
		pdf_file_ls=[]
		for item in txt_file_ls:
			for jtem in ['_fit.pdf', '_fit2.pdf', '_wb.pdf', '_gp.pdf', '_lg.pdf', '_llg.pdf', '_gl.pdf']:
				pdf_file_ls.append(item.replace('.txt', jtem))
	t1=time.time()
	#print k, 'for debug only'
	#to_do_list=to_do_list[:4]
	#_imap_obj=pool.imap(container, to_do_list)
	for idx in range(len(to_do_list)):
		#_imap_obj.next()
		container(to_do_list[idx])
		this_time=time.time()-t1
		if idx==0:
			time_str=os.path.join(new_folder, str(_lendata/cpuC*this_time))+'.temptime'
			f1=open(time_str,'w')
			f1.close()
		print this_time, "Time to finish:", this_time*(_lendata-idx-1)/(idx+1)/cpuC,"Just finished: '"+_treatls[idx]+"'"
	t2=time.time()
	os.chdir(new_folder)
	combin_txt(txt_file_ls, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.txt')
	combin_txt([item.replace(_versionName, _versionName+'var_') for item in txt_file_ls], \
		_versionName+'var_'+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.txt')
	q1, q2=z_test(_versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.txt', 'Ztest_'+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.txt')

	Ztitles=[r'logistic (2p) $\mu$', r'logistic (2p) $s$',\
                r'Gompertz (2p) $\alpha$', r'Gompertz (2p) $\beta$',\
		r'Weibull (2p) $k$', r'Weibull (2p) $\lambda$',\
		r'logistic (3p) $\mu$', r'logistic (3p) $s$', r'logistic (3p) $c$',\
		r'Gompertz (3p) $\alpha$', r'Gompertz (3p) $\beta$', r'Gompertz (3p) $c$', \
		r'Weibull (3p) $k$', r'Weibull (3p) $\lambda$', r'Weibull (3p) $c$',\
                r'Log-logistic (2p) $\alpha$', r'Log-logistic (2p) $\beta$',\
                r'Log-logistic (3p) $\alpha$', r'Log-logistic (3p) $\beta$', r'Log-logistic (3p) $c$',\
                r'Gompertz-Logistic (3p) $\alpha$', r'Gompertz-Logistic (3p) $\beta$', r'Gompertz-Logistic (3p) $d$', \
                r'Gompertz-Logistic (4p) $\alpha$', r'Gompertz-Logistic (4p) $\beta$', r'Gompertz-Logistic (4p) $c$', r'Gompertz-Logistic (4p) $d$']
	Ztreatstr=[[item.replace('\t', ' ') for item in q1],]*len(Ztitles)
	Zpv_keys=['lg2_u', 'lg2_s', 'gp2_a', 'gp2_b', 'wb2_k', 'wb2_l', \
                    'lg3_u', 'lg3_s', 'lg3_c', 'gp3_a', 'gp3_b', 'gp3_c', 'wb3_k', 'wb3_l', 'wb3_c', \
                    'llg2_a', 'llg2_b', 'llg3_a', 'llg3_b', 'llg2_c', \
                    'gl3_a', 'gl3_b', 'gl3_d', 'gl4_a', 'gl4_b', 'gl4_c', 'gl4_d']
	rezt1=pool.map_async(mZtest_heat_plt, zip([q2[item] for item in Zpv_keys], Ztreatstr, Ztitles, [os.path.join(new_folder, 'Z_test_'+item) for item in Zpv_keys]))
	rezt1.get()

	AICBIC_whole(_datals, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.txt', \
                     _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'_Model_Summary.txt')

#Combine PDF
#Attn: 1, if there are a lot of pdf files, first assemble three smaller files helps to avoid 'too many opened file' error.
	if len(pdf_file_ls)>3:
		spf=len(pdf_file_ls)/3
		rezt2=pool.map_async(mcombin_pdf, [\
			([os.path.join(new_folder, item) for item in pdf_file_ls[:spf]], \
			os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'0.pdf')),\
			([os.path.join(new_folder, item) for item in pdf_file_ls[spf:spf*2]], \
			os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'1.pdf')),\
			([os.path.join(new_folder, item) for item in pdf_file_ls[spf*2:]], \
			os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'2.pdf')),\
			([os.path.join(new_folder, 'Z_test_'+item+'.pdf') for item in Zpv_keys], \
			os.path.join(new_folder, 'Ztest_'+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.pdf'))])
		rezt2.get()
		mcombin_pdf(([os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'0.pdf'),\
				os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'1.pdf'),\
				os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'2.pdf')],
				os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.pdf'), True))
	else:
		rezt2=pool.map_async(mcombin_pdf, [\
			([os.path.join(new_folder, item) for item in pdf_file_ls], \
			os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.pdf')),\
			([os.path.join(new_folder, 'Z_test_'+item+'.pdf') for item in Zpv_keys], \
			os.path.join(new_folder, 'Ztest_'+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'.pdf'))])
		rezt2.get()
	
	import shutil
	os.mkdir(os.path.join(new_folder, 'MLEfit'))
	os.mkdir(os.path.join(new_folder, 'Variance'))
	for item in txt_file_ls:
		shutil.move(item, os.path.join(new_folder, 'MLEfit'))
		shutil.move(item.replace(_versionName, _versionName+'var_') , os.path.join(new_folder, 'Variance'))
		for jtem in ['_fit.pdf', '_fit2.pdf', '_wb.pdf', '_gp.pdf', '_lg.pdf', '_llg.pdf', '_gl.pdf', \
                             '_fit.png', '_fit2.png', '_wb.png', '_gp.png', '_lg.png', '_llg.png', '_gl.png', ]:
			shutil.move(item.replace('.txt', jtem), os.path.join(new_folder,'MLEfit'))
	os.mkdir(os.path.join(new_folder, 'Ztest'))
	for item in os.listdir(os.getcwd()):
		if 'Z_test_' in item:
			shutil.move(item, os.path.join(new_folder, 'Ztest'))
		if '.temptime' in item:
			os.remove(item)
	if len(pdf_file_ls)>3:
		for item in [os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'0.pdf'),\
				os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'1.pdf'),\
				os.path.join(new_folder, _versionName+list(os.path.split(arvg[0]))[-1].replace('.txt', '')+'2.pdf')]:
			os.remove(item)
	t3=time.time()
	print 'Total time:', t3-t1, 'All Done.'
#['/Users/user/Documents/python/MLE/dist/DeDayMLEcore.app/Contents/MacOS/DeDayMLEcore', '/Users/user/Documents/python/MLE/testdata', '2']
#['python', '/Users/user/Documents/MLE-current/DeDayMLEcore.py', '/Users/user/Documents/DeDay/4ebp_rapa/4ebp_rapa.txt', '2']

