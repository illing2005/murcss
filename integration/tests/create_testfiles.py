'''
Created on 15.05.2013

@author: Sebastian Illing
@contact: sebastian.illing@met.fu-berlin.de

Copyright (C) 2014  Sebastian Illing
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''

import os
from cdo import *
cdo = Cdo()

import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

folder = '/pf/b/b324057/goddard-metrics/test_files/test_data/'
prename = 'start_files.nc'
output = '/tmp/msss_test/test_data/'
if not os.path.isdir(output):
    os.makedirs(output)


def modifyStd(ts,std_goal=1):
    std = np.round(np.std(ts),2)
    if std == std_goal: return ts
    if std > std_goal:
        return modifyStd(0.9995*ts,std_goal)
    if std < std_goal:
        return modifyStd(1.0005*ts,std_goal)

def getTimeSeries(std1_goal, std2_goal, corr_goal, length=50):
    r=0
    while r!=corr_goal:
        val=np.random.multivariate_normal((0,0),[[std1_goal**2,corr_goal*(std1_goal*std2_goal)],[corr_goal*(std1_goal*std2_goal),std2_goal**2]],length)
        t = np.linspace(0, 1, length)
        T = 2*(t-np.mean(t))+val[:,0]
        T2 = 2*(t-np.mean(t))+val[:,1]
        r= np.round(np.corrcoef(T,T2),2)[0,1]
        #print r
        #print np.corrcoef(val[:,0],val[:,1])[0,1]
    T=modifyStd(T,std1_goal)
    T=T-np.mean(T)
    T2=modifyStd(T2,std2_goal)
    T2=T2-np.mean(T2)
    T=T+273
    T2=T2+273
    return (t,T,T2)

def createDecadal(startyear, value, flag):    
    tmp = cdo.setyear(startyear, input=cdo.setrtoc('-1000,1000,'+str(value), input=folder+prename), 
                      output=output+'test_files_'+flag+'_'+str(startyear))
    
    return tmp

def createMsssInput(startyear,std1_goal, std2_goal, corr_goal, length=50):
    
    (t,hind,obs) = getTimeSeries(std1_goal,std2_goal,corr_goal,length)    
    hindcasts = dict()
    observations = dict()
    for i,val in enumerate(hind):
        hindcasts[startyear+i] = [createDecadal(startyear+i, val, 'decadal')]
    for i,val in enumerate(obs):
        observations[startyear+i] = createDecadal(startyear+i, val, 'observation')
        
    return (hindcasts,observations)

def createCrpssInput(startyear,hind, obs, ensemble_size, length):
    
    #(t,hind,obs) = getCrpssTimeseries(mse_goal,ensspread_goal,ensemble_size,length)
    hindcasts = dict()
    observations = dict()
    
    for i in range(length):
        #print hind[i,:]
        ensemble_members = list()
        for j in range(ensemble_size):
            ensemble_members.append(createDecadal(startyear+i, hind[i,j], 'ens'+str(j+1)+'_decadal'))
        hindcasts[startyear+i] = ensemble_members
    for i,val in enumerate(obs):
        observations[startyear+i] = createDecadal(startyear+i, val, 'observation')
    print len(observations)
    return (hindcasts, observations)    
    

def getCrpssTimeseries(MSE_goal, Ensspread_goal = 2, ensemble_size=3, length=50):
        
    #val=np.random.multivariate_normal((0,0),[[std1_goal**2,corr_goal*(std1_goal*std2_goal)],[corr_goal*(std1_goal*std2_goal),std2_goal**2]],length)
    t = np.linspace(0,1,length)
    #Create Ensemble Mean
    val = np.random.normal(0,1,length)
    T = val + (t-np.mean(t))
    T = T - np.mean(T)
    factor = length/(length-2.) 
    MSE=0
    #Create Observations with fixed MSE
    while MSE != MSE_goal:
        Obs = T + np.random.normal(0,MSE_goal,length)
        Obs = Obs - np.mean(Obs)
        MSE = np.round((factor*np.mean((T-Obs)**2))**0.5,2)
    
    #create ensemble member
    MSE_ens=0
    Ensspread=0
    bias = 0
    bias_goal=1
    T_re=0
    T_re_tmp = 0
    j=0
    while MSE_ens != MSE_goal or Ensspread != Ensspread_goal:
        j+=1
        del T_re, T_re_tmp
        for i in range(ensemble_size):
            T_new = T + np.random.normal(0,Ensspread_goal,length)
            try:
                T_re = np.concatenate((T_re,T_new.reshape((T.size,1))),axis=1)
            except UnboundLocalError:
                T_re = T_new.reshape((T_new.size,1))
        T_ensmean = np.mean(T_re,axis=1)
        T_re_tmp = T_re
        
        #print Ensspread,MSE_ens
        #T_re = T_re - T_ensmean.reshape((T_ensmean.size,1))
        r = np.corrcoef(T_ensmean,Obs)[0,1]
        std_H = np.std(T_ensmean)
        std_O = np.std(Obs)
        bias = r * std_O/std_H
        T_ensmean = T_ensmean - np.mean(T_ensmean)
        #bias correction?
        T_re_tmp = T_re_tmp - np.mean(T_ensmean)
        T_re_tmp = T_re_tmp - T_ensmean.reshape((T_ensmean.size,1)) * (bias-1)
        T_ensmean2 = np.mean(T_re_tmp,axis=1)
        MSE_ens = np.round((factor*np.mean((T_ensmean2-Obs)**2))**0.5,2)
        Ensspread = np.round((np.mean(np.var(T_re_tmp,axis=1))**0.5),2)
              
    r = np.corrcoef(T_ensmean,Obs)[0,1]
    std_H = np.std(T_ensmean)
    std_O = np.std(Obs)
    
    bias = r * std_O/std_H
    
    return (t,T_re, Obs)

def crps(var,H,O):
    x = (H-O)/var
    crps = -var * (1./np.sqrt(np.pi) - 2. * stats.norm.pdf(x) - x * (2. * stats.norm.cdf(x) - 1.))
    return crps

def crpss(H,O):
    H_mean = np.mean(H,axis=1)
    ens_spread = (np.mean(np.var(H,axis=1)))**0.5
    Mse = (len(O)/(len(O)-2.)*np.mean((H_mean-O)**2))**0.5
    crps_Ens = list()
    crps_Mse = list()
    for i,val in enumerate(H_mean):
        #print i,val
        crps_Ens.append(crps(ens_spread,val,O[i]))
        crps_Mse.append(crps(Mse,val,O[i]))
    
    crpss = 1 - np.mean(crps_Ens)/np.mean(crps_Mse)
    
    return crpss

      
if __name__ == '__main__':
#    std1_goal = 1
#    std2_goal = 2
#    corr_goal = 0.5
#    startyear = 1960
#    fig = plt.figure(figsize=(12,8))
#    (t,T,T2) = getTimeSeries(std1_goal,std2_goal,corr_goal,50)
#    
#    ax = fig.add_subplot(111)
#    ax.plot(t,2*(t-np.mean(t))+273)
#    ax.plot(t,T)
#    ax.plot(t,T2)
#    
#    std1 = np.std(T)
#    std2 = np.std(T2)
#    print np.mean(T)
#    print np.mean(T2)
#    print std1,std2
#    print np.corrcoef(T,T2)[0,1]
    
    #create 10 decadal experiments (starting every 5 years
    #print createMurcssInput(startyear, T, T2)
    #createCrpssInput(1960,1,2,10,50)
    (t,T,Obs) = getCrpssTimeseries(3.14,2,10,50)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.plot(t,(t-np.mean(t)))
    ax.plot(t,T)
    ax.scatter(t,Obs)
    
    crpss(T,Obs)
    #print T.shape
    
    mse = 1
    ens_spread = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5]
    ens_sizes = [100,1000]
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    for ens_size in ens_sizes:
        print 'Ensemble Size',ens_size
        crpss_vs = list()
        for val in ens_spread:
            print val
            (t,H,O) = getCrpssTimeseries(mse, val, ens_size, 500)
            crpss_vs.append(crpss(H,O))
        ax.plot(ens_spread,crpss_vs)
    #print T[:,0]
    plt.show() 
    #createDecadal(1960, 280)
