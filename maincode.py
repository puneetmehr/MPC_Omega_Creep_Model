# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:54:57 2021

@author: pmehra
"""
import datetime as dt
import subroutine
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
#from IPython.display import display, Latex
apidata="api-data.xlsx"

omegadata = pd.read_excel(apidata, sheet_name=0, index_col=0, header=1)
omegadata = omegadata.drop(axis=1,columns='Material')
lmpdata = pd.read_excel(apidata, sheet_name=1, index_col=0, header=1)
lmpdata = lmpdata.drop(axis=1,columns='Material')
wrcdata = pd.read_excel(apidata, sheet_name=2, index_col=0, header=1)
wrcdata = wrcdata.drop(axis=1,columns='Material')
print("Please note that the subroutine is only available in ksi since the MPC Omega data in API 579 is only in ksi.")
project=input("Enter project number: ")
mats=int(input("Enter no. of creep materials defined in Abaqus: "))
mat=[]
matname=[]
peak=[]
temp=[]
#loop to input material details
i=0
for i in range(mats):
    while True:
        m=input("Pointer for material "+str(i+1)+": ")
        try:
            m=int(m)
        except:
            print("\nPlease enter an integer between 1 to 24")
            continue
        if int(m) <= 0 or int(m) > 24:
            print ("\nPlease enter the pointer from 1 to 24")
            continue
        break
    mat.append(m-1)    #to account for data index starting at 0
    mn=input("Enter material name to be used (like LCS_U)\n\nPlease use all caps. This should exactly match Abaqus material name: ")
    matname.append(mn)
    pk=input("Enter peak von Mises stress (in ksi) from 'static' step for "+mn+": ")
    peak.append(pk)
    tp=input("Enter temperature in \xb0F for "+mn+": ")
    temp.append(tp)
#########################################################################################################################
odata=omegadata.iloc[mat].to_numpy()
ldata=lmpdata.iloc[mat].to_numpy()
wdata=wrcdata.iloc[mat].to_numpy()
sindata=[]  
#setup materials for subroutine
#default omega model values, do not change
#############################
alpha=2.0                   #
beta_omega=0.33             #
delo=np.array([0.0,0.0])    #
#############################    
print("\nBy default, this code uses Δ_Ω_SR and Δ_Ω_CD = 0 ")
#display(Latex(f'$\Delta_\Omega^{{sr}}=0$'))
#display(Latex(f'$\Delta_\Omega^{{cd}}=0$'))
while True:
    delt=input("Do you want to change them? Enter 1 (yes) or 2 (No): ")
    try:
        delt=int(delt)
    except:
        print("Please enter 1 or 2 only.")
        continue
    if int(delt)<1 or int(delt)>2:
        print("Please enter 1 or 2 only.")
        continue
    break
i=0
if delt==1:
    delo=[]
    for i in range(mats):
        sr=float(input("Enter value for Strain Rate parameter for "+matname[i]+": "))
        cd=float(input("Enter value for Creep Ductility parameter for "+matname[i]+": "))
        del1=np.array([sr,cd])
        delo=np.append(delo,del1)
else:
    del1=delo
    for i in range(mats):
        delo=np.append(delo,del1)        
i=0
j=0
count=1
fig=plt.figure(figsize=(20,10))
for i in range(mats):
    del_sr=delo[j]
    j+=1
    del_cd=delo[j]
    j+=1
    a0=odata[i,0]
    a1=odata[i,1]
    a2=odata[i,2]
    a3=odata[i,3]
    a4=odata[i,4]
    b0=odata[i,5]
    b1=odata[i,6]
    b2=odata[i,7]
    b3=odata[i,8]
    b4=odata[i,9]
    
    qtild=np.linspace(1,float(peak[i]),50)
    tempr=float(temp[i])+460
    
    sl=np.log10(qtild)
    eps_dotco=10**(-((a0+del_sr)+(a1+a2*sl+a3*sl**2+a4*sl**3)/tempr))
    #assuming eps_dotc=eps_dotco
    #
    def sinh(x,a,s,n):
        return a*(np.sinh(x/s))**n
    guess=[0,2,1] #do not change, only needs to be changed if the sinh fit is not good.
    popt,pcov=optimize.curve_fit(sinh,qtild,eps_dotco,p0=guess,maxfev=100000)
    fig=plt.subplot(1,mats,count)
    count+=1
    plt.scatter (qtild,eps_dotco,label='Omega',color='red')
    plt.plot(qtild, sinh(qtild, popt[0], popt[1], popt[2]),label='Sinh',color='black')
    plt.xlabel(r'$\sigma$ (ksi)',fontsize=15) 
    plt.ylabel(r'$\dot\epsilon_{co}$',fontsize=15)
    plt.legend(loc='best',fontsize=13)
    plt.title(label=matname[i],fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(True)
    plt.savefig(project+'_materials.png')
    sindata=np.append(sindata,popt)
########################################################################################################################### 
x = dt.datetime.now()
date=x.strftime("%A, %d %B %Y")
subroutine.omega_model(project,date,mat,mats,matname,odata,ldata,wdata,sindata,delo,beta_omega,alpha)
subroutine.sinh_model(project,date,mat,mats,matname,odata,ldata,wdata,sindata,delo,beta_omega,alpha,peak,temp)