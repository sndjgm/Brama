# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 17:20:01 2022

@author: Guoqi Liu 642847452@qq.com
"""

import sys,os
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon,QTextCursor
from PyQt5.QtGui import QFont
from PyQt5.QtGui import QIntValidator,QDoubleValidator,QRegExpValidator
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal,QObject,QThread
import PyQt5.QtCore as qc
from PyQt5 import QtWidgets
from scipy import optimize
import time
import re
import csv
from math import *
import numpy as np
import pandas as pd
import xlwt,xlrd
import matplotlib.pyplot as plt
import logging
from xlutils.copy import copy
from mpl_toolkits.mplot3d import Axes3D
import statistics
import random
from scipy .interpolate import griddata
from scipy.stats import linregress



LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"

logging.basicConfig(filename='Brama.log', level=logging.INFO, format=LOG_FORMAT, datefmt=DATE_FORMAT)

def Cal_age(age):    
    a=exp(0.000000000155125*age*1000000)-1
    b=exp(0.00000000098485*age*1000000)-1
    c=exp(0.000000000049475*age*1000000)-1
    P382=1/137.88*(b/a)
            
    Q8=18.700-a*9.735
    R8=15.628-b*(9.735/137.88)            
    S8=38.630-c*36.837
    return R8,Q8,S8,P382,a,b,c







class Intercept_age:
    def __init__(self,A,Aerr,B,Berr,rho,method='W',normal=True):

        self.method=method
        self.normal=normal
        self.normal=True
        self.Aerr=A*Aerr/100
        self.Berr=B*Berr/100
        self.rho=rho
        self.A=A
        self.B=B
        self.Slope,self.Intercept,self.Slopeerr,self.Yinterr=self.SLOP_Cal(A,B)
        print('Slop:{0}+/-{2},Intercept:{1}+/-{3}'.format(self.Slope,self.Intercept,self.Slopeerr,self.Yinterr))  
        self.InitAge,self.InitAgeSig=self.Initial_Age()
        self.InitPbC,self.InitPbCSig=self.Initial_Pb()
        
        print('Age:{0}+/-{1},Pbc:{2}+/-{3}'.format(self.InitAge,self.InitAgeSig,self.InitPbC,self.InitPbCSig))  
        if normal==False:
            Slopi,Yint=self.SLOP_Cal(A,B)
            print(Slopi,Yint)
            Resid=B-Slopi*A-Yint
            Resid2=Resid*Resid
            s = 1.4826 * (1 + 5 / (A.size - 2)) * sqrt(np.median(Resid2))   
        
            xx=[]
            yy=[]
            ee=[]
            ff=[]
            gg=[]
            for i in range(A.size):
    
                u=abs(Resid[i]/s)
               
                if u<2.5:
                    xx.append(A[i])
                    yy.append(B[i])
                    ee.append(A[i]*Aerr[i]/100)
                    ff.append(B[i]*Berr[i]/100)
                    gg.append(rho[i])
            C=np.array(xx)
            D=np.array(yy)
            E=np.array(ee)
            F=np.array(ff)
            G=np.array(gg)
            self.Aerr=E
            self.Berr=F
            self.A=C
            self.B=D
            self.rho=G
            self.Slope,self.Intercept=self.SLOP_Cal(C,D)
            #reg=linregress(self.A, self.B)
            
            
            #self.Slope,self.Intercept=reg[0],reg[1]
            print(self.Slope,self.Intercept)  
            
            a,b=self.Err_Cal()
            print('Sloperro:{},Intercept erro:{}'.format(a,b))           
               
        
        
    def SLOP_Cal(self,A,B):
        
        Xinter=[]
        Slop=[]
        Yinter=[]
        
        for i in range(A.size-1):
            
            for j in range(i+1,A.size):
                if i!=j:
                    Vs =(B[j]-B[i])/(A[j]-A[i])+(0.5-random.random())*0.00000001
                    Slop.append(Vs)
                else:
                    Vs=0+(0.5-random.random())*0.00000001
                    Slop.append(Vs)
                #Slop.append(Vs)
                Vy = B[i]-Vs*A[i]+(0.5-random.random())*0.00000001
                Yinter.append(Vy)
                Xinter.append(-Vy / Vs)
        Slopi=statistics.median(Slop)  
        Yint=statistics.median(Yinter)
        #Slop=sorted(Slop)
        #Slopeerr=abs((Slop[int(0.75*(len(Slop)))]-Slop[int(0.25*(len(Slop)))])/2)
        #Slopeerr=statistics.pstdev(Slop)
        InA=sorted(Yinter)[int(0.95*(len(Yinter)))]
        InB=sorted(Yinter)[int(0.05*(len(Yinter)))]
        Yinterr=(InA-InB)/2
        ss=[]
        for inter in [InA, InB]:
            for i in range(A.size):
                ss.append((B[i]-inter)/A[i])
        Slopeerr=(sorted(ss)[int(0.95*(A.size*2))]-sorted(ss)[int(0.05*(A.size*2))])/2
        Xint = statistics.median(Xinter)  
        return Slopi,Yint,Slopeerr,Yinterr
        

    def Initial_Pb(self):
        self.InitPbC=1/(self.Slope*137.88)
        #self.InitPbCErr=self.InitPbC - (1 / (137.88 * (self.Slope + self.InitSlopeErr)))
        
        MSWD,n=self.Calmswd()
        self.InitPbCSig = self.Yinterr / (2 * sqrt(MSWD))
        return self.InitPbC,self.InitPbCSig
    def Initial_Age(self):
        
        self.InitAge= self.age_cal()[0]

       
        
        self.InitAgeErr= self.Slopeerr/self.Slope
        MSWD,n=self.Calmswd()
        self.InitAgeSig = self.InitAgeErr / (2 * sqrt(MSWD))
        return  self.InitAge,self.InitAgeSig 
    def Calmswd(self):
        if self.A.size>3:
            self.T=self.age_cal()[0]
            if self.T:              
                                   
                SWD=((self.B-self.Slope*self.A-self.Intercept)**2)/(self.Berr**2+self.Slope**2*self.Aerr**2-2*self.Slope*self.Berr*self.Aerr*self.rho)
                MSWD=(SWD.sum()/(self.A.size-2))

            else:
               MSWD=None  
        else:
            MSWD=None
        
        return MSWD,self.A.size
    def ConcSlope(self,T):
        
        Eterm=(0.000155125-0.00098485)*T
        if abs(Eterm)>1E+308:
            ConcSlop=False
        
        cs =  0.000155125* exp(Eterm) /0.00098485
        if self.normal:
            ConcSlope=cs
        else:
            e5=0.00098485*T
            if abs(e5)>1E+308: ConcSlop=False
                 
            temp = exp(e5) - 1 - (exp(0.00015512 * T) - 1) / cs
            ConcSlope = temp /137.88    
        
        return ConcSlope
    def ConcXage(self,R):
        if self.method=='TW':
            Lterm = 1 + 1 / R
            if Lterm<1E-307 or Lterm>1E+308:
                print('WRONG')
            Lterm=log(Lterm)
            ConcXage = Lterm / 0.000155125
        else:
            Lterm = 1 + R
            if Lterm<1E-307 or Lterm>1E+308:
                print('WRONG')
                Lterm=1E-100
                Lterm=log(Lterm)
                ConcXage = Lterm / 0.00098485  
            else:
                
                Lterm=log(Lterm)
            
                ConcXage = Lterm / 0.00098485   
            
        return  ConcXage 
        

        
    
    def ConcYage(self,R):
        if self.method=='TW':
            ConcYage=Age76Pb(R)
        else:
            Lterm = 1 + R
            
            if Lterm<1E-307 or Lterm>1E+308:
                print('WRONG')
            ConcYage = log(Lterm) / 0.000155125
    
          
        return  ConcYage 
        
    def ConcX(self,T):
         
        
        if self.method=='TW':
            Eterm =0.000155125*T
            if abs(Eterm)>1E+308:
              ConcX=False
              
            elif Eterm==0:
                ConcX = 1E+32
            else:
                ConcX =1/(exp(Eterm)-1)
        else:
            Eterm=0.00098485*T
            if abs(Eterm)>1E+308:
                ConcX=False
            
            ConcX =exp(Eterm)-1
         
        return ConcX
    
    def ConcY(self,T):
        
        e8=T*0.000155125
        if self.method=='TW':
            e5=T*0.00098485
            if abs(e5)>1E+308:
                ConcY=Falsee
            elif e5==0:
                ConcY=0.00098485/0.000155125/137.88
            else:
                ConcY=(exp(e5)-1)/(exp(e8)-1)/137.88
        else:
            if abs(e8)>1E+308:
                ConcY=False
            ConcY=exp(e8)-1
        
        return  ConcY 

    def f(self,TT):
        cs= self.ConcSlope(TT)   
        Numer = self.Intercept +cs* self.ConcX(TT) - self.ConcY(TT)
        Denom = (cs - self.Slope)
        x = Numer / Denom
        return x
    def ConcordiaCurve(self,output,age):
        
        plt.scatter(self.A,self.B,3,'red')

        plt.plot(self.A,self.Slope*self.A+self.Intercept)
        Start_t=max([self.A.min(),0])
        End_t=min([self.A.max(),367])
        X=np.arange(Start_t,End_t,(End_t-Start_t)/100).tolist()
        
        Y1=[]
        for i in X:
            
            Y1.append((i+1)**(0.000155125/0.00098485)-1)

        plt.plot(X,Y1)
        if age[0] and age[1]:
            plt.title(str(round(age[0],2))+'Ma'+' & '+str(round(age[1],2))+'Ma')
        else:
            if age[0] :
                plt.title(str(round(age[0],2))+'Ma')
            else:
                plt.title(str(round(age[1],2))+'Ma')
        plt.ylabel('206Pb/238U')
        plt.xlabel('207Pb/235U')
        for i in range(10):
            plt.text(X[i*10],Y1[i*10],str(int(log(X[i*10]+1)/0.00098485)), fontdict={'size': 6, 'color': 'g'})
        plt.savefig(output+'ConcordiaCurve',dpi=200, bbox_inches='tight', transparent=False)
        plt.show()
    def age_cal(self):      
        
        JustYoung=False
        JustOld=False
        if  JustYoung and JustOld:
            print('Error in call to ConcordiaIntercepts')
            
        age={} 
        
        for i in [0-JustOld,1+JustYoung]:
            j=i+1
            if j==1:
                TT=-500            
            else:             
                TT=5500        
            NoInter = 0 
            
            
            x = self.f(TT)
            
            if x<-1:
                NoInter = True
                age[i]=None
                print('NoInter')
                
            else:
                
                
                T = self.ConcXage(x)            
                Delta = abs(T - TT)
                while Delta>0.01:
                    TT=T
                    x = self.f(TT)   
                    T = self.ConcXage(x)            
                    Delta = abs(T - TT)
                    
                
            if not NoInter:
                ConcInt0=T
                age[i]=ConcInt0
            else:
                pass     
        
        return age
        
class Bayesian_Regression(QThread):
    
    signalForText=pyqtSignal(str)
    
    def write(self,text):
        self.signalForText.emit(str(text))
    def __init__(self,instrument=2,rawdata=None,name=None,       
                 Del206=0.025,Del207=0.04,Del232=0.025,Del238=0.025,Std76Bias=1,Std68Bias=0.7575,AblCorr=1,\
                 StartRow=120,EndRow=330,BslErr=0.04,ComPbFilter=False,Interpol=True,LowestUPb=4,bcg0=2,bcg1=10,MinU=50,MinCt=25,\
                     MaxSWD=10,N=5,Delspeak=False,Mult=2,InitAge=130,InitAgeErr=9.5,InitSlope=0.00753,InitSlopeErr=0.00052,WethMSWD=1,AgeIncr=1,PbCIncr=0.01,outputpath=None):        
        # Seting the initial parameter
        super(Bayesian_Regression,self).__init__()
        self.rawdata=rawdata
        if self.rawdata:            
            if name==None:
                self.filename=(list(rawdata.split('\\'))[-1]).split('.')[0]
            else:
                self.filename=name
        
        self.instrument=instrument
        self.EndRow=EndRow
        self.Del206=Del206
        self.Del207=Del207
        self.Del238=Del238
        self.Del232=Del232
        
        self.Corr76=1/Std76Bias
        self.Corr86 = Std68Bias/AblCorr
        self.LowestUPb=LowestUPb*Del238/Del206
        self.StartRow=StartRow
        self.bcg0=bcg0
        self.AblCorr=AblCorr
        self.bcg1=bcg1
        self.MinU=MinU
        self.MinCt=MinCt
        self.BslErr=BslErr
        self.ComPbFilter=ComPbFilter
        self.Interpol=Interpol

        self.InitAge=InitAge
        self.InitAgeErr=InitAgeErr
        self.InitSlope=InitSlope
        self.InitSlopeErr=InitSlopeErr
        self.WethMSWD=WethMSWD
        self.MaxSWD=MaxSWD
        self.N=int(N)
        self.Delspeak=Delspeak
        self.Mult=Mult
        self.InitPbC =1/(137.88*InitSlope)
        self.InitPbCErr = self.InitPbC - (1 / (137.88 * (self.InitSlope + self.InitSlopeErr)))
        
        self.InitAgeSig = self.InitAgeErr / (2 * sqrt(self.WethMSWD))
        self.InitPbCSig = self.InitPbCErr / (2 * sqrt(self.WethMSWD))
        
        self.AgeIncr=AgeIncr
        self.PbCIncr=PbCIncr
        
        self.MaxAge = False
        self.MaxPbC = False
        
        self.outputpath=outputpath
        
        self.L238 = 0.000155125
        self.L235 = 0.00098485
        self.Uic = 137.88
        
        if self.rawdata:                
            if self.outputpath:
                self.cwd=self.outputpath+'/'+self.filename+'_result'
            else:
                self.cwd=os.getcwd()+'/'+self.filename+'_result'
            if not os.path.exists(self.cwd):
                os.makedirs(self.cwd)
            self.start_time=time.time() 
            #print('Start  time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(self.start_time))))
            
            
            with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
                print('*'*50,file=f)
                print('#Start  time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(self.start_time))),file=f)  
                print('-'*50,file=f)
                print('#Instrument Parameters',file=f)
                print('Ms file:',self.rawdata,file=f)
                print('Background Start Row:{}\nBackground End Row:{}'.format(self.bcg0,self.bcg1),file=f)
                print('Single Start Row:{}\nSingle End Row:{}'.format(self.StartRow,self.EndRow),file=f)            
                print('#Dwell Times(s)\n206Pb:{0}\n207Pb:{1}\n238U:{3}\nOthers:{2} '.format(self.Del206,self.Del207,self.Del232,self.Del238),file=f)
                print('ComPbFilter:',self.ComPbFilter,file=f)
                print('Interpol:',self.Interpol,file=f)

    def pikes(self,A,multi):
        B=[]
        B.append(A[0])
        for i in range(A.size-1):
            #print(A[i])
            if A[i]>(((A[i-1]+A[i+1])*0.5)*multi) :
                B.append((A[i-1]+A[i+1])/2)
            else:
                B.append(A[i])
        B.append(A[-1])
        F=np.array(B) 
        return F  

    def flush(self):
        pass
    
  
        
    def couts(self):
        multi=self.Mult
        data_all=np.loadtxt(self.rawdata,dtype=str,delimiter=',',skiprows=self.instrument,comments='     ')[:self.EndRow]
        
        #data_all=np.loadtxt(self.rawdata,dtype=str,delimiter=',',skiprows=2,comments='     ')[:self.EndRow]
        col_data_name=data_all[0,:].tolist()
        self.data_Num=data_all.shape[0]
        col_name=[]
        for i in range(len(col_data_name)):
            col_name.append(''.join(list(filter(str.isdigit, col_data_name[i]))))
        if self.instrument==13:
            Time=data_all[2:,0].astype(float)
        else:
            Time=data_all[1:,0].astype(float)
        
        self.Tc=(Time[self.bcg1]-Time[self.bcg0])/(self.bcg1-self.bcg0)
        self.Ts=(self.Tc-(self.Del206+self.Del207+self.Del232+self.Del238))/(self.N-1)
        

        
        try:
            if self.instrument==13:
    

                if self.Delspeak :
                    Pb206=data_all[2:,col_name.index('206')].astype(float)*self.Del206
                    Pb207=data_all[2:,col_name.index('207')].astype(float)*self.Del207                   
                    U238=data_all[2:,col_name.index('238')].astype(float)*self.Del238
                    #Th232=data_all[1:,col_name.index('232')].astype(float)*self.Del232
                    Pb206=self.pikes(Pb206,multi)
                    Pb207=self.pikes(Pb207,multi)
                    U238=self.pikes(U238,multi)
                else:
                    Pb206=data_all[2:,col_name.index('206')].astype(float)*self.Del206
                    Pb207=data_all[2:,col_name.index('207')].astype(float)*self.Del207                   
                    U238=data_all[2:,col_name.index('238')].astype(float)*self.Del238
            else:
                if self.Delspeak :
                    Pb206=data_all[1:,col_name.index('206')].astype(float)*self.Del206
                    Pb207=data_all[1:,col_name.index('207')].astype(float)*self.Del207                   
                    U238=data_all[1:,col_name.index('238')].astype(float)*self.Del238
                    Pb206=self.pikes(Pb206,multi)
                    Pb207=self.pikes(Pb207,multi)
                    U238=self.pikes(U238,multi)
                else:
                    Pb206=data_all[1:,col_name.index('206')].astype(float)*self.Del206
                    Pb207=data_all[1:,col_name.index('207')].astype(float)*self.Del207                   
                    U238=data_all[1:,col_name.index('238')].astype(float)*self.Del238
                    #Th232=data_all[1:,col_name.index('232')].astype(float)*self.Del232
 
        
            self.Bsl206=Pb206[self.bcg0:self.bcg1].mean()*self.Del206
            self.Bsl207=Pb207[self.bcg0:self.bcg1].mean()*self.Del207
            self.Bsl238=U238[self.bcg0:self.bcg1].mean()*self.Del238
        
            Time_c=Time[self.StartRow:]
            Pb206_c=Pb206[self.StartRow:]
            Pb207_c=Pb207[self.StartRow:]
            U238_c=U238[self.StartRow:]
            TimeStamp=Time_c[0]
        except Exception as e:
            
            Time_c=np.zeros(1)
            Pb206=np.zeros(1)
            Th232=np.zeros(1)
            Pb207=np.zeros(1)
            U238=np.zeros(1)
            TimeStamp=Time_c[0]
            print(str(e))
            logging.info('An error occurred when loading data, please check the data format.%s',str(e))
        
        U238_Int=[]
        Pb206_Int=[]
        Pb207_Int=[]
        Time_Int=[]
        
        
        for i in range(1,len(Time_c)):
            
            if self.Interpol:
                
                U238Int=U238_c[i-1] + (U238_c[i] -U238_c[i-1]) * (Time_c[i]-Time_c[i-1] - (self.Tc - (self.Del238 + self.Del206) / 2 - self.Ts)) / (Time_c[i]-Time_c[i-1])-self.Bsl238
               
                Pb207Int=Pb207_c[i-1]+ (Pb207_c[i]-Pb207_c[i-1]) * (Time_c[i]-Time_c[i-1]-((self.Del207 + self.Del206)/2+self.Ts))/(Time_c[i]-Time_c[i-1])-self.Bsl207
                Pb206Int =Pb206_c[i]-self.Bsl206
                Time_Int.append(Time_c[i-1])

                    
                
                
                Pb206_Int.append(Pb206Int)
                Pb207_Int.append(Pb207Int)
                U238_Int.append(U238Int)
                
                
                
            else:
                Pb206Int =Pb206_c[i-1]-self.Bsl206              
                Pb207Int =Pb207_c[i-1]-self.Bsl207
                U238Int = U238_c[i-1]-self.Bsl238
                Time_Int.append(Time_c[i-1])
                Pb206_Int.append(Pb206Int)
                Pb207_Int.append(Pb207Int)
                U238_Int.append(U238Int)         
        data_couts=[]
        for i in range(len(U238_Int)):
            if U238_Int[i]>=self.MinU:
                data_couts.append([Time_Int[i],Pb206_Int[i],Pb207_Int[i],U238_Int[i]])
            else:
                pass
          
        data_couts.append([None]*4)
        self.data_couts=np.array(data_couts)
        return self.data_couts
    def fliterdata(self,logfile=None):
        
        
        JJ=[]
        J206={}
        J207={}
        J238={}
        JK={}
        JS206={}
        JS207={}
        JS238={}
      

        
        
        i=0
        j=1
        K=1
        NumLines = 0
        NumData = 1

        workbook = xlwt.Workbook()
        worksheet=workbook.add_sheet('U-Pb Isotope results')
        style = xlwt.XFStyle() # ?????
        font = xlwt.Font() # ???????
        font.name = 'Times New Roman' 
        font.bold = True # ??
        font.underline = True # ???
        font.italic = True # ???
        style.font = font # ???? 
        worksheet.write(0, 0, 'No.')
        worksheet.write(0, 1, 'Time(s)') 
        worksheet.write(0, 2, 'XX') 
        worksheet.write(0, 3, 'YY') 
        worksheet.write(0, 4, 'Pb206')
        worksheet.write(0, 5, 'Pb207')
        worksheet.write(0, 6, 'U238')
        worksheet.write(0, 7,'Num.')
        worksheet.write(0, 8,'Sig206')
        worksheet.write(0, 9,'Sig207')
        worksheet.write(0, 10,'Sig238')
        worksheet.write(0, 11,'Sq Wt Dev')
        worksheet.write(0, 12,'Ed. Lines')
        
        
        worksheet.write(0, 15,'238U/206Pb')
        worksheet.write(0, 16,'Sig8/6(%)')
        worksheet.write(0, 17,'207Pb/206Pb')
        worksheet.write(0, 18,'Sig7/6(%)')
        worksheet.write(0, 19,'Rho')
        
        worksheet.write(0, 21,'207Pb/235U')
        worksheet.write(0, 22,'Sig7/5(%)')
        worksheet.write(0, 23,'206Pb/238U')
        worksheet.write(0, 24,'Sig6/8(%)')
        worksheet.write(0, 25,'Rho')
        
        worksheet.write(0, 27,'Age(Ma)')
        worksheet.write(0, 28,'206Pb/238U')
        worksheet.write(0, 29,'Err(1s)')
        worksheet.write(0, 30,'207Pb/235U')
        worksheet.write(0, 31,'Rrr(1s)')
        result_all=[]
        wetherilldata=[]
        
        while self.data_couts[i,1]!=None:                
                          
            TimeStamp=self.data_couts[i,0]
            NextTimeStamp = TimeStamp
        
            Pb206=self.data_couts[i,1]
            Pb207=self.data_couts[i,2]
            U238=self.data_couts[i,3] 
            
            if self.ComPbFilter:
                i=i+1
                if self.data_couts[i-1,1]<U238/self.LowestUPb and self.data_couts[i-1,2]<U238/self.LowestUPb:
                    
                    if self.data_couts[i-1,1]<self.MinCt or self.data_couts[i-1,2]<self.MinCt:
                    
                        K=K+1
                        NextTimeStamp=self.data_couts[i,0]
                        NextPb206=self.data_couts[i,1]
                        NextPb207=self.data_couts[i,2]
                        NextU238 =self.data_couts[i,3]
                        
                        if NextPb206!=None:
                            
                            self.data_couts[i,1]=Pb206+NextPb206
                            self.data_couts[i,2]=Pb207+NextPb207
                            self.data_couts[i,3]=U238+NextU238
                        else :
                            pass
                        
                    
                    else:
                        TimeStamp= (TimeStamp + NextTimeStamp + self.Tc) / 2
                        NumData = NumData + K
                        
                        Sig206 = sqrt(Pb206 + (K + 1) * self.Bsl206)
                        Sig207 = sqrt(Pb207 + (K + 1) * self.Bsl207)
                        Sig238=sqrt(U238 + (K + 1) * self.Bsl238)
                        SqrtK = sqrt(K)
                        if self.BslErr > 0:
                            Sig206 = sqrt(Sig206*Sig206 + (self.BslErr * Pb206 / SqrtK) *(self.BslErr * Pb206 / SqrtK))
                            Sig207 = sqrt(Sig207*Sig207 + (self.BslErr * Pb207 / SqrtK) *(self.BslErr * Pb207 / SqrtK))
                            Sig238 = sqrt(Sig238*Sig238 + (self.BslErr * U238 / SqrtK) *(self.BslErr * U238 / SqrtK))
                        worksheet.write(j, 0, NumLines+1)    
                        worksheet.write(j, 1, TimeStamp)
                        worksheet.write(j, 4, Pb206)
                        worksheet.write(j, 5, Pb207)
                        worksheet.write(j, 6, U238)
                        worksheet.write(j, 7,K)
                        worksheet.write(j, 8,Sig206)
                        worksheet.write(j, 9,Sig207)
                        worksheet.write(j, 10,Sig238)
                        if logfile and self.data_couts[i-1,0]:
                            workbooklog=pd.read_excel(logfile,sheet_name='Sheet1')
                            
                            
                            XX0=list(workbooklog[workbooklog['FileName']==self.filename]['X0'])[0]
                            YY0=list(workbooklog[workbooklog['FileName']==self.filename]['Y0'])[0]
                            XX1=list(workbooklog[workbooklog['FileName']==self.filename]['X1'])[0]
                            YY1=list(workbooklog[workbooklog['FileName']==self.filename]['Y1'])[0]
                            if XX0==XX1 or YY1==YY0:
                                if XX0==XX1:
                                    XX=XX0
                                if YY1==YY0:
                                    YY=YY0                           
                            else: 
                                XX=XX0+(self.data_couts[i,0]-self.data_couts[0,0])*(XX1-XX0)/(self.data_couts[-2,0]-self.data_couts[0,0])                               
                            
                                YY=YY0+(self.data_couts[i,0]-self.data_couts[0,0])*(YY1-YY0)/(self.data_couts[-2,0]-self.data_couts[0,0])
                            
                            
                            worksheet.write(j, 2,XX)
                            worksheet.write(j, 3,YY)
                        
                        
                        PC86Err=100*sqrt(pow((Sig206/Pb206),2)+(pow((Sig238/U238),2)))
                        PC76Err=100*sqrt(pow((Sig206/Pb206),2)+(pow((Sig207/Pb207),2)))
                        PC75Err=100*sqrt(pow((Sig238/U238),2)+(pow((Sig207/Pb207),2)))
                        RhoTW=(PC86Err*PC86Err+PC76Err*PC76Err-PC76Err*PC76Err)/(2*PC86Err*PC75Err)
                        
                        RhoWeth=(PC86Err*PC86Err+PC75Err*PC75Err-PC75Err*PC75Err)/(2*PC86Err*PC76Err)
                        worksheet.write(j, 15,self.Corr86 * (self.Del206 / self.Del238) * U238 / Pb206)
                        worksheet.write(j, 16,PC86Err)
                        worksheet.write(j, 17,self.Corr76*(self.Del206/self.Del207)*Pb207/Pb206)
                        worksheet.write(j, 18,PC76Err)
                        worksheet.write(j, 19,RhoTW)
                        
                        worksheet.write(j, 21,(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)
                        worksheet.write(j, 22,PC75Err)
                        worksheet.write(j, 23,(self.Del238/self.Del206)*Pb206/U238/self.Corr86)
                        worksheet.write(j, 24,PC86Err)
                        worksheet.write(j, 25,RhoWeth)
                        
                        worksheet.write(j, 28,log(1+(self.Del238/self.Del206)*Pb206/U238/self.Corr86)/0.000000000155125/1000000)
                        worksheet.write(j, 29,(PC86Err/100)*log(1+(self.Del238/self.Del206)*Pb206/U238/self.Corr86)/0.000000000155125/1000000)
                        worksheet.write(j, 30,log(1+(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)/0.00000000098485/1000000)
                        worksheet.write(j, 31,(PC75Err/100)*log(1+(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)/0.00000000098485/1000000)
                        
                        result=[NumLines+1,TimeStamp,Pb206,Pb207,U238,Sig206,Sig207,Sig238,self.Corr86 * (self.Del206 / self.Del238) * U238 / Pb206,\
                                PC86Err,self.Corr76*(self.Del206/self.Del207)*Pb207/Pb206,PC76Err,RhoTW,(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238\
                                    ,PC75Err,(self.Del238/self.Del206)*Pb206/U238/self.Corr86,PC86Err,RhoWeth]                   
                        
                        
                        
                        
                        
                        K = 1
                        j=j+1
                        NumLines = NumLines + 1
                        result_all.append(result)
                        wetherilldata.append(result[13:])
                       
            
            else:
            
                i=i+1
                if self.data_couts[i-1,1]<self.MinCt or self.data_couts[i-1,2]<self.MinCt:
                    
                    K=K+1
                    NextTimeStamp=self.data_couts[i,0]
                    NextPb206=self.data_couts[i,1]
                    NextPb207=self.data_couts[i,2]
                    NextU238 =self.data_couts[i,3]
                    if NextPb206!=None:
                        self.data_couts[i,1]=Pb206+NextPb206
                        self.data_couts[i,2]=Pb207+NextPb207
                        self.data_couts[i,3]=U238+NextU238
                    else:
                        pass
                    
                
                else:
                    TimeStamp= (TimeStamp + NextTimeStamp + self.Tc) / 2
                    NumData = NumData + K
                    Sig206 = sqrt(Pb206 + (K + 1) * self.Bsl206)
                    Sig207 = sqrt(Pb207 + (K + 1) * self.Bsl207)
                    Sig238=sqrt(U238 + (K + 1) * self.Bsl238)
                    SqrtK = sqrt(K)
                    if self.BslErr > 0:
                        Sig206 = sqrt(Sig206*Sig206 + (self.BslErr*Pb206 / SqrtK) *(self.BslErr * Pb206 / SqrtK))
                        Sig207 = sqrt(Sig207*Sig207 + (self.BslErr*Pb207 / SqrtK) *(self.BslErr * Pb207 / SqrtK))
                        Sig238 = sqrt(Sig238*Sig238 + (self.BslErr*U238 / SqrtK) *(self.BslErr * U238 / SqrtK))
                    worksheet.write(j, 0, NumLines+1)     
                    worksheet.write(j, 1, TimeStamp)
                    worksheet.write(j, 4, Pb206)
                    worksheet.write(j, 5, Pb207)
                    worksheet.write(j, 6, U238)
                    worksheet.write(j, 7,K)
                    worksheet.write(j, 8,Sig206)
                    worksheet.write(j, 9,Sig207)
                    worksheet.write(j, 10,Sig238)
                    
                    if logfile and self.data_couts[i-1,0]:
                        workbooklog=pd.read_excel(logfile,sheet_name='Sheet1')
                        
                        XX0=list(workbooklog[workbooklog['FileName']==self.filename]['X0'])[0]
                        YY0=list(workbooklog[workbooklog['FileName']==self.filename]['Y0'])[0]
                        XX1=list(workbooklog[workbooklog['FileName']==self.filename]['X1'])[0]
                        YY1=list(workbooklog[workbooklog['FileName']==self.filename]['Y1'])[0]
                        try:
                            XX=XX0+(self.data_couts[i,0]-self.data_couts[0,0])*(XX1-XX0)/(self.data_couts[-2,0]-self.data_couts[0,0])
                            YY=YY0+(self.data_couts[i,0]-self.data_couts[0,0])*(YY1-YY0)/(self.data_couts[-2,0]-self.data_couts[0,0])
                        except TypeError:
                            XX=XX1
                            YY=YY1
                            
                        worksheet.write(j, 2,XX)
                        worksheet.write(j, 3,YY)
                        
                    
                    
                    PC86Err=100*sqrt(pow((Sig206/Pb206),2)+(pow((Sig238/U238),2)))
                    PC76Err=100*sqrt(pow((Sig206/Pb206),2)+(pow((Sig207/Pb207),2)))
                    PC75Err=100*sqrt(pow((Sig238/U238),2)+(pow((Sig207/Pb207),2)))
                    RhoTW=(PC86Err*PC86Err+PC76Err*PC76Err-PC76Err*PC76Err)/(2*PC86Err*PC75Err)
                    
                    RhoWeth=(PC86Err*PC86Err+PC75Err*PC75Err-PC75Err*PC75Err)/(2*PC86Err*PC76Err)
                    
                    worksheet.write(j, 15,self.Corr86 * (self.Del206 / self.Del238)*U238/Pb206)
                    worksheet.write(j, 16,PC86Err)
                    worksheet.write(j, 17,self.Corr76*(self.Del206/self.Del207)*Pb207/Pb206)
                    worksheet.write(j, 18,PC76Err)
                    worksheet.write(j, 19,RhoTW)
                    
                    worksheet.write(j, 21,(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)
                    worksheet.write(j, 22,PC75Err)
                    worksheet.write(j, 23,(self.Del238/self.Del206)*Pb206/U238/self.Corr86)
                    worksheet.write(j, 24,PC86Err)
                    worksheet.write(j, 25,RhoWeth)      
                    
                    worksheet.write(j, 28,log(1+(self.Del238/self.Del206)*Pb206/U238/self.Corr86)/0.000000000155125/1000000)
                    worksheet.write(j, 29,(PC86Err/100)*log(1+(self.Del238/self.Del206)*Pb206/U238/self.Corr86)/0.000000000155125/1000000)
                    worksheet.write(j, 30,log(1+(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)/0.00000000098485/1000000)
                    worksheet.write(j, 31,(PC75Err/100)*log(1+(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238)/0.00000000098485/1000000)
                    
                    result=[NumLines+1,TimeStamp,Pb206,Pb207,U238,Sig206,Sig207,Sig238,self.Corr86 * (self.Del206 / self.Del238) * U238 / Pb206,\
                            PC86Err,self.Corr76*(self.Del206/self.Del207)*Pb207/Pb206,PC76Err,RhoTW,(self.Corr76/self.Corr86)*self.Uic*(self.Del238/self.Del207)*Pb207/U238\
                                ,PC75Err,(self.Del238/self.Del206)*Pb206/U238/self.Corr86,PC86Err,RhoWeth]                   
                        
                    
                    
                    K = 1
                    j=j+1
                    NumLines = NumLines + 1
                    result_all.append(result)
                    wetherilldata.append(result[13:])
        
        
        
        
        '''
        t='2023-12-6 18:34:00'
        s_t=time.strptime(t,'%Y-%m-%d %H:%M:%S')
        mkt=int(time.mktime(s_t))
        
        if np.sin(15*pi*mkt**(10))*np.exp(0*mkt**(6)-3)*np.cos(2*pi*mkt**(9))==0.037272471111608585:
         '''              
        workbook.save(self.cwd+'/'+self.filename+'_'+'Result.xls')
        
        
        

        

        #aerr,berr=intercept_age.Slop_Intercept_err()
        try:
            outputwetherill=np.array(wetherilldata)
       

            intercept_age=Intercept_age(outputwetherill[:,0],outputwetherill[:,1],outputwetherill[:,2],outputwetherill[:,3],outputwetherill[:,4])
        
            age=intercept_age.age_cal()
            MSWD,Number=intercept_age.Calmswd()    
            
            intercept_age.ConcordiaCurve(self.cwd+'/'+self.filename,age)
            self.Agestart,self.Agesig=intercept_age.Initial_Age()
            self.Pbstart,self.Pbsig=intercept_age.Initial_Pb()
            print('Age(Wetherill)={}+/-{}'.format(self.Agestart,self.Agesig))
            print('Pbc(Wetherill)={}+/-{}'.format(self.Pbstart,self.Pbsig))
            if MSWD:
                #print('slop:+/-{}\nIntercept:+/-{}'.format(aerr,berr))
                print('MSWD={:.1f}(N={})'.format(MSWD,Number))
            else:
                print('MSWD=',MSWD)
            with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
                
                print('#Number of Data:{}'.format(NumLines),file=f)
                print('#Raw Data:{}'.format(NumData),file=f)
                print('#Compressed Data:{:.1%}'.format( 1 - NumLines / NumData),file=f)
                           
                print('Baseline  Error:{:.1%}'.format(self.BslErr),file=f)
                print('Min. U  (cps):{}'.format(self.MinU),file=f)
                print('Min. U Pb (counts):{}'.format(self.MinCt),file=f)
                print('Std Bias 207/206:{:.4f}'.format(1/self.Corr76),file=f)
                print('Std Bias 206/238:{:.4f}'.format(self.Corr86*self.AblCorr),file=f)
                print('Alblatio Corr. factors:{:.4f}'.format(self.AblCorr),file=f)
                if age[0]==None or age[1]==None:
                    if age[0]==None and age[1]:
                        print('No Lower Intercept Age',file=f)
                        print('Upper Intercept Age:{:.2f} Ma.'.format(age[1]),file=f)
                    elif age[1]==None and age[0]:
                        print('No Upper Intercept Age',file=f)
                        print('Lower Intercept Age:{:.2f} Ma.'.format(age[0]),file=f)
                else:                
                    print('Lower Intercept Age:{:.2f} Ma;\nUpper Intercept Age:{:.2f} Ma.'.format(age[0],age[1]),file=f)


                print('-'*50,file=f)
        except Exception as e:
            print(str(e))

            MSWD=0
            Number=0
            self.Agestart=110
            self.Agesig=5
            self.Pbstart=0.85
            self.Pbsig=0.2
        
        
        if self.data_couts[0][0]:
            
            print('-'*28)
            print('#Number of Data:{}'.format(NumLines))
            print('#Raw Data:{}'.format(NumData))
            print('#Compressed Data:{:.1%}'.format( 1 - NumLines / NumData))
            print('-'*28)
            
            print('Baseline  Error:{:.1%}'.format(self.BslErr))
            print('Min. U (cps):{}'.format(self.MinU))
            print('Min. U Pb (counts):{}'.format(self.MinCt))
            print('Std Bias 207/206:{:.4f}'.format(1/self.Corr76))
            print('Std Bias 206/238:{:.4f}'.format(self.Corr86*self.AblCorr))
            print('Alblatio Corr. factors:{:.4f}'.format(self.AblCorr))
            
            
            
            self.result=np.array(result_all).astype(float)
        else:
           print('Unknown error!')
        
        try:
            return self.Agestart,self.Agesig,self.Pbstart,self.Pbsig
        except:
            return 0,0,0,0


        
    def XintValue(self,PeakAge, PeakPbC):
        XConcVal = 1 / (exp(self.L238 * PeakAge) - 1) 
        YConcVal = (exp(self.L235 * PeakAge) - 1) / (exp(self.L238 * PeakAge) - 1) / self.Uic 
        XintValue = PeakPbC * XConcVal / (PeakPbC - YConcVal)
        return XintValue
    
    def  regress(self,Pb76Start=0.8132,Pb76End=1.115,AgeStart=170,AgeEnd=216.2,DelPb76=0.003,DelAge=0.46,pltshowT=False):

        xint_all=[]
        NumAges = int(abs(AgeEnd - AgeStart) / DelAge) + 1
        NumYints =int(abs(Pb76End - Pb76Start) / DelPb76) + 1
        
        LogProb_all=[]
        
        self.Ages=list(np.linspace(AgeStart,AgeEnd,NumAges))
        
        
        self.Pb76=list(np.linspace(Pb76Start,Pb76End,NumYints))
        
        print('Age start:{0}\nAge end:{1}\nNum. Ages:{2} '.format(AgeStart,AgeEnd,NumAges))
        print('Pbc76 start:{0}\nPbc76 end:{1}\nNum. Pbc76:{2} '.format(Pb76Start,Pb76End,NumYints))
        #print(list(np.linspace(AgeStart,AgeEnd,NumAges)))
        #print(list(np.linspace(Pb76Start,Pb76End,NumYints)))
        with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
            print('-'*50,file=f)
            print('#Regression Parameters',file=f)
            print('Age start:{0}\nAge end:{1}\nNum. Ages:{2} '.format(AgeStart,AgeEnd,NumAges),file=f)
            print('Pbc76 start:{0}\nPbc76 end:{1}\nNum. Pbc76:{2} '.format(Pb76Start,Pb76End,NumYints),file=f)
            print('-'*50,file=f)

        #numi=0
        
        #print("\r", end="")
        #print(" {}%: ".format(numi), "█" * (numi // 2), end="",flush=True)
        
        for AgeVal in self.Ages:
           
            #numi+=1
            for Pb76Val in self.Pb76: 
                Xint_Ix_Jy =self.XintValue(AgeVal,Pb76Val)
                LogProb=self.RPD(AgeVal,Pb76Val)
                LogProb_all.append(LogProb)
                xint_all.append(Xint_Ix_Jy)
         
        print()        
        MaxLogProb=max(LogProb_all)  
        self.xint_all=np.array(xint_all).reshape(NumAges,NumYints)
        LogProb_all=np.power(10,np.array(LogProb_all)-MaxLogProb).reshape(NumAges,NumYints)
        self.LogProb_all=LogProb_all/LogProb_all.max()
        
        
        self.AgeProbSum=self.LogProb_all.sum(axis=1)
        self.YintProbSum=self.LogProb_all.sum(axis=0)
        
        
        self.PeakAge=self.Ages[self.AgeProbSum.argmax()]
        self.PeakPbC=self.Pb76[self.YintProbSum.argmax()]
        
                    
        
        
        TailProb1=0
        TailProb2=0
        
        for i in range(self.AgeProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2 = TailProb1 +self.AgeProbSum[i]

            if TailProb2>=0.025*self.AgeProbSum.sum():
                AbsAgeErrLow=self.Ages[i]
                break
            else:
                pass
        TailProb1=0
        TailProb2=0
        
        for i in range(self.AgeProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2 = TailProb1 +self.AgeProbSum[i]

            if TailProb2>=0.5*self.AgeProbSum.sum():
                MedianAge=self.Ages[i]
                break
            else:
                pass        
        TailProb1=0
        TailProb2=0
            
        for i in range(self.AgeProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2 = TailProb1 +self.AgeProbSum[i]
            if TailProb2>=0.975*self.AgeProbSum.sum():
                AbsAgeErrHigh=self.Ages[i]
                break
            else:
                pass
                
        
        TailProb1=0
        TailProb2=0
        a=0
        b=0
        for i in range(self.YintProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2=TailProb1+self.YintProbSum[i]
             
            
            if TailProb2>=0.025*self.YintProbSum.sum():
                AbsYintErrLow=self.Pb76[i]
                break
        TailProb1=0
        TailProb2=0
        for i in range(self.YintProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2=TailProb1+self.YintProbSum[i]
                     
            if TailProb2>=0.975*self.YintProbSum.sum():
                AbsYintErrHigh=self.Pb76[i]
                break
            else:
                pass
        TailProb1=0
        TailProb2=0 
        for i in range(self.YintProbSum.shape[0]):
            TailProb1=TailProb2
            TailProb2=TailProb1+self.YintProbSum[i]
            if TailProb2>=0.5*self.YintProbSum.sum():
                MedianPbC=self.Pb76[i]
                break
            else:
                pass
            
        XConcVal=1/(exp(self.L238*self.PeakAge)-1) 
        YConcVal=(exp(self.L235 *self.PeakAge)-1)/(exp(self.L238*self.PeakAge)- 1)/self.Uic
        PeakXint=self.PeakPbC*XConcVal/(self.PeakPbC-YConcVal)
        NumSWDRejects = 0
        SWD_all=[]
        for i in range(self.Datum.shape[0]):
            Pb206Y_Datum=self.Datum[i,0]
            Pb207Z_Datum=self.Datum[i,1]
            U238X_Datum=self.Datum[i,2]
            SigY_Datum=self.Datum[i,3]
            SigZ_Datum=self.Datum[i,4]
            SigX_Datum=self.Datum[i,5]

            
            Kst = Pb206Y_Datum * self.PeakPbC - U238X_Datum * self.PeakPbC /PeakXint - Pb207Z_Datum
            A = (SigZ_Datum/SigX_Datum)**2+(self.PeakPbC/PeakXint)**2
            B = -2*(self.PeakPbC**2)/PeakXint
            C = (SigZ_Datum/SigY_Datum)**2+self.PeakPbC**2
            D = -2*Kst*self.PeakPbC/PeakXint
            E = 2*Kst*self.PeakPbC
                
            Disc =B**2-4*A*C
                   
            f = sqrt((Kst**2-(C*D**2+A* E**2-B*E*D)/Disc)/SigZ_Datum**2)
            SWD=(f**2)/4
            if SWD > self.MaxSWD:
                NumSWDRejects=NumSWDRejects+1
            
            SWD_all.append(SWD)
        UsedDataNum=self.Datum.shape[0]-NumSWDRejects
        self.SWD=np.array(SWD_all) 
        
        
        
        result_excele=xlrd.open_workbook(self.cwd+'/'+self.filename+'_'+'Result.xls',formatting_info=True)
        wb = copy(result_excele)
        ws = wb.get_sheet(0)
        
        
        for i in range(self.SWD.shape[0]):
            ws.write(i+1,11,self.SWD[i])                
               
        wb.save(self.cwd+'/'+self.filename+'_'+'Result.xls')
        
        
        MSWD_S=self.SWD[self.SWD<self.MaxSWD]

        Datum_selected=np.column_stack((self.Datum,self.SWD))
        self.Datum_selected=Datum_selected[np.where(Datum_selected[:,-1]<self.MaxSWD)]
        #np.savetxt(self.cwd+'/'+self.filename+'_output(CPS).csv',self.Datum_selected,delimiter=',')
        
        
        Datum_selected=self.Datum_selected
        try:
            
            Pb206o=Datum_selected[:,0]
            Pb207o=Datum_selected[:,1]
            U238o=Datum_selected[:,2]
            Sig206o=Datum_selected[:,3]
            Sig207o=Datum_selected[:,4]
            Sig238o=Datum_selected[:,5]
            mswdo=Datum_selected[:,6]
        except Exception as e:
            logging.error('Unknown erro:%s',str(e))
            print(str(e))    
            

        
        
        #wetherilldata=[]
        outputdata=[]
        for i in range(Datum_selected.shape[0]):                         
            PC86Err=100*sqrt(pow((Sig206o[i]/Pb206o[i]),2)+(pow((Sig238o[i]/U238o[i]),2)))
            PC76Err=100*sqrt(pow((Sig206o[i]/Pb206o[i]),2)+(pow((Sig207o[i]/Pb207o[i]),2)))
            PC75Err=100*sqrt(pow((Sig238o[i]/U238o[i]),2)+(pow((Sig207o[i]/Pb207o[i]),2)))
            RhoTW=(PC86Err*PC86Err+PC76Err*PC76Err-PC76Err*PC76Err)/(2*PC86Err*PC75Err)
            
            RhoWeth=(PC86Err*PC86Err+PC75Err*PC75Err-PC75Err*PC75Err)/(2*PC86Err*PC76Err)
            #wetherilldata.append([137.88*Pb207o[i]/U238o[i],PC75Err,Pb206o[i]/U238o[i],PC86Err,RhoWeth])
            outputdata.append([i+1,Pb206o[i],Pb207o[i],U238o[i],Sig206o[i],Sig207o[i],Sig238o[i],mswdo[i],'',U238o[i]/Pb206o[i],PC86Err,Pb207o[i]/Pb206o[i],PC76Err,RhoTW,
                               '',137.88*Pb207o[i]/U238o[i],PC75Err,Pb206o[i]/U238o[i],
                               PC86Err,RhoWeth,'',log(1+Pb206o[i]/U238o[i])/0.000000000155125/1000000,
                               (PC86Err/100)*(log(1+Pb206o[i]/U238o[i])/0.000000000155125/1000000),
                               log(1+137.88*Pb207o[i]/U238o[i])/0.00000000098485/1000000,
                               (PC75Err/100)*log(1+137.88*Pb207o[i]/U238o[i])/0.00000000098485/1000000])
        
        '''outputwetherill=np.array(wetherilldata)
       
        try:
            intercept_age=Intercept_age(outputwetherill[:,0],outputwetherill[:,1],outputwetherill[:,2],outputwetherill[:,3],outputwetherill[:,4])
            age=intercept_age.age_cal()
            MSWD,Number=intercept_age.Calmswd()
            intercept_age.ConcordiaCurve(self.cwd+'/'+self.filename,age)
            if MSWD:
                print('MSWD={:.1f}(N={})'.format(MSWD,Number))
            else:
                print('MSWD=',MSWD)
        except:
            age=[None,None]
        '''    
        workbook = xlwt.Workbook()
        worksheet=workbook.add_sheet('regress results')
        style = xlwt.XFStyle() # ?????
        font = xlwt.Font() # ???????
        font.name = 'Times New Roman' 
        font.bold = True # ??
        font.underline = True # ???
        font.italic = True # ???
        style.font = font # ???? 
        worksheet.write(0, 0, 'No.')
       
        worksheet.write(0, 1, 'Pb206')
        worksheet.write(0, 2, 'Pb207')
        worksheet.write(0, 3, 'U238')

        worksheet.write(0, 4,'Sig206')
        worksheet.write(0, 5,'Sig207')
        worksheet.write(0, 6,'Sig238')
        worksheet.write(0, 7,'Sq Wt Dev')
        
        
       
        worksheet.write(0, 8,'Tera-Wasserburg')
        worksheet.write(0, 9,'238U/206Pb')
        worksheet.write(0, 10,'Sig8/6(%)')
        worksheet.write(0, 11,'207Pb/206Pb')
        worksheet.write(0, 12,'Sig7/6(%)')
        worksheet.write(0, 13,'Rho')
        
        
        worksheet.write(0, 14,'Wetherill')
        worksheet.write(0, 15,'207Pb/235U')
        worksheet.write(0, 16,'Sig7/5(%)')
        worksheet.write(0, 17,'206Pb/238U')
        worksheet.write(0, 18,'Sig6/8(%)')
        worksheet.write(0, 19,'Rho')
       
        worksheet.write(0, 20,'Age(Ma)')
        worksheet.write(0, 21,'206Pb/238U')
        worksheet.write(0, 22,'Err(1s)')
        worksheet.write(0, 23,'207Pb/235U')
        worksheet.write(0, 24,'Rrr(1s)')     
        for row in range(len(outputdata)):
            for col in range(25):
                worksheet.write(row+1,col,outputdata[row][col])
        
        workbook.save(self.cwd+'/'+self.filename+'_regress_'+'Result.xls')
        
        
        
        #np.savetxt(self.cwd+'/'+self.filename+'_output.csv',self.Datum,delimiter=',')
        self.MSWD=MSWD_S.sum()/UsedDataNum
        with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
            print('-'*50,file=f)
            print('#Regression Results',file=f)
            print('Peak age:{0:.1f} +{2:.2f}/-{1:.2f} Ma'.format(self.PeakAge,self.PeakAge-AbsAgeErrLow,AbsAgeErrHigh-self.PeakAge),file=f)                   
            print('Peak PbC:{0:.4f} +{2:.4f}/-{1:.4f} '.format(self.PeakPbC,self.PeakPbC-AbsYintErrLow,AbsYintErrHigh-self.PeakPbC),file=f)
            print('MSWD={:.2f} (Num.Data={})'.format(self.MSWD,UsedDataNum),file=f)                  
            print('Mean age:{0:.1f} +{2:.2f}/-{1:.2f} Ma'.format(MedianAge,MedianAge-AbsAgeErrLow,AbsAgeErrHigh-MedianAge),file=f)           
            print('Mean PbC:{0:.4f} +{2:.4f}/-{1:.4f} '.format(MedianPbC,MedianPbC-AbsYintErrLow,AbsYintErrHigh-MedianPbC),file=f) 
            '''if age[0]==None or age[1]==None:
                if age[0]==None and age[1]:
                    print('No Lower Intercept Age',file=f)
                    print('Upper Intercept Age:{:2f} Ma.'.format(age[1]),file=f)
                elif age[1]==None and age[0]:
                    print('No Upper Intercept Age',file=f)
                    print('Lower Intercept Age:{:2f} Ma.'.format(age[0]),file=f)
            else:                
                print('Lower Intercept Age:{:.2f} Ma;\nUpper Intercept Age:{:.2f} Ma.'.format(age[0],age[1]),file=f)
            '''
            print('-'*50,file=f)
            
        print()  
        

        print('Peak age:{0:.1f} +{2:.2f}/-{1:.2f} Ma'.format(self.PeakAge,self.PeakAge-AbsAgeErrLow,AbsAgeErrHigh-self.PeakAge))
                   
        print('Peak PbC:{0:.4f} +{2:.4f}/-{1:.4f} '.format(self.PeakPbC,self.PeakPbC-AbsYintErrLow,AbsYintErrHigh-self.PeakPbC))
        print('MSWD={:.2f} (Num.Data={})'.format(self.MSWD,UsedDataNum))
        print('\n')
        
        print('Mean age:{0:.1f} +{2:.2f}/-{1:.2f} Ma'.format(MedianAge,MedianAge-AbsAgeErrLow,AbsAgeErrHigh-MedianAge))
           
        print('Mean PbC:{0:.4f} +{2:.4f}/-{1:.4f} '.format(MedianPbC,MedianPbC-AbsYintErrLow,AbsYintErrHigh-MedianPbC))
        '''if age[0]==None or age[1]==None:
            if age[0]==None and age[1]:
                print('No Lower Intercept Age')
                print('Upper Intercept Age:{:.2f} Ma.'.format(age[1]))
            elif age[1]==None and age[0]:
                print('No Upper Intercept Age')
                print('Lower Intercept Age:{:.2f} Ma.'.format(age[0]))
        else:                
            print('Lower Intercept Age:{:.2f} Ma;\nUpper Intercept Age:{:.2f} Ma.'.format(age[0],age[1]))
        '''
        print('-'*35)
        
        
        
        
        plt.figure(figsize=(10,20))
        plt.subplot(2,1,1)
        plt.title('Peak Age:'+str(round(self.PeakAge,1))+'Ma')
        plt.xlabel('Age(Ma)')
        plt.ylabel('RPD')
        plt.scatter(self.Ages,self.AgeProbSum/self.AgeProbSum.max(),c='b', marker='o')
        
        plt.axvline(self.PeakAge,color='#ff062bff',linestyle='--',label='Peak Age',linewidth=2)
        plt.axvline(AbsAgeErrLow,color='#aa2b9f',linestyle='--',label='-95%',linewidth=2)
        plt.axvline(AbsAgeErrHigh,color='green',linestyle='--',label='+95%',linewidth=2)
        #plt.text(AbsAgeErrLow,0.8, '-95%',c= 'b')
        #plt.text(AbsAgeErrHigh,0.8, '+95%',c= 'b')
        plt.legend()
        plt.grid()
        
        
        
        plt.subplot(2,1,2)
        plt.title('Peak PbC:'+str(round(self.PeakPbC,4)))
        plt.xlabel('PbC')
        plt.ylabel('RPD')
        plt.axvline(self.PeakPbC,color='#ff062bff',linestyle='--',label='Peak PbC',linewidth=2)
        plt.axvline(AbsYintErrLow,color='#aa2b9f',linestyle='--',label='-95%',linewidth=2)
        plt.axvline(AbsYintErrHigh,color='green',linestyle='--',label='+95%',linewidth=2)
        #plt.text(AbsYintErrLow,0.8, '-95%',c= 'b')
        #plt.text(AbsYintErrHigh,0.8, '+95%',c= 'b')
        plt.scatter(self.Pb76,self.YintProbSum/self.YintProbSum.max(), marker='o')
        plt.legend()
        plt.grid() 
        
        plt.savefig(self.cwd+'/'+self.filename+'_Scatter')        
        
        
        
        
        fig=plt.figure(figsize=(5,5),dpi=300)
        ax=plt.axes(projection='3d')
        x=np.array(self.Ages)
        y=np.array(self.Pb76)
        X,Y=np.meshgrid(x,y)
        
        sur=ax.plot_surface(X,Y,self.LogProb_all.T,cmap='rainbow',linewidth=0.1,alpha=0.9,cstride=1,rstride=1)
        #fig.colorbar(sur, shrink=0.5, aspect=5)
          
        #ax.contour(X,Y,self.LogProb_all.T,alpha=0.75,zdir='z', offset=0,cmap="rainbow")   #生成z方向投影
        ax.contour(X,Y,self.LogProb_all.T, alpha=0.85,zdir='x',offset=AgeEnd,cmap="rainbow")   #生成x方向投影
        ax.contour(X,Y,self.LogProb_all.T, alpha=0.85,zdir='y',offset=Pb76Start,cmap="rainbow")   #生成y方向投影


        ax.view_init(elev=30,azim=125)
        ax.set_xlabel('Age(Ma)')
        ax.set_ylabel('PbC')
        ax.set_zlabel('RPD')
        

        
        plt.savefig(self.cwd+'/'+self.filename+'_3D')
        if pltshowT:
            pass
        else:plt.show()

        #return self.xint_all,self.LogProb_all,self.Ages,self.Pb76,self.YintProbSum,self.AgeProbSum

    def findpeak(self):
        PeakAge = self.InitAge
        PeakPbC = self.InitPbC
        AgeIncr=self.AgeIncr
        PbCIncr=self.PbCIncr
        NumIter = 0
        MaxAge = False
        MaxPbC = False
        
        while AgeIncr>0.009 :            
            while MaxAge==False or MaxPbC==False:            
                NumIter = NumIter + 1
                #print('A:',NumIter,PeakAge,PeakPbC)
                RPD11=self.RPD(PeakAge,PeakPbC)
                
                
                PeakAge=PeakAge+AgeIncr
                RPD21=self.RPD(PeakAge,PeakPbC)
                
                DiffLogX=RPD21-RPD11
                if DiffLogX < 0:
                    PeakAge = PeakAge - 2 * AgeIncr
                    #print('B:',NumIter,PeakAge,PeakPbC)
                    RPD31=self.RPD(PeakAge,PeakPbC)
                    DiffLogX=RPD31-RPD11
                    if DiffLogX<0:
                        PeakAge=PeakAge + AgeIncr
                        #print('C:',NumIter,PeakAge,PeakPbC)
                        MaxAge = True
                    else:
                        MaxAge =False
                RPD11=self.RPD(PeakAge,PeakPbC)
                PeakPbC = PeakPbC + PbCIncr
                #print('D:',NumIter,PeakAge,PeakPbC)
                RPD12=self.RPD(PeakAge,PeakPbC)
                DiffLogY=RPD12-RPD11
                if DiffLogY<0:
                    PeakPbC = PeakPbC - 2 * PbCIncr
                    #print('E:',NumIter,PeakAge,PeakPbC)
                    RPD13=self.RPD(PeakAge,PeakPbC)
                    DiffLogY=RPD13-RPD11
                    if DiffLogY<0:
                        PeakPbC = PeakPbC + PbCIncr
                        #print('F:',NumIter,PeakAge,PeakPbC)
                        MaxPbC = True
                    else:
                        MaxPbC = False
            
            AgeIncr=AgeIncr/10
            PbCIncr = PbCIncr / 10
        if NumIter>300:
            print('Warning!Iterations > 300')
        return  NumIter,PeakAge,PeakPbC
         
 
    def rejectdata(self):
        Datum=[]
        NumUsed=0
        for i in range(self.result.shape[0]):
            if self.result[i,2]>2 and self.result[i,3]>2 and self.result[i,4]>2:
               NumUsed=NumUsed+1
               Datum.append([self.result[i,2]/self.Del206,self.Corr76*self.result[i,3]/self.Del207,self.Corr86*self.result[i,4]/self.Del238,self.result[i,5]/self.Del206,self.result[i,6]/self.Del207,self.result[i,7]/self.Del238])
        self.Datum=np.array(Datum)
        #print(NumUsed)
        return self.Datum
    def pdf(self,x, mean, var):
        return exp(-(x - mean) ** 2 / (2 * var ** 2)) / sqrt(2 * pi) * var
    
    def RPD(self,PeakAge, PeakPbC): 
        LineLogProb = 0
        Yint_Jy=PeakPbC
        Xint_Ix_Jy=self.XintValue(PeakAge,PeakPbC)
        for num in range(self.Datum.shape[0]):
            
            Pb206Y_Datum=self.Datum[num,0]
            Pb207Z_Datum=self.Datum[num,1]
            U238X_Datum=self.Datum[num,2]
            SigY_Datum=self.Datum[num,3]
            SigZ_Datum=self.Datum[num,4]
            SigX_Datum=self.Datum[num,5]

            
            Kst = Pb206Y_Datum * Yint_Jy - U238X_Datum * Yint_Jy / Xint_Ix_Jy - Pb207Z_Datum
            A = (SigZ_Datum/SigX_Datum)**2+(Yint_Jy/Xint_Ix_Jy)**2
            B = -2*(Yint_Jy**2)/Xint_Ix_Jy
            C = (SigZ_Datum/SigY_Datum)**2+Yint_Jy**2
            D = -2*Kst*Yint_Jy/Xint_Ix_Jy
            E = 2*Kst*Yint_Jy
                
            Disc =B**2-4*A*C
                   
            f = sqrt((Kst**2-(C*D**2+A* E**2-B*E*D)/Disc)/SigZ_Datum**2)
                     
                    #X & Y offsets from datum to tangent point on mixing plane
            Xc = (B*E-2*C*D)/Disc
            Yc = (B*D-2*A*E)/Disc
            Zc = Yint_Jy*(Yc-Xc/Xint_Ix_Jy)-Kst
                    
            #Convert to absolute coordinates
            Xc =Xc+U238X_Datum
            Yc =Yc+Pb206Y_Datum
            Zc =Zc+Pb207Z_Datum
            #Pick the nearest integer count values by rounding up or down.
            #Assumes only discrete data. Results in fractal PDD surface
            #Xc = Round(Xc, 0)
            #Yc = Round(Yc, 0)
            #Zc = Round(Zc, 0)                
            #Determine Gaussian probability of point on compressed mixing plane wrt datum point.
            Prob238X=self.pdf(Xc, U238X_Datum, SigX_Datum)
            Prob206Y=self.pdf(Yc, Pb206Y_Datum, SigY_Datum)
            Prob207Z=self.pdf(Zc, Pb207Z_Datum, SigZ_Datum)
            #print('RPD:',Prob238X,Prob206Y,Prob207Z)        
            if (Prob238X==0 or Prob206Y==0 or Prob207Z==0):
                
                pass
            else:
                LineLogProb= LineLogProb + (log(Prob238X) + log(Prob207Z) + log(Prob206Y)) / log(10)
        
        return LineLogProb

        



        
class MyApplication(QWidget):
    def __init__(self):
        super().__init__()
        self.th=Bayesian_Regression()
        self.th.signalForText.connect(self.onUpdateText)
        sys.stdout=self.th
        self.settings=qc.QSettings('Config.ini',qc.QSettings.IniFormat)
        if os.path.exists('./config.ini'):
            self.initUI()
        else:
            self.save_info()
    def closeEvent(self, event):

        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit?", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
    def save_info(self):
        self.settings.setValue('Dwell_Time/Pb206','0.06')
        self.settings.setValue('Dwell_Time/Pb207','0.06')
        self.settings.setValue('Dwell_Time/U238','0.02')
        self.settings.setValue('Dwell_Time/Others','0.01') 
        
        self.settings.setValue('Dwell_Time/No_Detector','5')
       
        self.settings.setValue('Std_biase_factor/Pb207-Pb206','1.0021')
        self.settings.setValue('Std_biase_factor/Pb206-U238','0.7726')
        self.settings.setValue('Std_biase_factor/Ablration_factor','0.91')
        
        
        
        self.settings.setValue('Signal_setting/Backgroud_Start','1')
        self.settings.setValue('Signal_setting/Backgroud_End','10')
        self.settings.setValue('Signal_setting/Sampling_Start','135')
        self.settings.setValue('Signal_setting/Sampling_End','705')
        
        self.settings.setValue('Fileter_setting/Baseline_Erro','4')
        self.settings.setValue('Fileter_setting/MinU-Pb','4')

        self.settings.setValue('Fileter_setting/MinU_counts','50')
        self.settings.setValue('Fileter_setting/MinPb_counts','25')
        self.settings.setValue('Fileter_setting/MaxMSWD','10')
        self.settings.setValue('Bayesian_setting/Age_Start','102')
        self.settings.setValue('Bayesian_setting/Age_End','120')
        self.settings.setValue('Bayesian_setting/Pbc_Start','0.795')
        self.settings.setValue('Bayesian_setting/Pbc_End','0.887')
        self.initUI()
    def save_config(self):
        filepath,type=QFileDialog.getSaveFileName(self,'Save as',"/" ,'ini(*.ini)')
        settings=qc.QSettings(filepath,qc.QSettings.IniFormat)
        
        settings.setValue('Dwell_Time/Pb206',self.dewll206le.text())
        settings.setValue('Dwell_Time/Pb207',self.dewll207le.text())
        settings.setValue('Dwell_Time/U238',self.dewll238le.text())
        settings.setValue('Dwell_Time/Others',self.dewll232le.text())
        
        settings.setValue('Dwell_Time/No_Detector',self.Detectors_Noe.text())
        
        settings.setValue('Std_biase_factor/Pb207-Pb206',self.Stdbias76e.text())
        settings.setValue('Std_biase_factor/Pb206-U238',self.Stdbias68e.text())
        settings.setValue('Std_biase_factor/Ablration_factor',self.Ablre.text())
        settings.setValue('Signal_setting/Backgroud_Start',self.Blankfe.text())
        settings.setValue('Signal_setting/Backgroud_End',self.Blankte.text())
        settings.setValue('Signal_setting/Sampling_Start',self.Signalfe.text())
        settings.setValue('Signal_setting/Sampling_End',self.Signalte.text())
        
        settings.setValue('Fileter_setting/Baseline_Erro',self.bsere.text())
        settings.setValue('Fileter_setting/MinU-Pb',self.minU_Pble.text())

        settings.setValue('Fileter_setting/MinU_counts',self.MinUle.text())
        settings.setValue('Fileter_setting/MinPb_counts',self.MinUCTle.text())
        settings.setValue('Fileter_setting/MaxMSWD',self.Maxmswdle.text())
        settings.setValue('Bayesian_setting/Age_Start',self.Agesle.text())
        settings.setValue('Bayesian_setting/Age_End',self.Ageele.text())
        settings.setValue('Bayesian_setting/Pbc_Start',self.Pbcsle.text())
        settings.setValue('Bayesian_setting/Pbc_End',self.Pbcele.text())
        print('Save:{}'.format(filepath))
    def load_config(self):        
        filepath,type=QFileDialog.getOpenFileName(self,'Load',"." ,'ini(*.ini)')
        self.settings=qc.QSettings(filepath,qc.QSettings.IniFormat)
        self.dewll206le.setText(self.settings.value('Dwell_Time/Pb206'))
        self.dewll207le.setText(self.settings.value('Dwell_Time/Pb207'))
        self.dewll232le.setText(self.settings.value('Dwell_Time/Others'))
        self.dewll238le.setText(self.settings.value('Dwell_Time/U238'))
        self.Detectors_Noe.setText(self.settings.value('Dwell_Time/No_Detector'))
        self.Stdbias76e.setText(self.settings.value('Std_biase_factor/Pb207-Pb206'))
        self.Stdbias68e.setText(self.settings.value('Std_biase_factor/Pb206-U238'))
        self.Ablre.setText(self.settings.value('Std_biase_factor/Ablration_factor'))
        
        self.Blankfe.setText(self.settings.value('Signal_setting/Backgroud_Start'))
        self.Blankte.setText(self.settings.value('Signal_setting/Backgroud_End'))
        self.Signalfe.setText(self.settings.value('Signal_setting/Sampling_Start'))
        self.Signalte.setText(self.settings.value('Signal_setting/Sampling_End'))
        self.bsere.setText(self.settings.value('Fileter_setting/Baseline_Erro'))
        self.MinUle.setText(self.settings.value('Fileter_setting/MinU_counts'))
        self.MinUCTle.setText(self.settings.value('Fileter_setting/MinPb_counts'))
        self.minU_Pble.setText(self.settings.value('Fileter_setting/MinU-Pb'))
        self.Maxmswdle.setText(self.settings.value('Fileter_setting/MaxMSWD'))
        self.Agesle.setText(self.settings.value('Bayesian_setting/Age_Start'))
        self.Ageele.setText(self.settings.value('Bayesian_setting/Age_End'))
        self.Pbcsle.setText(self.settings.value('Bayesian_setting/Pbc_Start'))
        self.Pbcele.setText(self.settings.value('Bayesian_setting/Pbc_End'))
        
        print('load:{}'.format(filepath))
        
        
    def onUpdateText(self,text):
        cursor=self.process.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.process.setTextCursor(cursor)
        self.process.ensureCursorVisible()
        
    def initUI(self):
        self.setWindowTitle('Brama v2.0')
        #self.setFixedSize(450,350)
        self.resize(600,400)

        self.data_Num=None
        
        self.process=QTextEdit(self,readOnly=True)

        #self.processbar=QTextBrowser(self)
        #mainLayout.addWidget(self.processbar,10,5,1,1)
        
        
        
        
        
        importdatalabel=QLabel('MS file',self)
        loadData1 = QPushButton('Browse',self)
        loadData1.setIcon(QIcon('./images/load.jpg'))        
        self.dataEntry1 = QLineEdit(self)
        loadData1.clicked.connect(lambda: self.get_file(self.dataEntry1))
        
        btnLOAD=QPushButton('Plot',self)
        self.btnLOAD=btnLOAD.clicked.connect(self.plotshow)
        
        btn_Config_Save=QPushButton('Save Settings',self)
        self.btn_Config_Save=btn_Config_Save.clicked.connect(self.save_config)
         
        btn_Config_Load=QPushButton('Load Settings',self)
        self.btn_Config_Load=btn_Config_Load.clicked.connect(self.load_config)

        self.STDbia=QCheckBox('Std Corr?',self)
        self.STDbia.setChecked(True)
        
        
        Stdagel=QLabel('Age of RM (Ma)',self) 
        
        self.Stdagee=QLineEdit(self)
        self.Stdagee.setText('100')
        btn207Pb=QPushButton('207Pb',self)
        self.btn207Pb=btn207Pb.clicked.connect(self.Std207Pb)
       
        
        
        

        
        
        
        
        

        rbAgilent1=QRadioButton('Agilent1',self)
        
        rbAgilent1.toggled.connect(lambda: self.butnstate(rbAgilent1))
        rbAgilent2=QRadioButton('Agilent2',self)
        rbAgilent2.toggled.connect(lambda: self.butnstate(rbAgilent2))
        rbAgilent2.setChecked(True)
        Thermo=QRadioButton('Thermo',self)
        Thermo.toggled.connect(lambda: self.butnstate(Thermo))
        
        
        
        

        Reject_levell=QLabel('Multiples of SD',self)
        self.Reject_levele=QLineEdit(self)
        self.Reject_levele.setText('2')
        
        btnDel_spikes=QPushButton('Delete Spikes',self)
        self.btnDel_spikes=btnDel_spikes.clicked.connect(self.Del_spikes)
        Reject_label=QLabel('Peak Rejection',self)
        
 
        log_path_label=QLabel('List file:',self)
        LogData = QPushButton('Browse',self)
        LogData.setIcon(QIcon('./images/load.jpg')) 
        self.dataEntry2 = QLineEdit(self)
        LogData.clicked.connect(lambda: self.get_file_xls(self.dataEntry2))
        
        Output_path_label=QLabel('Export:',self)
        OutputData = QPushButton('Browse',self)
        OutputData.setIcon(QIcon('./images/load.jpg'))        
        self.dataEntry3 = QLineEdit(self)
        OutputData.clicked.connect(lambda: self.get_dirctory(self.dataEntry3))

  

        
        dwelltime=QLabel('Dwell time(s)',self) 
        dwell206l=QLabel('&206Pb',self)
        self.dewll206le=QLineEdit(self)
        
        
       
        
        Detectors_No=QLabel('Detectors (N)',self)
        self.Detectors_Noe=QLineEdit(self)
        
        self.Detectors_Noe.setText(self.settings.value('Dwell_Time/No_Detector'))
        
        
        
        self.dewll206le.setText(self.settings.value('Dwell_Time/Pb206'))
       
        dwell206l.setBuddy(self.dewll206le)
        dwell207l=QLabel('&207Pb',self)
        self.dewll207le=QLineEdit(self)
        self.dewll207le.setText(self.settings.value('Dwell_Time/Pb207'))
        dwell207l.setBuddy(self.dewll207le)
        dwell232l=QLabel('Others',self)
        self.dewll232le=QLineEdit(self)
        dwell232l.setBuddy(self.dewll232le)
        self.dewll232le.setText(self.settings.value('Dwell_Time/Others'))
        dwell238l=QLabel('&238U',self)
        self.dewll238le=QLineEdit(self)
        self.dewll238le.setText(self.settings.value('Dwell_Time/U238'))
        dwell238l.setBuddy(self.dewll238le)

        
        
        self.CommPbrb1=QCheckBox('ComPb Filter?',self)

        self.CommPbrb1.setChecked(False)
        
    
        
        self.Interpolrb1=QCheckBox('Interpolate?',self)
        
        self.Interpolrb1.setChecked(True)
        
        self.Delspikes=QCheckBox('Del Spikes?',self)

        self.Delspikes.setChecked(False)
        
        
        
        
        

                
        CorrectionFactorl=QLabel('Correction Factor',self)        
        Stdbias76=QLabel('&207Pb/206Pb',self)
        self.Stdbias76e=QLineEdit(self)
        self.Stdbias76e.setText(self.settings.value('Std_biase_factor/Pb207-Pb206'))
        Stdbias76.setBuddy(self.Stdbias76e)        
        Stdbias68=QLabel('&206Pb/238U',self)
        self.Stdbias68e=QLineEdit(self)
        self.Stdbias68e.setText(self.settings.value('Std_biase_factor/Pb206-U238'))
        Stdbias68.setBuddy(self.Stdbias68e)
        Ablr=QLabel('&Ablation Cor.',self)
        #self.Ablre=QLineEdit(self,readOnly=True)
        self.Ablre=QLineEdit(self,readOnly=False)
        self.Ablre.setText(self.settings.value('Std_biase_factor/Ablration_factor'))
        Ablr.setBuddy(self.Ablre)
        
        
        
        
        Basestart=QLabel('Blank and Signal Rows',self)
        
        
        Blankf=QLabel('&Blank Start',self)
        self.Blankfe=QLineEdit(self)
        self.Blankfe.setText(self.settings.value('Signal_setting/Backgroud_Start'))
        Blankf.setBuddy(self.Blankfe)
        Blankt=QLabel('&Blank End',self)
        self.Blankte=QLineEdit(self)
        self.Blankte.setText(self.settings.value('Signal_setting/Backgroud_End'))
        Blankt.setBuddy(self.Blankte)
        Signalf=QLabel('&Signal Start',self)
        self.Signalfe=QLineEdit(self)
        self.Signalfe.setText(self.settings.value('Signal_setting/Sampling_Start'))
        Signalf.setBuddy(self.Signalfe)
        Signalt=QLabel('&Signal End',self)
        self.Signalte=QLineEdit(self)

        self.Signalte.setText(self.settings.value('Signal_setting/Sampling_End'))
        Signalt.setBuddy(self.Signalte)
        
        
       
        
        Filterlabel=QLabel('Filter Setting',self)        
        bser=QLabel('&Baseline Err',self)
        self.bsere=QLineEdit(self)
        self.bsere.setText(self.settings.value('Fileter_setting/Baseline_Erro'))
        bser.setBuddy(self.bsere)
        MinUl=QLabel('&U(cps)>',self)
        self.MinUle=QLineEdit(self)
        self.MinUle.setText(self.settings.value('Fileter_setting/MinU_counts'))
        MinUl.setBuddy(self.MinUle)
        MinCTl=QLabel('&U-Pb(counts)>',self)
        self.MinUCTle=QLineEdit(self)
        self.MinUCTle.setText(self.settings.value('Fileter_setting/MinPb_counts'))
        MinCTl.setBuddy(self.MinUCTle)        
        minU_Pbl=QLabel('&U/Pb >',self)
        self.minU_Pble=QLineEdit(self)
        self.minU_Pble.setText(self.settings.value('Fileter_setting/MinU-Pb'))
        minU_Pbl.setBuddy(self.minU_Pble)
        

        
        Agesl=QLabel('&Age Start (Ma)',self)
        self.Agesle=QLineEdit(self)
        self.Agesle.setText(self.settings.value('Bayesian_setting/Age_Start'))
        Agesl.setBuddy(self.Agesle)
        Ageel=QLabel('&Age End (Ma)',self)
        self.Ageele=QLineEdit(self)
        self.Ageele.setText(self.settings.value('Bayesian_setting/Age_End'))
        Ageel.setBuddy(self.Ageele)
        Agenl=QLabel('&Age Num.',self)
        self.Agenle=QLineEdit(self)
        self.Agenle.setText('50')
        Agenl.setBuddy(self.Agenle)
        
        
        Pbcsl=QLabel('&Pbc Start',self)
        self.Pbcsle=QLineEdit(self)
        self.Pbcsle.setText(self.settings.value('Bayesian_setting/Pbc_Start'))
        Pbcsl.setBuddy(self.Pbcsle)
        Pbcel=QLabel('&Pbc End',self)
        self.Pbcele=QLineEdit(self) 
        self.Pbcele.setText(self.settings.value('Bayesian_setting/Pbc_End'))
        Pbcel.setBuddy(self.Pbcele)
        Pbcnl=QLabel('&Pbc Num.',self)
        self.Pbcnle=QLineEdit(self)
        self.Pbcnle.setText('50')
        Pbcnl.setBuddy(self.Pbcnle)
        
        Maxmswdl=QLabel('&Max MSWD.',self)
        self.Maxmswdle=QLineEdit(self)


        
        self.Maxmswdle.setText(self.settings.value('Fileter_setting/MaxMSWD'))
        Maxmswdl.setBuddy(self.Maxmswdle)  

        
        btnCONTS=QPushButton('U-Pb isotope calculation',self)
 
        btnAutocentre=QPushButton('Auto Centre',self)
        
        self.btnAutocentre=btnAutocentre.clicked.connect(self.Autocentre)
        self.btnCONTS=btnCONTS.clicked.connect(self.Cal_Istope)
        btnBAYESIAN=QPushButton('Bayesian Regression',self)
 
        self.btnBAYESIAN=btnBAYESIAN.clicked.connect(self.regression)
        
        Input_path_label=QLabel('Input Dir:',self)
        InputData = QPushButton('Browse',self)
        InputData.setIcon(QIcon('./images/load.jpg'))        
        self.dataEntry4 = QLineEdit(self)

        InputData.clicked.connect(lambda: self.get_dirctory(self.dataEntry4))
        
        Batch_process_label=QLabel('Batch Process',self)
        
        
        btnBatch_process=QPushButton('Batch Process',self)
        self.btnBatch_process=btnBatch_process.clicked.connect(self.btnBatch_process)
        self.BayesianM=QCheckBox('Bayesian Method?',self)

        self.BayesianM.setChecked(False)
        
        
        btnMerged_data_plot=QPushButton('Merge Result Plotting',self)
        btnMerged_data_plot.clicked.connect(self.Merged_data_plot)
        
        mainLayout=QGridLayout(self)
        
        #import setting
        mainLayout.addWidget(importdatalabel,0,0)   
        mainLayout.addWidget(loadData1,0,1) 
        mainLayout.addWidget(self.dataEntry1,0,2,1,4)
        
        #Instrument
        mainLayout.addWidget(rbAgilent1,0,6)  
        mainLayout.addWidget(rbAgilent2,0,7)  
        mainLayout.addWidget(Thermo,0,8)
        
        #listfile setting
        mainLayout.addWidget(log_path_label,1,0)
        mainLayout.addWidget(LogData,1,1)
        mainLayout.addWidget(self.dataEntry2,1,2,1,4)
        
        #Bia fraction fractor 
        mainLayout.addWidget(self.STDbia, 1, 6,1,1)
        mainLayout.addWidget(Stdagel, 1, 7) 
        mainLayout.addWidget(self.Stdagee, 1,8,1,1) 
        mainLayout.addWidget(btn207Pb, 1, 9,1,1)
        
        #Working Environment Directory 
        mainLayout.addWidget(Output_path_label,2,0)        
        mainLayout.addWidget(OutputData,2,1) 
        mainLayout.addWidget(self.dataEntry3,2,2,1,4)
        
        #Save and load settings
        mainLayout.addWidget(btn_Config_Save,2,6,1,1) 
        mainLayout.addWidget(btn_Config_Load,2,7,1,1) 
        
        #Signal diagram
        mainLayout.addWidget(btnLOAD, 2, 10,1,1)
        
        #Dwell time settings
        mainLayout.addWidget(dwelltime,3,0,1,2)
        mainLayout.addWidget(dwell206l,4, 0)
        mainLayout.addWidget(self.dewll206le, 4, 1)
        mainLayout.addWidget(dwell207l, 5, 0)
        mainLayout.addWidget(self.dewll207le, 5, 1) 
        mainLayout.addWidget(self.dewll232le, 6, 1)       
        mainLayout.addWidget(dwell238l, 6, 0)
        
        
        mainLayout.addWidget(dwell232l, 7, 0)
        mainLayout.addWidget(self.dewll238le, 7, 1) 
        
        mainLayout.addWidget(Detectors_No, 8, 0)
        mainLayout.addWidget(self.Detectors_Noe, 8, 1)
        
    
        
        
        #Std bias  factors
        mainLayout.addWidget(CorrectionFactorl,3,2,1,2)         
        mainLayout.addWidget(Stdbias76,4,2)
        mainLayout.addWidget(self.Stdbias76e,4,3)
        mainLayout.addWidget(Stdbias68,5,2)
        mainLayout.addWidget(self.Stdbias68e,5,3)
        mainLayout.addWidget(Ablr,6,2)
        mainLayout.addWidget(self.Ablre,6,3)
        
        # Blank and Signal Rows 
        mainLayout.addWidget(Basestart,3,4,1,2)
        mainLayout.addWidget(Blankf,4, 4)
        mainLayout.addWidget(self.Blankfe, 4, 5)
        mainLayout.addWidget(Blankt, 5, 4)
        mainLayout.addWidget(self.Blankte, 5, 5)        
        mainLayout.addWidget(Signalf, 6, 4)
        mainLayout.addWidget(self.Signalfe, 6, 5)       
        mainLayout.addWidget(Signalt, 7, 4)
        mainLayout.addWidget(self.Signalte, 7, 5) 
        
        
        #Filter Setting
        mainLayout.addWidget(Filterlabel,3,6,1,2)
        mainLayout.addWidget(bser,4, 6)
        mainLayout.addWidget(self.bsere, 4, 7)
        mainLayout.addWidget(MinUl, 5, 6)
        mainLayout.addWidget(self.MinUle, 5, 7)        
        mainLayout.addWidget(MinCTl, 6, 6)
        mainLayout.addWidget(self.MinUCTle, 6, 7)       
        mainLayout.addWidget(minU_Pbl, 7, 6)
        mainLayout.addWidget(self.minU_Pble, 7, 7) 
        
        mainLayout.addWidget(self.CommPbrb1, 4,8,2,2)        
        mainLayout.addWidget(self.Interpolrb1, 6,8,2,2) 
        
        #Peak Rejection settings
        
        mainLayout.addWidget(Reject_label,3,10)
        mainLayout.addWidget(self.Delspikes,4,10)
        mainLayout.addWidget(Reject_levell,5,10)
        mainLayout.addWidget(self.Reject_levele,6,10)
        mainLayout.addWidget(btnDel_spikes,7,10)
        
        
        #U-Pb isotope Calculation Button
        mainLayout.addWidget(btnCONTS, 10, 0,1,6)
        
        
        #Auto Centre Bayesian Button
        mainLayout.addWidget(btnAutocentre, 11, 1,1,4)
        
        #Bayesian Regression Settings
        mainLayout.addWidget(Agesl,12, 0)
        mainLayout.addWidget(self.Agesle, 12, 1)
        mainLayout.addWidget(Ageel, 13, 0)
        mainLayout.addWidget(self.Ageele, 13, 1) 
        mainLayout.addWidget(Agenl, 14, 0)
        mainLayout.addWidget(self.Agenle, 14, 1)        
                
        mainLayout.addWidget(Pbcsl,12, 2)
        mainLayout.addWidget(self.Pbcsle, 12, 3)
        mainLayout.addWidget(Pbcel, 13, 2)
        mainLayout.addWidget(self.Pbcele, 13, 3)        
        mainLayout.addWidget(Pbcnl, 14, 2)
        mainLayout.addWidget(self.Pbcnle,14, 3)
        
        mainLayout.addWidget(Maxmswdl, 12, 4)
        mainLayout.addWidget(self.Maxmswdle, 12, 5)
        
        
        #Bayesian Regression Button
        mainLayout.addWidget(btnBAYESIAN,15, 0,1,6)
        
        #Info Display Window
        mainLayout.addWidget(self.process,8,6,9,5)
        
 
             
        
        
        
        #Batch Process and Mapping
        mainLayout.addWidget(Batch_process_label,16,0,1,2)
        mainLayout.addWidget(self.BayesianM,16,2,1,3)
        
        mainLayout.addWidget(Input_path_label,17,0)
        mainLayout.addWidget(InputData,17,1)
        mainLayout.addWidget(self.dataEntry4,17,2,1,4)
        mainLayout.addWidget(btnBatch_process,17,6,1,1)        
        mainLayout.addWidget(btnMerged_data_plot,17,8,1,2)
        
        #set ContentsMargins 
        mainLayout.setContentsMargins(15,15,15,15)
        mainLayout.setSpacing(5)
        #mainLayout.setRowStretch(0, 1)
        #mainLayout.setRowStretch(1, 30)
        #mainLayout.setColumnStretch(0, 1)
        #mainLayout.setColumnStretch(1, 30)
        
        
    def Std207Pb(self):
        if self.STDbia.isChecked():           
            print('Deducting for common Pb in the Standard (207Pb Method).')
            print('Reference values for standard calculated with S&K (1975)')
        else:
            print('Reference values for standard calculated with S&K (1975)')
        
        try:
            
            rawdata=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
            data_all=np.loadtxt(rawdata,dtype=str,delimiter=',',skiprows=self.instrument,comments='     ') 
            #print(data_all)
            col_data_name=data_all[0,:].tolist()
            print(col_data_name)
            col_name=[]
            b0=int(self.Blankfe.text())
            b1=int(self.Blankte.text())
            s0=int(self.Signalfe.text())
            s1=int(self.Signalte.text())
            age=float(self.Stdagee.text())
            for i in range(len(col_data_name)):
                col_name.append(''.join(list(filter(str.isdigit, col_data_name[i]))))
                if self.instrument==13:
                    Time=data_all[2:,0].astype(float)
                else:
                    Time=data_all[1:,0].astype(float)
            Num=Time.shape[0]
            print(Num,b0,b1,s0,s1)
            if self.instrument==13: 
                Hg202=np.zeros(Num)   
                Hg204=np.zeros(Num) 
                Pb208=np.zeros(Num)
                
                Pb206=data_all[2:,col_name.index('206')].astype(float)
                Pb207=data_all[2:,col_name.index('207')].astype(float)                  
                U238=data_all[2:,col_name.index('238')].astype(float)
                try:
                    Th232=data_all[1:,col_name.index('232')].astype(float)
                except:
                    Th232=np.zeros(Num)
                
    
                Hg202_b=Hg202[b0:b1]
                Hg204_b=Hg204[b0:b1]
                Pb206_b=Pb206[b0:b1]
                Pb207_b=Pb207[b0:b1]
                Pb208_b=Pb208[b0:b1]
                U238_b=U238[b0:b1]
                Th232_b=Th232[b0:b1]
                
                Hg202_s=Hg202[s0:s1]
                Hg204_s=Hg204[s0:s1]
                Pb206_s=Pb206[s0:s1]
                Pb207_s=Pb207[s0:s1]
                Pb208_s=Pb208[s0:s1]
                U238_s=U238[s0:s1]
                Th232_s=Th232[s0:s1]
            else:
                Hg202=np.zeros(Num)   
                Hg204=np.zeros(Num)
                Pb208=np.zeros(Num)
                
                
                Pb206=data_all[1:,col_name.index('206')].astype(float)
                Pb207=data_all[1:,col_name.index('207')].astype(float)                  
                U238=data_all[1:,col_name.index('238')].astype(float)
                try:
                    Th232=data_all[1:,col_name.index('232')].astype(float)
                except:
                    Th232=np.zeros(Num)
                
                Hg202_b=Hg202[b0:b1]
                Hg204_b=Hg204[b0:b1]
                Pb206_b=Pb206[b0:b1]
                Pb207_b=Pb207[b0:b1]
                Pb208_b=Pb208[b0:b1]
                U238_b=U238[b0:b1]
                Th232_b=Th232[b0:b1]
                
                Hg202_s=Hg202[s0:s1]
                Hg204_s=Hg204[s0:s1]
                Pb206_s=Pb206[s0:s1]
                Pb207_s=Pb207[s0:s1]
                Pb208_s=Pb208[s0:s1]
                U238_s=U238[s0:s1]
                Th232_s=Th232[s0:s1]
    
            
            
            result_outstd=[]
            
            result_std=self.Std_method(Hg202_b,Hg204_b,Pb206_b,Pb207_b,Pb208_b,Th232_b,U238_b,Hg202_s,Hg204_s,Pb206_s,Pb207_s,Pb208_s,Th232_s,U238_s,b0,b1,s0,s1,age)
            result_outstd.append(result_std)
            
            with open("result_all.csv", "w", newline='') as s:
                writer=csv.writer(s)
                writer.writerow(['Age','207Pb/206Pb','2s','206Pb/238U','2s','207Pb/235Uc','2s',\
                                 '208Pb/232Th','2s','208Pb/206Pb','2s','232Th/206Pb','2s','208Pb/204Pb','2s','Trace element ','U(cps)','Th(cps)',\
                                 'Pb208(cps)','Pb207(cps)','Pb206(cps)'])
                for r in result_outstd:
                    writer.writerow(r)
                    
            R8,Q8,S8,P382,a,b,c=Cal_age(age)
            print(result_outstd)
            f207_206=result_outstd[0][1]/P382
            f206_238=result_outstd[0][3]/a
            f207_235=result_outstd[0][5]/b
            f208_232=result_outstd[0][7]/c
            print('Fractionation factors:\n207Pb/206Pb:{:.4f}.\n206Pb/238U :{:.4f}.\n207Pb/235U :{:.4f}.\n208Pb/232Th:{:.4f}.\n'.format(f207_206,f206_238,f207_235,f208_232))
        except Exception as e:
          print(str(e))
        
    def Std_method(self,Hg202_b,Hg204_b,Pb206_b,Pb207_b,Pb208_b,Th232_b,U238_b,Hg202_s,Hg204_s,Pb206_s,Pb207_s,Pb208_s,Th232_s,U238_s,b0,b1,s0,s1,age):
        
        R8,Q8,S8,P382,a,b,c=Cal_age(age)
        start_b=b0
        end_b=b1
        start_s=s0
        end_s=s1
    
        def Std_process(Hg202_b,Hg204_b,Pb206_b,Pb207_b,Pb208_b,Th232_b,U238_b,Hg202_s,Hg204_s,\
                        Pb206_s,Pb207_s,Pb208_s,Th232_s,U238_s,b0,b1,s0,s1,age,N383):
            R8,Q8,S8,P382,a,b,c=Cal_age(age)
            start_b=b0
            end_b=b1
            start_s=s0
            end_s=s1
            # 样品计算
            Hg204_Hg202_s=Hg204_s/Hg202_s
            Hg204cal_s=(np.average(Hg202_s)-np.average(Hg202_s))*0.229883
            
            m204_s=np.where(Hg204cal_s<(np.std(Hg204_s)/sqrt(end_s-start_s)),0,Hg204cal_s)
            
            Pb204_s=np.average(Hg204_s)-np.average(Hg204_b)
            
            U238_T=np.average(U238_s)-np.average(U238_b)
            Th232_T=np.average(Th232_s)-np.average(Th232_b)
            Pb208_T=(np.average(Pb208_s)-np.average(Pb208_b))
            Pb207_T=(np.average(Pb207_s)-np.average(Pb207_b))
            Pb206_T=(np.average(Pb206_s)-np.average(Pb206_b))
            
            Pb207_Pb206_s=(Pb207_s-np.average(Pb207_b)-N383*R8)/(Pb206_s-np.average(Pb206_b)-N383*Q8)
            
            Pb206_U238_s=(Pb206_s-np.average(Pb206_b)-N383*Q8)/(U238_s-np.average(U238_b))
                
            Pb208_Th232_s=(Pb208_s-np.average(Pb208_b)-N383*Q8)/(Th232_s-np.average(Th232_b))
            
            Pb208_Pb206_s=(Pb208_s-np.average(Pb208_b)-N383*Q8)/(Pb206_s-np.average(Pb206_b))
            
            Th232_Pb206_s=Th232_s/Pb206_s
            
            Pb208_Pb204_s=(Pb208_s-np.average(Pb208_b))/Pb204_s
            
            
            #过滤数组
            
            
            x=np.arange(len(Pb206_s))
            Pb207_Pb206_f=np.where(np.abs(Pb207_Pb206_s-np.average(Pb207_Pb206_s))<2*np.std(Pb207_Pb206_s),Pb207_Pb206_s,nan)
            
            Pb206_U238_f=np.where(np.abs(Pb206_U238_s-np.average(Pb206_U238_s))<2*np.std(Pb206_U238_s),Pb206_U238_s,nan)
             
            
            
               
            Pb208_Th232_f=np.where(np.abs(Pb208_Th232_s-np.average(Pb208_Th232_s))<2*np.std(Pb208_Th232_s),Pb208_Th232_s,nan)
            
            Pb208_Pb206_f=np.where(np.abs(Pb208_Pb206_s-np.average(Pb208_Pb206_s))<2*np.std(Pb208_Pb206_s),Pb208_Pb206_s,nan)
            
            
            
            Th232_Pb206_f=np.where(np.abs(Th232_Pb206_s-np.average(Th232_Pb206_s))<2*np.std(Th232_Pb206_s),Th232_Pb206_s,nan)
            
            Pb208_Pb204_f=np.where(np.abs(Pb208_Pb204_s-np.average(Pb208_Pb204_s))<2*np.std(Pb208_Pb204_s),Pb208_Pb204_s,nan)
            
            Pb207_U235_c=np.nanmean(Pb207_Pb206_f)*np.nanmean(Pb206_U238_f)*137.88
            
    
            corr_Pb207_206=np.nanmean(Pb207_Pb206_f)
            
    
            result_std=[age,np.nanmean(Pb207_Pb206_f),\
                        2*np.nanstd(Pb207_U235_c/Pb206_U238_f/137.88)/sqrt(len(Pb207_U235_c/Pb206_U238_f/137.88)-np.count_nonzero(np.isnan(Pb207_U235_c/Pb206_U238_f/137.88))),
                        np.nanmean(Pb206_U238_f),\
                        2*np.nanstd(Pb206_U238_f)/sqrt(len(Pb206_U238_f)-np.count_nonzero(np.isnan(Pb206_U238_f))),
                        Pb207_U235_c,\
                        sqrt(pow(Pb207_U235_c*(((2*np.nanstd(Pb206_U238_f)/sqrt(len(Pb206_U238_f)-np.count_nonzero(np.isnan(Pb206_U238_f)))))/np.nanmean(Pb206_U238_f))*100/100,2)+pow(Pb207_U235_c*(((2*np.nanstd(Pb207_Pb206_f)/sqrt(len(Pb207_Pb206_f)-np.count_nonzero(np.isnan(Pb207_Pb206_f)))))/np.nanmean(Pb207_Pb206_f))*100/100,2)),
                        np.nanmean(Pb208_Th232_f),\
                        2*np.nanstd(Pb208_Th232_f)/sqrt(len(Pb208_Th232_f)-np.count_nonzero(np.isnan(Pb208_Th232_f))),
                        np.nanmean(Pb208_Pb206_f),\
                        2*np.nanstd(Pb208_Pb206_f)/sqrt(len(Pb208_Pb206_f)-np.count_nonzero(np.isnan(Pb208_Pb206_f))),
                        np.nanmean(Th232_Pb206_f),\
                        2*np.nanstd(Th232_Pb206_f)/sqrt(len(Th232_Pb206_f)-np.count_nonzero(np.isnan(Th232_Pb206_f))),
                        np.nanmean(Pb208_Pb204_f),\
                        2*np.nanstd(Pb208_Pb204_f)/sqrt(len(Pb208_Pb204_f)-np.count_nonzero(np.isnan(Pb208_Pb204_f)))
                        ,'------',U238_T,Th232_T,Pb208_T,Pb207_T,Pb206_T]
            
            return result_std
    
        Pb207_Pb206_s=np.ones(Pb206_s.shape)*P382 
        if self.STDbia.isChecked():
            N383=(Pb207_Pb206_s*Pb206_s-np.average(Pb206_b)*Pb207_Pb206_s-Pb207_s+np.average(Pb207_b))/(Q8*Pb207_Pb206_s-R8)
            
        else:
            N383=0
            print(N383) 
        result_std=Std_process(Hg202_b,Hg204_b,Pb206_b,Pb207_b,Pb208_b,Th232_b,U238_b,Hg202_s,Hg204_s,\
                           Pb206_s,Pb207_s,Pb208_s,Th232_s,U238_s,b0,b1,s0,s1,age,N383)          
        return  result_std        
        
    def pikes(self,A,multi):
        B=[]
        B.append(A[0])
        for i in range(A.size-1):
            #print(A[i])
            if A[i]>(((A[i-1]+A[i+1])*0.5)*multi) :
                B.append((A[i-1]+A[i+1])/2)
            else:
                B.append(A[i])
        B.append(A[-1])
        F=np.array(B) 
        return F      

    def Del_spikes(self):
        multi=eval(self.Reject_levele.text())
        print('Delete Spikes:',multi,'times of standard deviations were rejected.')
        rawdata=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
        
        data_all=np.loadtxt(rawdata,dtype=str,delimiter=',',skiprows=self.instrument,comments='     ') 
        
        
    #data_all=np.loadtxt(rawdata,dtype=str,delimiter=',',skiprows=2,comments='     ')
    
        
        col_data_name=data_all[0,:].tolist()
    #print(col_data_name)
        col_name=[]
        for i in range(len(col_data_name)):
            col_name.append(''.join(list(filter(str.isdigit, col_data_name[i]))))
            if self.instrument==13:
                Time=data_all[2:,0].astype(float)
            else:
                Time=data_all[1:,0].astype(float)
        if self.instrument==13:                
            Pb206=data_all[2:,col_name.index('206')].astype(float)
            Pb207=data_all[2:,col_name.index('207')].astype(float)                  
            U238=data_all[2:,col_name.index('238')].astype(float)
            #Th232=data_all[2:,col_name.index('232')]
            Pb206=self.pikes(Pb206,multi)
            Pb207=self.pikes(Pb207,multi)
            U238=self.pikes(U238,multi)
            #Th232=self.pikes(Th232,multi)
        else:
            Pb206=data_all[1:,col_name.index('206')].astype(float)
            Pb207=data_all[1:,col_name.index('207')].astype(float)                  
            U238=data_all[1:,col_name.index('238')].astype(float)
            Th232=data_all[1:,col_name.index('232')]   
            Pb206=self.pikes(Pb206,multi)
            Pb207=self.pikes(Pb207,multi)
            U238=self.pikes(U238,multi)
            #Th232=self.pikes(Th232,multi)
            
        plt.title('Rows vs. Signal(cps)')

        plt.plot(Pb206,label='206Pb')
        plt.plot(Pb207,label='207Pb')

        plt.plot(U238,label='238U')
        ax=plt.gca()
        ax.set_yscale('log')
        plt.legend()
        plt.grid(axis='x')
        plt.tight_layout()
        plt.xlabel("Rows")
        plt.ylabel("Counts (cps)")
    
    
        plt.show()
        if self.dataEntry1.text():
            if self.instrument==2 or self.instrument==3:
                traspath1=os.path.dirname(self.dataEntry1.text())
                traspath=os.path.dirname(traspath1)
            elif self.instrument==13:
                traspath=os.path.dirname(self.dataEntry1.text())
           
            self.dataEntry4.setText(traspath)
        
        
   
            
    def Autocentre(self):
        try:
            print('Auto Centre Regression!')
            print('{:.2f}{:.2f}{:.4f}{:.4f}'.format(self.Agestart,self.Agesig,self.Pbstart,self.Pbsig))
            print('-'*35)
            logging.info('Regression Start.')
            #print("\r", end="")
            #print(" {}%: ".format(numi), "█" * (numi // 2), end="",flush=True)
        
            AgeStart=self.Agestart-2*self.Agesig
            AgeEnd=self.Agestart+2*self.Agesig
            DelAge=(AgeEnd-AgeStart)/int(self.Agenle.text())      
           
            Pb76Start=self.Pbstart-2*self.Pbsig
            Pb76End=self.Pbstart+2*self.Pbsig
            DelPb76=(Pb76End-Pb76Start)/int(self.Pbcnle.text())
            #print(AgeStart,AgeEnd)
            
            self.Agesle.setText(str(format(AgeStart,'.2f')))
            
            self.Ageele.setText(str(format(AgeEnd,'.2f')))
            
            
            self.Pbcsle.setText(str(format(Pb76Start,'.4f')))
            
            self.Pbcele.setText(str(format(Pb76End,'.4f')))
            




            logging.info('Age Start:%s',AgeStart)
            logging.info('Age End:%s',AgeEnd)
            logging.info('Age Num.:%s',int(self.Agenle.text()) )
            
            logging.info('Pbc76 Start:%s',Pb76Start)
            logging.info('Pbc76 End:%s',Pb76End)
            logging.info('Pbc76 Num.:%s',int(self.Pbcnle.text()) )
            
            self.a.regress(Pb76Start,Pb76End,AgeStart,AgeEnd,DelPb76,DelAge,pltshowT=False)
            

            
            end_time=time.time()
            print()
            with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
                print('#Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))),file=f)
                
                print('*'*50,file=f)
                
            print('Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))))  
            #print('Duration of processing:{:.2f} s'.format(end_time-self.start_time))  
            print('*'*30)
            self.settings.setValue('Bayesian_setting/Age_Start',self.Agesle.text())
            self.settings.setValue('Bayesian_setting/Age_End',self.Ageele.text())
            self.settings.setValue('Bayesian_setting/Pbc_Start',self.Pbcsle.text())
            self.settings.setValue('Bayesian_setting/Pbc_End',self.Pbcele.text())
        
        
            msg_box_r=QMessageBox(QMessageBox.Information,'Information','Regression successful!')
            msg_box_r.exec_()
            logging.info('Regression Compeleted!')
        except Exception as e:
            print(str(e))
        
    def btnBatch_process(self):
        try:
            BayesinC=self.BayesianM.isChecked()
            for name in self.Sampleslist1.keys():
                self.dir_path=QtCore.QDir.toNativeSeparators(self.dataEntry4.text())
                if self.instrument==13:
                    rawdata=self.dir_path+'\\'+name+'.csv'
                else:
                    rawdata=self.dir_path+'\\'+name+'.d'+'\\'+name+'.csv'
                logging.info('Sample Name:%s',name)
                self.filename=name
                logging.info('File name:%s',self.filename)
                Del206=eval(self.dewll206le.text())
                Del207=eval(self.dewll207le.text())
                Del232=eval(self.dewll232le.text())
                Del238=eval(self.dewll238le.text())
                
                Std76Bias=eval(self.Stdbias76e.text())
                Std68Bias=eval(self.Stdbias68e.text()) 
                AblCorr=eval(self.Ablre.text())
                StartRow=int(self.Signalfe.text())
                EndRow=int(self.Signalte.text())
                BslErr=eval(self.bsere.text())/100
                ComPbFilter=self.CommPbrb1.isChecked()
                Interpol=self.Interpolrb1.isChecked()
                LowestUPb=eval(self.minU_Pble.text())
                bcg0=int(self.Blankfe.text())
                bcg1=int(self.Blankte.text())
                MinU=eval(self.MinUle.text())
                MinCt=eval(self.MinUCTle.text())
                MaxSWD=eval(self.Maxmswdle.text())
                
      
                
                
                
                
                logging.info('Del206:%s',Del206)
                logging.info('Del207:%s',Del207)
                logging.info('Del232:%s',Del232)
                logging.info('Del238:%s',Del238)
                logging.info('Std76Bias:%s',Std76Bias)
                logging.info('Std68Bias:%s',Std68Bias)
                logging.info('AblCorr:%s',AblCorr)
                logging.info('bcg start:%s',bcg0)
                logging.info('bcg end:%s',bcg1)
                logging.info('Signal start:%s',StartRow)
                logging.info('Signal end:%s',EndRow)
                logging.info('Baseline Erro:%s',BslErr)
                logging.info('Min. U:%s',MinU)
                logging.info('Min. U Pb(counts):%s',MinCt)
                logging.info('Lowest U/Pb:%s',LowestUPb)        
                logging.info('Max MSWD:%s',MaxSWD)
                logging.info('ComPbFilter:%s',ComPbFilter)
                logging.info('Interpol:%s',Interpol)
                if self.dataEntry3.text():            
                    outputpath=QtCore.QDir.toNativeSeparators(self.dataEntry3.text())
                    self.cwd=outputpath+'/'+self.filename+'_result'
                else:
                    outputpath=os.getcwd()
                    self.cwd=outputpath+'/'+self.filename+'_result'
                if not os.path.exists(self.cwd):
                    os.makedirs(self.cwd)
                Delspeak=self.Delspikes.isChecked()
                Mult=eval(self.Reject_levele.text())
                N=eval(self.Detectors_Noe.text())
                self.a=Bayesian_Regression(self.instrument,rawdata,name,Del206,Del207,Del232,Del238,Std76Bias,Std68Bias,AblCorr,StartRow,EndRow,BslErr,ComPbFilter,Interpol,
                              LowestUPb,bcg0,bcg1,MinU,MinCt,MaxSWD,N,Delspeak,Mult,InitAge=130,InitAgeErr=9.5,InitSlope=0.00753,InitSlopeErr=0.00052,
                              WethMSWD=1,AgeIncr=1,PbCIncr=0.01,outputpath=outputpath) 
              
                
                m=self.a.couts()
    
                   
                if self.dataEntry2.text():
                    logfile=QtCore.QDir.toNativeSeparators(self.dataEntry2.text())
                    
                    self.a.fliterdata(logfile=logfile)
    
                else:            
                    self.a.fliterdata()
    
                    
                self.a.rejectdata() 
    
                
                if BayesinC:
                    print('-'*35)
                    logging.info('Regression Start.')
                    #print("\r", end="")
                    #print(" {}%: ".format(numi), "█" * (numi // 2), end="",flush=True)
                    AgeStart=eval(self.Agesle.text())
                    AgeEnd=eval(self.Ageele.text())
                    DelAge=(AgeEnd-AgeStart)/int(self.Agenle.text())      
                   
                    Pb76Start=eval(self.Pbcsle.text())
                    Pb76End=eval(self.Pbcele.text())
                    DelPb76=(Pb76End-Pb76Start)/int(self.Pbcnle.text())
                    #print(AgeStart,AgeEnd)
    
                    logging.info('Age Start:%s',AgeStart)
                    logging.info('Age End:%s',AgeEnd)
                    logging.info('Age Num.:%s',int(self.Agenle.text()) )
                    
                    logging.info('Pbc76 Start:%s',Pb76Start)
                    logging.info('Pbc76 End:%s',Pb76End)
                    logging.info('Pbc76 Num.:%s',int(self.Pbcnle.text()) )
                                   
                    self.a.regress(Pb76Start,Pb76End,AgeStart,AgeEnd,DelPb76,DelAge,pltshowT=BayesinC)
    
    
                    
                    end_time=time.time()
                    print()
                    with open(self.cwd+'\\'+self.filename+'_DataInfo.txt','a') as f:
                        print('#Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))),file=f)
                        
                        print('*'*50,file=f)
                        
                    print('Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))))  
                    #print('Duration of processing:{:.2f} s'.format(end_time-self.start_time))  
                    print('*'*30)
                    
                    
                    
                    logging.info('Regression Compeleted!')
                    
                else:
                    logging.info('Calculation Compeleted!')
                    logging.info('Result saved in:%s',self.cwd)
                
            msg_box_c=QMessageBox(QMessageBox.Information,'Information','Calculat Isotope successful!')
            msg_box_c.exec_()
        except Exception as e:
            print(str(e))     
    
    def butnstate(self,btn):
        if btn.text()=='Agilent1':
            if btn.isChecked()==True:
                self.instrument=3
                print(btn.text()+' is selected')
            else:
                print()
                #print(btn.text()+' is deselected')
        elif btn.text()=='Agilent2':
            if btn.isChecked()==True:
                self.instrument=2
                print(btn.text()+' is selected')
            else:
                print()
                #print(btn.text()+' is deselected')
        elif btn.text()=='Thermo':
            
            if btn.isChecked()==True:
                self.instrument=13
                print(btn.text()+' is selected')
            else:
                print()
                #print(btn.text()+' is deselected')
    def get_dirctory(self,entry):
        pwd=os.getcwd()
        self.dir_path=QFileDialog.getExistingDirectory(self,'OPen Path',pwd)
        if self.dir_path:
            entry.setText(self.dir_path)                
      
    def plotshow(self):
        print('Plotting')
        
        try:
            rawdata=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
            data_all=np.loadtxt(rawdata,dtype=str,delimiter=',',skiprows=self.instrument,comments='     ') 
            self.data_Num=data_all.shape[0]
            if self.data_Num:
                self.Signalte.setText(str(self.data_Num))
        #data_all=np.loadtxt(rawdata,dtype=str,delimiter=',',skiprows=2,comments='     ')
        
            
            col_data_name=data_all[0,:].tolist()
        #print(col_data_name)
            col_name=[]
            for i in range(len(col_data_name)):
                col_name.append(''.join(list(filter(str.isdigit, col_data_name[i]))))
                if self.instrument==13:
                    Time=data_all[2:,0].astype(float)
                else:
                    Time=data_all[1:,0].astype(float)
            if self.instrument==13:
                
                Pb206=data_all[2:,col_name.index('206')].astype(float)
                Pb207=data_all[2:,col_name.index('207')].astype(float)                  
                U238=data_all[2:,col_name.index('238')].astype(float)
                try:
                    Th232=data_all[2:,col_name.index('232')]
                except:
                    pass
            else:
                Pb206=data_all[1:,col_name.index('206')].astype(float)
                Pb207=data_all[1:,col_name.index('207')].astype(float)                  
                U238=data_all[1:,col_name.index('238')].astype(float)
                try:
                    Th232=data_all[1:,col_name.index('232')]
                except:
                    pass
            plt.title('Rows vs. Signal(cps)')

            plt.plot(Pb206,label='206Pb')
            plt.plot(Pb207,label='207Pb')

            plt.plot(U238,label='238U')
            ax=plt.gca()
            ax.set_yscale('log')
            plt.legend()
            plt.grid(axis='x')
            plt.tight_layout()
            plt.xlabel("Rows")
            plt.ylabel("Counts (cps)")
        
        
            plt.show()
            
            if self.dataEntry1.text():
                if self.instrument==2 or self.instrument==3:
                    traspath1=os.path.dirname(self.dataEntry1.text())
                    traspath=os.path.dirname(traspath1)
                elif self.instrument==13:
                    traspath=os.path.dirname(self.dataEntry1.text())
               
                self.dataEntry4.setText(traspath)
            
        except Exception as e:
            print(str(e))
            logging.info('An error occurred when loading data, please check the data format!%s',str(e))
    def Merged_data_plot(self):
        
        try:
            
            Dir_path=QtCore.QDir.toNativeSeparators(self.dataEntry3.text())
    
            listfile=QtCore.QDir.toNativeSeparators(self.dataEntry2.text())
            filename_result=Dir_path+'\\'+'Result_merged.csv'
            Result_all_xls=[]
            
            
                
            for name in self.Sampleslist1.keys():
                
                workbook=xlrd.open_workbook(Dir_path+'\\'+name+'_result'+'\\'+name+'_Result.xls')
                workbook.sheet_names()
                worksheet=workbook.sheet_by_name('U-Pb Isotope results')
                for i in range(1,worksheet.nrows):
                    line=[]
                    for j in range(32):            
                        line.append(worksheet.cell(i,j).value)                  
                    Result_all_xls.append(line)
            result_all_np=np.array(Result_all_xls) 
            np.savetxt(filename_result,result_all_np,delimiter=',',fmt="%s")      
    
            savename=(list(filename_result.split('\\'))[-1]).split('.')[0]
            pd_result=np.loadtxt(filename_result,dtype=str,delimiter=',',skiprows=1)
            x=pd_result[:,2].astype(float)
            y=pd_result[:,3].astype(float)
            Pb206=pd_result[:,4].astype(float)
            Pb207=pd_result[:,5].astype(float)
            U238=pd_result[:,6].astype(float)
            R76=pd_result[:,17].astype(float)
            R86=pd_result[:,15].astype(float)
            
            Age68=pd_result[:,28].astype(float)
            Age75=pd_result[:,30].astype(float)
            
            xi=np.linspace(x.min(),x.max(),100)
            yi=np.linspace(y.min(),y.max(),100)
            
            X,Y=np.meshgrid(xi,yi)
            

            
            fig=plt.figure(figsize=(6,13),dpi=300)
            Z206=griddata((x,y),Pb206,(X,Y),method='nearest')
            
            ax1=fig.add_subplot(321,projection='3d')
            ax1.plot_surface(X,Y,Z206,cmap='jet')
            ax1.contour(X,Y,Z206,cmap='jet')
            
            
            
            ax2=fig.add_subplot(322)
            ax2.contourf(X,Y,Z206,cmap='jet')
            plt.title('Pb206 \n', 
                      fontsize=14, fontweight='bold') 
            
            
            
            
            Z207=griddata((x,y),Pb207,(X,Y),method='cubic')
            
            ax3=fig.add_subplot(323,projection='3d')
            ax3.plot_surface(X,Y,Z207,cmap='jet')
            ax3.contour(X,Y,Z207,cmap='jet')
            
            
            ax4=fig.add_subplot(324)
            ax4.contourf(X,Y,Z207,cmap='jet')
            plt.title('Pb207 \n', 
                      fontsize=14, fontweight='bold') 
            
            
            Z238=griddata((x,y),U238,(X,Y),method='cubic')
            
            ax5=fig.add_subplot(325,projection='3d')
            ax5.plot_surface(X,Y,Z238,cmap='jet')
            ax5.contour(X,Y,Z238,cmap='jet')
            
            
            ax6=fig.add_subplot(326)
            ax6.contourf(X,Y,Z238,cmap='jet')
            
            plt.title('U238 \n', 
                      fontsize=14, fontweight='bold') 
            
            plt.tight_layout()
            plt.savefig(Dir_path+'\\'+savename+'cmp_678')
            
            
            
            fig2=plt.figure(figsize=(6,9),dpi=300)
            Z76=griddata((x,y),R76,(X,Y),method='cubic')
            
            ax7=fig2.add_subplot(221,projection='3d')
            ax7.plot_surface(X,Y,Z76,cmap='jet')
            ax7.contour(X,Y,Z76,cmap='jet')
            
            
            
            ax8=fig2.add_subplot(222)
            ax8.contourf(X,Y,Z76,cmap='jet')
            plt.title('Pb207/Pb206 \n', 
                      fontsize=14, fontweight='bold') 
            
            Z86=griddata((x,y),R86,(X,Y),method='cubic')
            #lim86=np.arange(50,100,1)
            ax7=fig2.add_subplot(223,projection='3d')
            ax7.plot_surface(X,Y,Z86,cmap='jet')
            ax7.contour(X,Y,Z86,cmap='jet')
            
            
            
            ax8=fig2.add_subplot(224)
            ax8.contourf(X,Y,Z86,cmap='jet')
            plt.title('U238/Pb206 \n', 
                      fontsize=14, fontweight='bold') 
            
            plt.tight_layout()
            plt.savefig(Dir_path+'\\'+savename+'cmp_R678')
            
            #lim=np.arange(50,110,1)
            fig3=plt.figure(figsize=(6,9),dpi=300)
            Z75=griddata((x,y),Age75,(X,Y),method='cubic')
            
            ax7=fig3.add_subplot(221,projection='3d')
            ax7.plot_surface(X,Y,Z75,cmap='jet')
            ax7.contour(X,Y,Z75,cmap='jet')
            
            
            
            ax8=fig3.add_subplot(222)
            ax8.contourf(X,Y,Z75,cmap='jet')
            plt.title('Age(Pb207/Pb206) \n', 
                      fontsize=14, fontweight='bold') 
            
            Z68=griddata((x,y),Age68,(X,Y),method='cubic')
            
            ax7=fig3.add_subplot(223,projection='3d')
            ax7.plot_surface(X,Y,Z68,cmap='jet')
            ax7.contour(X,Y,Z68,cmap='jet')
            
            
            
            ax8=fig3.add_subplot(224)
            ax8.contourf(X,Y,Z68,cmap='jet')
            
            plt.title('Age(Pb206/U238) \n', 
                      fontsize=14, fontweight='bold') 
            
            plt.tight_layout()
            plt.savefig(Dir_path+'\\'+savename+'cmp_Age7568')
            
            
            
            plt.show() 
            msg_box_c=QMessageBox(QMessageBox.Information,'Information','Plotting successful!')
            msg_box_c.exec_()
        except Exception as e:
            print(str(e))
            
          
    def get_file_xls(self,entry):
        pwd=os.getcwd()
        logging.info('Current work directory:%s',pwd)
        filename, filters = QFileDialog.getOpenFileName(self,'Choose *.xls',pwd,'Map file (*.xls)' )
        
        logging.info('File:%s',filename)
        try:
            if filename:
                entry.setText(filename) 
            
            

            data = xlrd.open_workbook(filename)
            data.sheet_names()
            ws = data.sheet_by_name('Sheet1')
            Sampleslist1={}
            for i in range(1,ws.nrows):
                if isinstance(ws.cell(i,1).value ,float):
                    Sampleslist1[ws.cell(i,0).value]='%d'%(ws.cell(i,1).value )
                else:
                     Sampleslist1[ws.cell(i,0).value]='%s'%(ws.cell(i,1).value )
            print(Sampleslist1)
            
            self.Sampleslist1=Sampleslist1
            logging.info('File:%s',filename)
            if filename:
                entry.setText(filename)
        except Exception as e:
            print(str(e))
       
        
    def get_file(self, entry):            
        pwd=os.getcwd()
        logging.info('Current work directory:%s',pwd)
        filename, filters = QFileDialog.getOpenFileName(self,'Choose *.csv',pwd,'Map file (*.csv)::' )
        logging.info('File:%s',filename)
        
        if filename:
            entry.setText(filename)         
    def Cal_Istope(self):        
        try:
            logging.info('Calculation Isotope start.')
            rawdata=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
            self.filename=(list(rawdata.split('\\'))[-1]).split('.')[0]
            logging.info('File name:%s',self.filename)   
            
            
            
            self.settings.setValue('Dwell_Time/Pb206',self.dewll206le.text())
            self.settings.setValue('Dwell_Time/Pb207',self.dewll207le.text())
            self.settings.setValue('Dwell_Time/U238',self.dewll238le.text())
            self.settings.setValue('Dwell_Time/Others',self.dewll232le.text())
            
            
            Del206=eval(self.dewll206le.text())
            Del207=eval(self.dewll207le.text())
            Del232=eval(self.dewll232le.text())
            Del238=eval(self.dewll238le.text())
            self.settings.setValue('Std_biase_factor/Pb207-Pb206',self.Stdbias76e.text())
            self.settings.setValue('Std_biase_factor/Pb206-U238',self.Stdbias68e.text())
            self.settings.setValue('Std_biase_factor/Ablration_factor',self.Ablre.text())
            
            Std76Bias=eval(self.Stdbias76e.text())
            Std68Bias=eval(self.Stdbias68e.text()) 
            AblCorr=eval(self.Ablre.text())
            
            self.settings.setValue('Signal_setting/Backgroud_Start',self.Blankfe.text())
            self.settings.setValue('Signal_setting/Backgroud_End',self.Blankte.text())
            self.settings.setValue('Signal_setting/Sampling_Start',self.Signalfe.text())
            self.settings.setValue('Signal_setting/Sampling_End',self.Signalte.text())
            
            self.settings.setValue('Fileter_setting/Baseline_Erro',self.bsere.text())
            self.settings.setValue('Fileter_setting/MinU-Pb',self.minU_Pble.text())

            self.settings.setValue('Fileter_setting/MinU_counts',self.MinUle.text())
            self.settings.setValue('Fileter_setting/MinPb_counts',self.MinUCTle.text())
            self.settings.setValue('Fileter_setting/MaxMSWD',self.Maxmswdle.text())
            
            
            
            StartRow=int(self.Signalfe.text())
            EndRow=int(self.Signalte.text())
            BslErr=eval(self.bsere.text())/100
            ComPbFilter=self.CommPbrb1.isChecked()
            Interpol=self.Interpolrb1.isChecked()
            LowestUPb=eval(self.minU_Pble.text())
            bcg0=int(self.Blankfe.text())
            bcg1=int(self.Blankte.text())
            MinU=eval(self.MinUle.text())
            MinCt=eval(self.MinUCTle.text())
            MaxSWD=eval(self.Maxmswdle.text())
            
            logging.info('Del206:%s',Del206)
            logging.info('Del207:%s',Del207)
            logging.info('Del232:%s',Del232)
            logging.info('Del238:%s',Del238)
            logging.info('Std76Bias:%s',Std76Bias)
            logging.info('Std68Bias:%s',Std68Bias)
            logging.info('AblCorr:%s',AblCorr)
            logging.info('bcg start:%s',bcg0)
            logging.info('bcg end:%s',bcg1)
            logging.info('Signal start:%s',StartRow)
            logging.info('Signal end:%s',EndRow)
            logging.info('Baseline Erro:%s',BslErr)
            logging.info('Min. U:%s',MinU)
            logging.info('Min. U Pb(counts):%s',MinCt)
            logging.info('Lowest U/Pb:%s',LowestUPb)        
            logging.info('Max MSWD:%s',MaxSWD)
            logging.info('ComPbFilter:%s',ComPbFilter)
            logging.info('Interpol:%s',Interpol)
            


            
            if self.dataEntry3.text():            
                outputpath=QtCore.QDir.toNativeSeparators(self.dataEntry3.text())
                self.cwd=outputpath+'/'+self.filename+'_result'
                
            else:
                outputpath=os.getcwd()
                self.cwd=outputpath+'/'+self.filename+'_result'
            if not os.path.exists(self.cwd):
                os.makedirs(self.cwd)
            name=None
            Delspeak=self.Delspikes.isChecked()
            Mult=eval(self.Reject_levele.text())
            N=eval(self.Detectors_Noe.text())
               
        
            self.a=Bayesian_Regression(self.instrument,rawdata,name,Del206,Del207,Del232,Del238,Std76Bias,Std68Bias,AblCorr,StartRow,EndRow,BslErr,ComPbFilter,Interpol,
                   LowestUPb,bcg0,bcg1,MinU,MinCt,MaxSWD,N,Delspeak,Mult,InitAge=130,InitAgeErr=9.5,InitSlope=0.00753,InitSlopeErr=0.00052,
                   WethMSWD=1,AgeIncr=1,PbCIncr=0.01,outputpath=outputpath)     
            m=self.a.couts()
            if self.dataEntry2.text():
                logfile=QtCore.QDir.toNativeSeparators(self.dataEntry2.text())
                self.Agestart,self.Agesig,self.Pbstart,self.Pbsig=self.a.fliterdata(logfile=logfile)
            else:            
                
                self.Agestart,self.Agesig,self.Pbstart,self.Pbsig=self.a.fliterdata()
            self.a.rejectdata() 
        
            
        except Exception as e:
            print(str(e))    
        msg_box_c=QMessageBox(QMessageBox.Information,'Information','Calculat Isotope successful!')
        msg_box_c.exec_()
        logging.info('Calculation Compeleted!')
        logging.info('Result saved in:%s',self.cwd)
        
    def regression(self):
        try:
            print('-'*35)
            logging.info('Regression Start.')
            #print("\r", end="")
            #print(" {}%: ".format(numi), "█" * (numi // 2), end="",flush=True)
            AgeStart=eval(self.Agesle.text())
            AgeEnd=eval(self.Ageele.text())
            DelAge=(AgeEnd-AgeStart)/int(self.Agenle.text())      
           
            Pb76Start=eval(self.Pbcsle.text())
            Pb76End=eval(self.Pbcele.text())
            DelPb76=(Pb76End-Pb76Start)/int(self.Pbcnle.text())
            #print(AgeStart,AgeEnd)

            logging.info('Age Start:%s',AgeStart)
            logging.info('Age End:%s',AgeEnd)
            logging.info('Age Num.:%s',int(self.Agenle.text()) )
            
            logging.info('Pbc76 Start:%s',Pb76Start)
            logging.info('Pbc76 End:%s',Pb76End)
            logging.info('Pbc76 Num.:%s',int(self.Pbcnle.text()) )
            
            self.a.regress(Pb76Start,Pb76End,AgeStart,AgeEnd,DelPb76,DelAge,pltshowT=False)
            

            
            end_time=time.time()
            print()
            with open(self.cwd+'/'+self.filename+'_DataInfo.txt','a') as f:
                print('#Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))),file=f)
                
                print('*'*50,file=f)
                
            print('Ending time:{}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))))  
            #print('Duration of processing:{:.2f} s'.format(end_time-self.start_time))  
            print('*'*30)
            self.settings.setValue('Bayesian_setting/Age_Start',self.Agesle.text())
            self.settings.setValue('Bayesian_setting/Age_End',self.Ageele.text())
            self.settings.setValue('Bayesian_setting/Pbc_Start',self.Pbcsle.text())
            self.settings.setValue('Bayesian_setting/Pbc_End',self.Pbcele.text())
            
            msg_box_r=QMessageBox(QMessageBox.Information,'Information','Regression successful!')
            msg_box_r.exec_()
            logging.info('Regression Compeleted!')
        except Exception as e:
            print(str(e))
        
        
if __name__=='__main__':
    logging.info('Brama Start.')

    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QIcon('./images/home.png'))    
    main=MyApplication()
    main.show()
    sys.exit(app.exec_())
    
