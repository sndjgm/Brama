import re
import csv
import os
from math import *
import time
import pandas as pd
import numpy as np
import xlwt,xlrd
import random 
import matplotlib.pyplot as plt
import datetime
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon,QTextCursor
from PyQt5 import QtCore
from xlrd import *




class DataReader():
    """
    Reads data into pandas Dataframe from a file.

    Parameters
    ----------
    filename: str
        Path to file to read.
    filetype : str
        Type of the file to read. If not specified, csv is used. 
        Possible options are csv, xlsx and asc.
    instrument : str
        Type of the instrument used for measurement. If not specified, raw data is expected. 
        Possible options are Agilent, and Element.
    """

    def __init__(self, filename, filetype=None, instrument=None):

        self.filename = filename
        self.filetype = filetype
        self.instrument = instrument
        self.data = self.read(filename, filetype, instrument)

    def __call__(self):
        return self.data

    def read(self, filename, filetype, instrument):
        """
        Reads data into pandas Dataframe from a file.

        Parameters
        ----------
        filename: str
            Path to file to read.
        filetype : str
            Type of the file to read. If not specified, csv is used. 
            Possible options are csv, xlsx and asc.
        instrument : str
            Type of the instrument used for measurement. If not specified, raw data is expected. 
            Possible options are Agilent, and Element.

        Returns
        -------
        data : dataframe
            `data` as a dataframe, which can be passed to MSData.
        """

        #if instrument == 'Element':
            #skipfooter = 4
            #header = 1
            #drop = 9
            #skiprows=0
        if instrument == 'Agilent1':
            
            skipfooter = 4
            header = 3
            drop = 3
            skiprows=0
            data = pd.read_csv(filename, sep=',', index_col=0, skipfooter=skipfooter,skiprows=skiprows,
                                header=header, engine='python')

            
        elif instrument == 'Agilent2':
            skipfooter = 3
            header = 3
            drop = 3
            skiprows=0
            data = pd.read_csv(filename, sep=',', index_col=0, skipfooter=skipfooter,skiprows=skiprows,
                               header=header, engine='python')
        elif instrument=='Thermo':
            skipfooter = 12
            header=13
            skiprows=1
            drop=12
            data = pd.read_csv(filename, sep=',', index_col=0, skipfooter=skipfooter,skiprows=skiprows,
                               header=header, engine='python')
        elif instrument=='Guanzhou':
            skipfooter = 12
            header=13
            skiprows=1
            drop=12
            data = pd.read_csv(filename, sep=',', index_col=0, skipfooter=skipfooter,skiprows=skiprows,
                               header=header, engine='python')
            
                
        else:
            skipfooter = 0
            header = 0
            drop = 0

        #if filetype == 'xlsx':
            #imported = pd.ExcelFile(filename)
            #data = imported.parse(
                #0, index_col=0, skipfooter=skipfooter, header=header)
            #data = data.drop(data.index[:drop], axis=0)

        #elif filetype == 'csv':
            

            
        #elif filetype == 'asc':
            #data = pd.read_csv(filename, sep='\t', index_col=0, skipfooter=skipfooter,
                               #header=header, engine='python')
            #data = data.drop(data.index[:drop], axis=0)
            #data.dropna(axis=1, how='all', inplace=True)
            #data = data.apply(pd.to_numeric, errors='coerce')

        #else:
            #warnings.warn('File type not supported.')

        return data


class Selector():
    """
    Class for peak and background identification.

    Parameters
    ----------
    ms_data : MSData
        Mass spectrometry data class, for which the selector will be used.
    s : float (Optional)
        Start of first peak in seconds for synchronisation with Iolite data if using 
        method Iolite or for calculation  of backround standard deviation if using 
        method treshold. Default is 60s.
    sdmul: int (Optional)
        Multiplier for backround standard deviation in method treshold. Default is 10. 
    iolite: MSEval.Iolite (Optional)
        Iolite class, only necessary if method is Iolite. Default is None.
    logger: logger class (optional)
        If logger is pssed all methods of Selector will log in the activity.
    """

    def __init__(self, ms_data,times, s=60, sdmul=10, iolite=None, logger=None):
        self.possible_methods = ('treshold', 'iolite')
        self.ms_data = ms_data
        self.times=times
        self.logger = logger
        self.filter_line = self.ms_data.data.sum(1)
        self.method = 'treshold'
        self.iolite = iolite
        if s>0:
            self.start = self.ms_data.data.index[get_index(self.ms_data.data, s)]
        else:
            self.start =int(s/(self.ms_data.data.index[1]-self.ms_data.data.index[0]))
            #print(self.start)
            
            
        #print(self.start)
        if self.logger is not None:
            self.logger.info(f'Starting point found at {self.start}s.')
        self.sdmul = sdmul
        self.skip = {'bcg_start': 0,
                     'bcg_end': 0,
                     'sample_start': 0,
                     'sample_end': 0}    # time in seconds to skip from each bcg and sample

    def __call__(self):
        if self.method == 'treshold':
            return self.create_selector_treshold()
        elif self.method == 'iolite':
            return self.create_selector_iolite()

    def create_selector_iolite(self):
        """
        Select starts and ends of ablation using iolite file.
        """



        

        starts = []
        ends = []
        XX=[]
        YY=[]
         
        for i, t in enumerate(self.times[1:]):            
            timeindex = get_difference(get_timestamp(self.times[0]), get_timestamp(t)) + self.start            
            if i % 2 == 0:
                try:
                    starts.append(get_index(self.ms_data.data, timeindex))
                except:
                    starts.append(0)
            if i % 2 != 0:
                try:
                    ends.append(get_index(self.ms_data.data, timeindex)+1)
                except TypeError:
                    ends.append(get_index(self.ms_data.data, timeindex))
       
        return starts, ends

    def create_selector_treshold(self):
        """
        Select starts and ends of ablation based on selected element or sum of all using treshold
        calculated from background.
        """

        if self.logger is not None:
            self.logger.info('Selecting peak bounds by setting treshold.')
        bcg_nr = self.ms_data.time_to_number(self.start)
        bcg = self.filter_line.iloc[0:bcg_nr].mean()
        std = self.filter_line.iloc[0:bcg_nr].std()
        ind = [True if value > bcg+self.sdmul *
               std else False for value in self.filter_line]
        ind2 = ind[1:]
        ind2.append(False)
        index = [i for i in range(0, len(ind)) if ind[i] != ind2[i]]

        starts = [index[i] for i in range(len(index)) if i % 2 == 0]
        ends = [index[i] for i in range(len(index)) if i % 2 != 0]

        return starts, ends

    def set_skip(self, bcg_s=None, bcg_e=None, sig_s=None, sig_e=None):
        """
        Set time skipped on start and end of background and ablation in seconds.
        """

        if bcg_s is not None:
            self.skip['bcg_start'] = bcg_s
        if bcg_e is not None:
            self.skip['bcg_end'] = bcg_e
        if sig_s is not None:
            self.skip['sample_start'] = sig_s
        if sig_e is not None:
            self.skip['sample_end'] = sig_e

    def create_on_off(self, starts, ends):
        """
        From starts and ends of ablation create laser_on and laser_off with skipped values.
        """

        laser_off = []
        laser_on = []

        laser_off.append((0+self.ms_data.time_to_number(
            self.skip['bcg_start']), starts[0]-self.ms_data.time_to_number(self.skip['bcg_end'])))

        for i in range(len(starts)-1):
            laser_off.append((ends[i]+self.ms_data.time_to_number(self.skip['bcg_start']),
                              starts[i+1]-self.ms_data.time_to_number(self.skip['bcg_end'])))
            laser_on.append((starts[i]+self.ms_data.time_to_number(self.skip['sample_start']),
                             ends[i]-self.ms_data.time_to_number(self.skip['sample_end'])))

        laser_off.append((ends[-1]+self.ms_data.time_to_number(self.skip['bcg_start']), len(
            self.ms_data.time)-2-self.ms_data.time_to_number(self.skip['bcg_end'])))
        laser_on.append((starts[-1]+self.ms_data.time_to_number(self.skip['sample_start']),
                         ends[-1]-self.ms_data.time_to_number(self.skip['sample_end'])))
        return laser_on, laser_off
class Background():
    """
    Class for calculating and visualising background.

    Parameters
    ----------
    isotope : MSData.Isotope
        Isotope from which background is calculated
    laser_on : list
        List of sets, where each set has 2 values (start of the laser and
        end of the laser.) indicating when the laser was fireing.
    laser_off : list
        List of sets, where each set has 2 values (end of the laser and
        start of the laser.) indicating when the laser wasn't fireing.    
    offset : float (Optional)
        Parameter for outlier removal. 
        1-offset = upper treshold for percentile filtering. 
        Accepts values between 0 and 1. Default is 0.15.
    width : float (Optional)
        Parameter for outlier removal.  Width od returned values, 
        where 1-offset-width = lower treshold for percentile filtering.
        Accepts values between 0 and 1. Default is 0.8.
    """

    def __init__(self, isotope, laser_on, laser_off, width=0.8, offset=0.15):
        assert isotope.data is not None

        # background by means
        self.bcg_all = []
        for i in range(len(laser_off)):
            data_points = isotope.data[laser_off[i][0]:laser_off[i][1]]
            if len(data_points) == 0:
                data_points = [0]
            if len(data_points) > 0:
                data_points = remove_outliers(data_points, offset, width)
            self.bcg_all.append(data_points)

        self.bcg_means = [x.mean() for x in self.bcg_all]

        # background by interpolation
        bcg_mskd = np.copy(isotope.data)
        for (s, e) in laser_on:
            bcg_mskd[s-2:e+2] = np.nan
        not_nan = np.logical_not(np.isnan(bcg_mskd))
        indices = np.arange(len(bcg_mskd))
        self.bcg_interp = np.interp(
            indices, indices[not_nan], bcg_mskd[not_nan])

    def __call__(self):
        plt.plot(self.bcg_means)

    def __repr__(self):
        return f'{self.__class__.__name__}: {self.bcg_means} '
class Iolite():
    """
    Class holding Iolite data.

    Parameters
    ----------
    path : str
        Path to Iolite .csv file.
    """

    def __init__(self, path,laser='Resolution'):
        self.data = self.read_iolite(path)
        self.laser=laser
        self.peak_names = self.names_from_iolite(laser)

    def read_iolite(self, path):
        """
        Read Iolite file.

        Parameters
        ----------
        path : str
            Path to Iolite .csv file.
        """
        iolite = pd.read_csv(path, sep=",", engine='python')
        #print(iolite)
        return iolite

    def names_from_iolite(self,laser):
        """
        Return peak names from iolite. 
        """
        
        #names = list(self.data[' Comment'].dropna())
        if laser=='Resolution':
            names = list(self.data[' Comment'].dropna())                   
            
        elif laser=='Tandam':
            print('Tandam  log')
            names = list(self.data[self.data[' Laser State']=='On'][' Comment'])
        return names

    def on_and_off_times(self):
       
            
        times = self.data['Timestamp'][self.data.ne(self.data.shift()).apply(lambda x: x.index[x].tolist())[' Laser State']]
    
        return times

        
    
    def coord(self):
        """
        Return X Y Z position      

        """
        if self.laser=='Tandam':
            XX=self.data[' X(um)']
            YY=self.data[' Y(um)']
        else:
            try:            
                XX=self.data[self.data[' Laser State']=='On'][' X(um)']
                YY=self.data[self.data[' Laser State']=='On'][' Y(um)']
                #print(XX,YY)
            except Exception:
                XX=self.data[self.data[' Laser State']=='On'][' X']
                YY=self.data[self.data[' Laser State']=='On'][' Y']

            
        XX0=[]
        XX1=[]
        YY0=[]
        YY1=[]
        print()
        
        
        names=self.names_from_iolite(self.laser)
        Timeslen=len(XX)/len(names)
        restlen=len(XX)%len(names)
        print(Timeslen,restlen)
        if Timeslen==1:
            if self.laser=='Tandam':
                for i in range(len(XX)):
                    if i==0:
                        pass
                    else:
                        if i%2 ==1:                            
                             XX0.append(XX.iloc[i])
                             YY0.append(YY.iloc[i]) 
                        elif i % 2 ==0:
                             XX1.append(XX.iloc[i])
                             YY1.append(YY.iloc[i])

            else:
                XX0=XX.tolist()
                XX1=XX.tolist()
                YY0=YY.tolist()
                YY1=YY.tolist()
            
        elif Timeslen==2:
            for i in range(len(XX)):
                if i %2 ==0:
                    XX0.append(XX.iloc[i])
                    YY0.append(YY.iloc[i])
                if i % 2 !=0:
                    XX1.append(XX.iloc[i])
                    YY1.append(YY.iloc[i])
        elif Timeslen==3:
            
            for i in range(len(XX)):
                if i%3 ==0:
                    XX0.append(XX.iloc[i])
                    YY0.append(YY.iloc[i])
                if i%3==1:
                    pass
                else:
                    XX1.append(XX.iloc[i])
                    YY1.append(YY.iloc[i])
        else:
            if self.laser=='Tandam':
                for i in range(len(XX)):
                    if i==0:
                        pass
                    else:
                        if i%2 ==1:                            
                             XX0.append(XX.iloc[i])
                             YY0.append(YY.iloc[i]) 
                        elif i % 2 ==0:
                             XX1.append(XX.iloc[i])
                             YY1.append(YY.iloc[i])
            else:
                XX0=XX.tolist()
                XX1=XX.tolist()
                YY0=YY.tolist()
                YY1=YY.tolist()
                
                    
                
        
        #print(XX0,YY0,XX1,YY1)
        return XX0,YY0,XX1,YY1
    
    
def get_timestamp(strTime):
    
 
    time_new=time.mktime(time.strptime(strTime, '%Y-%m-%d %H:%M:%S.%f'))

        
    return time_new
def get_difference(start, now):
    
    diff = abs(now - start)
    
    return diff
def get_index(data, time):
    """return closest index of MS time given time in seconds"""
    try:    
        for i in range(len(data.index)-1):
            if (data.index[i] <= time) and (data.index[i+1] > time):
                
                return i
    except Exception:  
        return 0
class Iolite_to_Aglient:
    
    def __init__(self,laser='Resolution',log_path=r'D:\6.贝叶斯回归法\calcite line\20220722A_log_20220722_112023.csv',
                 ms_path=r"D:\6.贝叶斯回归法\calcite line\20220722A_Calcite_U-Pb_120um_3.5j_15HZ_20umS\20220722A.csv",
                 outputpath=r'D:\6.贝叶斯回归法\calcite line\output',diff_sec=1,width=30,Dir_file='A',instrument='Agilent1'):
        self.path=log_path
        self.filename=ms_path
        self.outputpath=outputpath
        self.Dir_file=Dir_file
        self.width=width
        self.s=diff_sec
        self.instrument=instrument
        self.laser=laser
    def creat_Aglient(self):
        a=Iolite(self.path,self.laser)
        iolite=a.read_iolite(self.path)
        names=a.names_from_iolite(self.laser)
        
        self.times=a.on_and_off_times()
        
        XX0,YY0,XX1,YY1=a.coord()
        
        
        datareader=DataReader(self.filename,filetype='csv',instrument=self.instrument)
        ms_data=datareader.read(self.filename,filetype='csv',instrument=self.instrument)
        
        
        
        b=Selector(datareader,times=self.times,s=self.s)
        starts,ends=b.create_selector_iolite()

        
        
        samplelist={}

        

        for i in range(len(names)):
            #print(names[i])
            
            samplelist[i+1]=[names[i],starts[i],ends[i],XX0[i],YY0[i],XX1[i],YY1[i]]  
            
             
        print(samplelist)   
        workbook = xlwt.Workbook()
        worksheet=workbook.add_sheet('Sheet1')
        style = xlwt.XFStyle() # ?????
        font = xlwt.Font() # ???????
        font.name = 'Times New Roman' 
        font.bold = True # ??
        font.underline = True # ???
        font.italic = True # ???
        style.font = font # ????          

        worksheet.write(0, 4, 'sequence Num')
        worksheet.write(0, 0, 'FileName')
        worksheet.write(0, 1, 'SampleName')

        worksheet.write(0, 2, 'Laser on')
        worksheet.write(0, 3, 'Laser off')

        worksheet.write(0, 5, 'X0')
        worksheet.write(0, 6, 'Y0')
        worksheet.write(0, 7, 'X1')
        worksheet.write(0, 8, 'Y1')

        for key in samplelist.keys():
            s1=samplelist[key][1]
            s2=samplelist[key][2]
            
            worksheet.write(int(key),0,self.Dir_file+str(key))
            worksheet.write(int(key),1,samplelist[key][0])
            worksheet.write(int(key),2,samplelist[key][1])
            worksheet.write(int(key),3,samplelist[key][2])
            
            worksheet.write(int(key),4,key)
            worksheet.write(int(key),5,int(samplelist[key][3]))
            worksheet.write(int(key),6,int(samplelist[key][4]))
            worksheet.write(int(key),7,int(samplelist[key][5]))
            worksheet.write(int(key),8,int(samplelist[key][6]))

            
            
            os.mkdir(os.path.join(self.outputpath,self.Dir_file+str(key)+'.D'))
            newpath=os.path.join(self.outputpath+'\\'+self.Dir_file+str(key)+'.D',self.Dir_file+str(key)+'.csv')
            h=open(newpath,'w',newline='')
            s=os.path.dirname(self.outputpath)
            s=[''.join(s)]            
            writer=csv.writer(h)            
            writer.writerow([str(samplelist[key][0])])
            
            
            
            if s1 and s2:
                datum=ms_data.iloc[abs(s1-self.width):s2+self.width]
            else:
                print('Index Out of range!')
                datum=ms_data.iloc[0:s2+self.width]
            
            Time_all=datum.index.tolist()
            if key==1:
                
                writer.writerow(['start:0'])
                
                writer.writerow(['Time']+list(datum.columns.values))
                datum.insert(loc=0,column='Time',value=Time_all)
                for i in range(datum.shape[0]):
                    writer.writerow(datum.iloc[i])                
                h.close()
            else:
                Time0=Time_all[0]
                Time_all=[i-Time0 for i in Time_all]
                writer.writerow(['start:'+str(Time0)])
                
                writer.writerow(['Time']+list(datum.columns.values))
                datum.insert(loc=0,column='Time',value=Time_all) 
                for i in range(datum.shape[0]):
                    writer.writerow(datum.iloc[i])                
                h.close()
        workbook.save(self.outputpath+'//'+self.Dir_file+'_LIST.xls')
        return samplelist 
        

        




class MyApplication(QWidget):
    def __init__(self):
        super().__init__()
        self.laser='Resolution'
        self.th=Iolite_to_Aglient(self.laser)
        
        

        self.initUI()
    def onUpdateText(self,text):
        cursor=self.process.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.process.setTextCursor(cursor)
        self.process.ensureCursorVisible()
        
    def initUI(self):
        self.setWindowTitle('Continuous Data Segmentation')
        #self.setFixedSize(450,350)
        self.resize(400,250)
        mainLayout=QGridLayout(self)
        
        
        Input_path_label=QLabel('MS file:',self)
        InputData = QPushButton('Browse',self)
        InputData.setIcon(QIcon('./images/load.jpg'))        
        self.dataEntry1 = QLineEdit(self)
        InputData.clicked.connect(lambda: self.get_file(self.dataEntry1))
        
        log_path_label=QLabel('Log file:',self)
        LogData = QPushButton('Browse',self)
        LogData.setIcon(QIcon('./images/load.jpg')) 
        self.dataEntry2 = QLineEdit(self)
        LogData.clicked.connect(lambda: self.get_file(self.dataEntry2))
        
        Output_path_label=QLabel('Export:',self)
        OutputData = QPushButton('Browse',self)
        OutputData.setIcon(QIcon('./images/load.jpg'))        
        self.dataEntry3 = QLineEdit(self)
        OutputData.clicked.connect(lambda: self.get_dirctory(self.dataEntry3))
        
        
        btnLOAD=QPushButton('Plot',self)
        self.btnLOAD=btnLOAD.clicked.connect(self.plotshow)
        
        btnConvert=QPushButton('Convert',self)
        self.Convert=btnConvert.clicked.connect(self.convertCMD)
        
        btnPOS=QPushButton('Position',self)
        self.Convert=btnPOS.clicked.connect(self.plotposition)
        
        start_timel=QLabel('&Start time(s)',self)
        self.start_timee=QLineEdit(self)
        self.start_timee.setText('1')
        start_timel.setBuddy(self.start_timee)
        extrowl=QLabel('&Extended Rows',self)
        self.extrowe=QLineEdit(self) 
        self.extrowe.setText('10')
        extrowl.setBuddy(self.extrowe)
        
        Dir_filel=QLabel('&Dir Filename',self)
        self.Dir_filee=QLineEdit(self)
        self.Dir_filee.setText('A')
        Dir_filel.setBuddy(self.Dir_filee)
        
        
        self.laser='Resolution'
        rblaser1=QRadioButton('Resolution',self)
        rblaser1.setChecked(True)
        rblaser1.toggled.connect(lambda: self.butlstate(rblaser1))
        
        rblaser2=QRadioButton('Tandam',self)
        rblaser2.toggled.connect(lambda: self.butlstate(rblaser2))
        lasergroup=QButtonGroup(self)
        lasergroup.addButton(rblaser1)
        lasergroup.addButton(rblaser2)
        
        
        rbAgilent1=QRadioButton('Agilent1',self)
        
        rbAgilent1.toggled.connect(lambda: self.butnstate(rbAgilent1))
        rbAgilent2=QRadioButton('Agilent2',self)
        rbAgilent2.setChecked(True)
        rbAgilent2.toggled.connect(lambda: self.butnstate(rbAgilent2))
        Thermo=QRadioButton('Thermo',self)
        Thermo.toggled.connect(lambda: self.butnstate(Thermo))
        instrgroup=QButtonGroup(self)
        instrgroup.addButton(rbAgilent1)
        instrgroup.addButton(rbAgilent2)
        instrgroup.addButton(Thermo)
        Instrument=QLabel('Instrument:',self)
        
        mainLayout.addWidget(rblaser1,0,4)  
        mainLayout.addWidget(rblaser2,0,5)  
        
        mainLayout.addWidget(Instrument,0,0)  
        mainLayout.addWidget(rbAgilent1,0,1)  
        mainLayout.addWidget(rbAgilent2,0,2)  
        mainLayout.addWidget(Thermo,0,3) 
        
        
        
        
        
        
        
        mainLayout.addWidget(Input_path_label,1,0)
        mainLayout.addWidget(log_path_label,2,0)
        mainLayout.addWidget(Output_path_label,3,0)
        
        mainLayout.addWidget(self.dataEntry1,1,1,1,3)
        mainLayout.addWidget(self.dataEntry2,2,1,1,3)
        mainLayout.addWidget(self.dataEntry3,3,1,1,3)
        
        
        mainLayout.addWidget(InputData,1,4,1,1)
        mainLayout.addWidget(LogData,2,4,1,1)
        mainLayout.addWidget(OutputData,3,4,1,1)
        
        mainLayout.addWidget(start_timel,4,0)
        mainLayout.addWidget(self.start_timee,4,1)
        
        mainLayout.addWidget(extrowl,4,2)
        mainLayout.addWidget(self.extrowe,4,3)
        
        mainLayout.addWidget(Dir_filel,4,4)
        mainLayout.addWidget(self.Dir_filee,4,5)
        
        mainLayout.addWidget(btnLOAD,1,5,3,1)
        
        mainLayout.addWidget(btnConvert,5,1,2,4)
        
        mainLayout.addWidget(btnPOS,6,1,2,4)
        self.instrument=2
    def butlstate(self,btn):
        if btn.text()=='Resolution':
            if btn.isChecked()==True:
                self.laser='Resolution'
                print(btn.text()+' is selected')
            else:
                print()
                #print(btn.text()+' is deselected')
        elif btn.text()=='Tandam':
            if btn.isChecked()==True:
                self.laser='Tandam'
                print(btn.text()+' is selected')
            else:
                print()
                #print(btn.text()+' is deselected')
                
        
        
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
    def plotposition(self):
        try:
            print(self.samplelist)
            Samplesinfo=pd.DataFrame(self.samplelist)
            
            plt.plot(Samplesinfo.iloc[3::2],Samplesinfo.iloc[4::2],label=Samplesinfo.iloc[0])
            for key in self.samplelist.keys():
                plt.text(self.samplelist[key][5],self.samplelist[key][6],self.samplelist[key][0])
            plt.show()
        except Exception as e :
             print(str(e))
    def plotshow(self):
        
       
        try:
            ms_path=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
            log_path=QtCore.QDir.toNativeSeparators(self.dataEntry2.text())
            
            if self.instrument==3:        
                datareader=DataReader(ms_path,filetype='csv',instrument='Agilent1')
                ms_data=datareader.read(ms_path,filetype='csv',instrument='Agilent1')
            elif self.instrument==2:
                datareader=DataReader(ms_path,filetype='csv',instrument='Agilent2')
                ms_data=datareader.read(ms_path,filetype='csv',instrument='Agilent2')
            elif self.instrument==13:
                datareader=DataReader(ms_path,filetype='csv',instrument='Thermo')
                ms_data=datareader.read(ms_path,filetype='csv',instrument='Thermo')
        except Exception as e:
            print(str(e))
            
        
        try:
            a=Iolite(log_path,self.laser)
            iolite=a.read_iolite(log_path)
            names=a.names_from_iolite(self.laser)
            
            times=a.on_and_off_times()
            b=Selector(datareader,times=times,s=eval(self.start_timee.text()))
            starts,ends=b.create_selector_iolite()
            
            fig, ax = plt.subplots()
    
            ax.cla()
            ax.clear()
           
            self.time=ms_data.index
            
            plt.plot(self.time,ms_data)
            
    
            ax.set_yscale('log')
    
            if starts and ends:
                # create lines for start and end of each ablation
                for i in range(0, len(starts)):
                    ax.axvline(x=self.time[starts[i]],
                               color='red', linewidth=2)
                for i in range(0, len(ends)):
                    ax.axvline(x=self.time[ends[i]],
                               color='blue', linewidth=2)
    
    
            plt.show()
        except Exception as e:
            print(str(e))
    def convertCMD(self):
        try:
            ms_path=QtCore.QDir.toNativeSeparators(self.dataEntry1.text())
            log_path=QtCore.QDir.toNativeSeparators(self.dataEntry2.text())
            outputpath=QtCore.QDir.toNativeSeparators(self.dataEntry3.text())
            diff_sec=eval(self.start_timee.text())
            width=eval(self.extrowe.text())
            Dir_file=self.Dir_filee.text()
            
            if self.instrument==3:
                convert_segmentation=Iolite_to_Aglient(self.laser,log_path=log_path,ms_path=ms_path,
                         outputpath=outputpath,diff_sec=diff_sec,width=width,Dir_file=Dir_file,instrument='Agilent1')
            elif self.instrument==2:
                convert_segmentation=Iolite_to_Aglient(self.laser,log_path=log_path,ms_path=ms_path,
                         outputpath=outputpath,diff_sec=diff_sec,width=width,Dir_file=Dir_file,instrument='Agilent2')
            elif self.instrument==13:
                convert_segmentation=Iolite_to_Aglient(self.laser,log_path=log_path,ms_path=ms_path,
                         outputpath=outputpath,diff_sec=diff_sec,width=width,Dir_file=Dir_file,instrument='Thermo')
            self.samplelist=convert_segmentation.creat_Aglient()
            #except Exception as e:
                #print(str(e))
            msg_box_c=QMessageBox(QMessageBox.Information,'Information','Convertion successful!')
            msg_box_c.exec_()
        except Exception as e:
            print(str(e))
        
    def get_dirctory(self,entry):
        pwd=os.getcwd()
        self.dir_path=QFileDialog.getExistingDirectory(self,'OPen Path',pwd)
        if self.dir_path:
            entry.setText(self.dir_path)
    def get_file(self, entry):            
        pwd=os.getcwd()
        filename, filters = QFileDialog.getOpenFileName(self,'Choose *.csv',pwd,'Map file (*.csv)::' )
        if filename:
            entry.setText(filename) 


if __name__=='__main__':
    app=QApplication(sys.argv)
    #app.setWindowIcon(QIcon('./images/home.png'))    
    main=MyApplication()
    main.show()
    sys.exit(app.exec_())

