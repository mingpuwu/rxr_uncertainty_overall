import os
from Utils.GetGap import GetGap
import numpy as np
import matplotlib.pyplot as plt

markes = ['o', 's', '^', 'p', '^', 'v', 'p', 'd', 'h', 
          '2', '8','g.','b^','yo','bp','gs','kx','rv',
          'bd','y1','|','d','D','x','+','H','h','*','p','s']

def Plot(filename, xlim=5,ylim_max=14.5, ylim_min=12, ncol=1):
    fd = open(filename,encoding="utf-8")

    alllines = fd.readlines()

    fitdata = []
    compareOccurrence = False
    comparedata = []
    comparename = []
    otherdata = []
    othername = []
    allname = []
    CurrentName = ''
    for i in alllines:
        OneLine = i.split()

        if(OneLine[2] == '0'):#unused跳过
            continue
        
        name = OneLine[1]
        A = float(OneLine[5])
        N = float(OneLine[6])
        EAR = float(OneLine[7])
        M = float(OneLine[15])
        lowTemp = int(OneLine[3])
        HighTemp = 0
        if(OneLine[4] == '0'):
            HighTemp = lowTemp
        else:
            HighTemp = int(OneLine[4])

        interval = int(OneLine[20])
        if(len(CurrentName) == 0):
            CurrentName = OneLine[0]
        
        if(CurrentName != OneLine[0] and OneLine[0] != '0'):
            CurrentName = OneLine[0]
        
        allname.append(CurrentName+';'+OneLine[1])
        #allname.append(name)

        if(OneLine[1] == 'fit'):
            t = GetGap(300,2500,100)
            t = np.array(t)
            fitdata.append(t)
            fitdata.append((np.log(A)+N*np.log(fitdata[0])+(-EAR/fitdata[0])+np.log(1/M))/np.log(10))
            continue

        if(OneLine[0] == 'Compare'):
            compareOccurrence = True
        
        if(compareOccurrence):
            oneitem = []
            t = GetGap(lowTemp,HighTemp,interval)
            t = np.array(t)
            oneitem.append(t)
            oneitem.append((np.log(A)+N*np.log(oneitem[0])+(-EAR/oneitem[0])+np.log(1/M))/np.log(10))
            comparedata.append(oneitem)
            comparename.append(name)
            continue

        oneitem = []
        t = GetGap(lowTemp,HighTemp,interval)
        t = np.array(t)
        oneitem.append(t)
        #print(oneitem[0])
        oneitem.append((np.log(A)+N*np.log(oneitem[0])+(-EAR/oneitem[0])+np.log(1/M))/np.log(10))
        otherdata.append(oneitem)
        othername.append(name)

    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内

    for index, plotone in enumerate(otherdata):
        plt.plot(1000./plotone[0],plotone[1],markes[index%len(markes)],markersize=5)

    plt.plot(1000./fitdata[0], fitdata[1], '-',color="red", linewidth=1.0)

    color=["aqua","m","g","b"]
    for index, plotone in enumerate(comparedata):
        plt.plot(1000./plotone[0],plotone[1],linewidth=1.5, color=color[index])

    
    #print(allname)
    plt.legend(allname,loc=2, bbox_to_anchor=(1.01,1.0),borderaxespad = 0.,  ncol = ncol)
    #plt.legend(allname)
    plt.xlim((0, xlim))
    plt.ylim((ylim_min, ylim_max))
    plt.xlabel("1000/T",fontsize=14)
    plt.ylabel("log$_{10}$(k)",fontsize=14)
    plt.title("Arrhenius plot")
    plt.grid()#添加网格

    plt.show()

if __name__ == '__main__':
    Plot("Data/99_CH3+HO2_CH4+O2/99_CH3+HO2_CH4+O2_all.txt")

