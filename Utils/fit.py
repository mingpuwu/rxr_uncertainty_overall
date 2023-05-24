import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from Utils.GetGap import GetGap

#fileName = "UQ_143_all.txt"
markes = ['o', 's', '^', 'p', '^', 'v', 'p', 'd', 'h', 
          '2', '8','g.','b^','yo','bp','gs','kx','rv',
          'bd','y1','-','|','d','D','x','+','H','h','*','p','s']

def func(T, A, n, EAR):
    return (np.log(A)+n*np.log(T)+(-EAR/T))/np.log(10)

def FitHandler(T,K):
    #非线性最小二乘法拟合
    T = np.array(T)
    K = np.array(K)
    bounds = [[0,-2,-np.inf],[np.inf,2,np.inf]]
    popt, pcov = curve_fit(func, T, K ,maxfev=100000,bounds = bounds)
    #获取popt里面是拟合系数
    print(popt)
    A = popt[0] 
    n = popt[1]
    EAR = popt[2]
    
    xvals = np.linspace(300,2500,23)
    xvals = np.array(xvals)
    yvals = func(xvals, A, n, EAR) #拟合y值

    mean = np.mean(K)  # 1.y mean
    ss_tot = np.sum((K - mean) ** 2)  # 2.total sum of squares
    ss_res = np.sum((K - func(T, *popt)) ** 2)  # 3.residual sum of squares
    r_squared = 1 - (ss_res / ss_tot)  # 4.r squared
    r_squared2 = r2_score(K, func(T, *popt))
    std = np.std(K - func(T, *popt))
    print("sum_error**2",ss_res)
    print('popt:', popt)
    print('系数A:', A)
    print('系数n:', n)
    print('系数EA:', EAR)
    print('系数协方差:', pcov)
    print('参数标准差:',np.sqrt(np.diag(pcov)))
    print('RSE:',r_squared)
    print('RSE2:',r_squared2)
    print('std:',std)
    print('系数yvals:', yvals)
    #绘图
    T_1000 = 1000/T
    xvals = 1000/xvals
    plot1 = plt.plot(T_1000, K, '.',label='original values')
    plot2 = plt.plot(xvals, yvals, '-r',label='polyfit values')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc = 4) #指定legend的位置右下角
    plt.title('curve_fit')
    plt.show()
    return A,n,EAR,xvals,yvals

def Fit(fileName, reactionnum):
    CompareDataCount = 0
    Comparex = []
    Comparey = []
    fd = open(fileName)
    allData = fd.readlines()
    Name = []
    TempLow = []
    TempHigh = []
    Used = []
    A = []
    N = []
    EAR = []
    M = []
    interval = []

    for i in allData:
        OneLine = i.split()

        #fit不参与拟合
        if(OneLine[1] == 'fit'):
            continue

        #当出现Comare时后面的都不需要fit,先记录下有几Compare个
        if(OneLine[0] == 'Compare'):
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            CompareDataCount = 1

        #后面都是compare,都需要+1
        if(CompareDataCount):
            CompareDataCount = CompareDataCount + 1

        if(OneLine[2] == '0'):#unused跳过
            continue

        #print(OneLine[1])
        Name.append(OneLine[1])
        TempLow.append(int(OneLine[3]))
        if(OneLine[4] == '0'):#没有High temp,用Low填写
            TempHigh.append(int(OneLine[3]))
        else:
            TempHigh.append(int(OneLine[4]))
        A.append(float(OneLine[5]))
        N.append(float(OneLine[6]))
        EAR.append(float(OneLine[7]))
        M.append(float(OneLine[15]))
        interval.append(int(OneLine[20]))

    #print(M)

    rangeGap = []
    allTemp = []

    intervafd = open('interval','w+')
    
    for i in range(len(TempLow)):
        OneGap = GetGap(TempLow[i],TempHigh[i],interval[i])
        OneGap = np.array(OneGap)

        #如果次时遍历到了compare,需要把他们单独存放，防止参与拟合
        if(i > len(TempLow)-CompareDataCount):
            Comparex.append(OneGap)
        else:
            allTemp.extend(OneGap)

        #rangeGap里面存放的是包含compare的数据
        rangeGap.append(OneGap)
    intervafd.write(str(Comparex))
    intervafd.close()
    #print(rangeGap)
    K = []
    koriginal=[]
    allK = []
    A = np.array(A)
    N = np.array(N)
    EAR = np.array(EAR)
    rangeGap = np.array(rangeGap)

    # fd = open('review','w+')

    for i in range(len(TempLow)):
        OneK = (np.log(A[i])+N[i]*np.log(rangeGap[i])+(-EAR[i]/rangeGap[i])+np.log(1/M[i]))/np.log(10)
        oriK=A[i]*rangeGap[i]**N[i]*np.exp(-EAR[i]/rangeGap[i])*(1/M[i])

        #如果遍历到了compare，需要把他们单独存放，防止参与拟合
        if(i > len(TempLow)-CompareDataCount):
            Comparey.append(OneK)
        else:
            allK.extend(OneK)
            koriginal.extend(oriK)
        
        K.append(OneK)

    T_K=pd.DataFrame()
    T_K["allTemp"]=allTemp
    T_K["koriginal"]=koriginal
    T_K.sort_values(by="allTemp",inplace=True,ascending=True)
    #print(T_K)
    T_K.to_excel(str(reactionnum)+".xlsx")
    # fd.write(str(K))
    # fd.close()
    
    A,n,EAR,xvals,yvals = FitHandler(allTemp, allK)
    Name.append('fit')
    for index, plotindex in enumerate(K):
        if(index > len(TempLow)-CompareDataCount):
            continue
        plt.plot(1000./rangeGap[index],plotindex,markes[index%len(markes)])

    plt.plot(xvals,yvals,color="red",linewidth=2.0)

    for i in range(CompareDataCount-1):
        plt.plot(1000./Comparex[i],Comparey[i],color="red",linewidth=2.0)

    plt.legend(Name,loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)

    plt.xlabel("1000/T")
    plt.ylabel("log_10(k)")
    #plt.legend('fit')
    plt.show()
    return A,round(n,2),round(EAR,6)

if __name__ == '__main__':
    Fit('Data/99_CH3+HO2_CH4+O2/99_CH3+HO2_CH4+O2.txt')