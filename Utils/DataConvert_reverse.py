import pandas as pd
import math
import os

needProcessParticularA = True  #是否要处理特殊的A

def DeleteNullColumsAndNoCheck(data):
    firstLine = data.columns
    delindex = []
    for index, i in enumerate(firstLine):
        if('Unnamed' in i):
            delindex.append(i)
    data = data.drop(delindex,axis=1)

    idname = data['Squib']
    dellineindex = []
    for index, i in enumerate(idname):
        if(pd.isna(i)):
            continue
        if('Reference reaction' in i):
            dellineindex.append(index)
            dellineindex.append(index-1)
    data = data.drop(dellineindex,axis=0)        
    return data

def ConvertDataReverse(excelname, outname):
    excelbaesname = os.path.basename(excelname)
    #reverse reaction
    reaction = excelbaesname.split('.')[0].split("_")[2]+"<=>"+excelbaesname.split('.')[0].split("_")[1]
    print(reaction)
    data = pd.read_excel(excelname)
    #删除空行或者带Reference reaction的反应(这种反应是nocheck的)
    data = DeleteNullColumsAndNoCheck(data)
    datasave = pd.DataFrame(columns=data.columns)

    #只要有用的数据
    for index, row in data.iterrows():
        if((pd.isna(row['Squib']) or pd.isna(row['Temp [K]'])) and (pd.isna(row['Plot']))):
            continue
        #print(row['Squib'])
        #print(row['A'])

        if(pd.notna(row['Squib']) and pd.isna(row['A'])):
            continue
        
        print(row['Squib'])
        datasave.loc[datasave.shape[0]] = row

    newData = pd.DataFrame(columns=['Plot','Squib','Use','low','high','A','n','E/R','k0','Order','11','12','13','14','15','16','17','18','19','20','21'])

    #print(datasave)

    #开始逐行计算各个值
    CurrentPlot = ''
    PlotHaveChange = False
    for index, row in datasave.iterrows():
        if(pd.isna(row['Temp [K]'])):#说明这行是Plot = Review/Experiment这样表头
            CurrentPlot = row['Plot']
            PlotHaveChange = True #Plot已经更改
            #print("current plot is ",CurrentPlot)
        else:
            Plot = ''
            #print(row['Squib'])

            #判断使用Review/Experiment还是使用0
            if(PlotHaveChange):
                Plot = CurrentPlot
                PlotHaveChange = False
            else:
                Plot = '0'

            N = row['n']
            if(pd.isna(N)):
                N = 0
            A_temp = row['A']
            if (type(A_temp) == str):
                if(needProcessParticularA == False):
                    continue
                A_temp = A_temp[1:]
                A_temp = float(A_temp)
            A = A_temp*6.023E+23*((1/298)**float(N))
            if('Ea [J/mole]' in datasave.columns):#reverse不需要生成EAR，只生成EA
                ER = (row['Ea [J/mole]']/1000)*238.9
            else:
                ER = (row['Ea [kJ/mole]'])*238.9

            if(pd.isna(ER)):
               ER = 0 
            LOW= 0
            HIGH = 0
            K0 = 0
            interval = 0
            if(pd.isna(row['interval'])):
                interval = 0
            else:
                interval = int(row['interval'])

            Temp = row['Temp [K]']
            if(type(Temp) == str):
                splitTemp = Temp.split('-')
                LOW = int(splitTemp[0])
                HIGH = int(splitTemp[1])
                if((LOW <= 289) and (HIGH >= 298)):
                    K0 = A*math.exp((-ER)/(298*1.987))
            else:
                LOW = int(Temp)
                if(LOW == 298):
                    K0 = A*math.exp((-ER)/(298*1.987))
            oneline =  [Plot,row['Squib'],'1',LOW, HIGH, A, N,ER, K0, row['Order'],0,0,0,0,0,1,0,0,0,0,interval]
            newData.loc[newData.shape[0]] = oneline

    fd = open(outname,'w+')
    for index, row in newData.iterrows():
        oneLine = reaction+' '+"{:.3e}".format(row['A'])+\
            " "\
            +"{:<5s}".format(str(row['n']))+\
            ' '+\
            "{:<5s}".format(str(round(row['E/R'],2)))+\
            '  '+\
            "{:<20s}".format('reverse')+\
            "{:<30s}".format(row['Squib'])+\
            "{:<5s}".format(str(row['Use']))+\
            "{:<10s}".format(str(row['low']))+\
            "{:<10s}".format(str(row['high']))+\
            "{:.5e}".format(row['A'])+\
            "    "+\
            "{:<10s}".format(str(row['n']))+\
            "{:<20s}".format(str(round(float(row['E/R']),4)))+\
            "{:.5e}".format(row['k0'])+\
            "    "+\
            "{:<5s}".format(str(int(row['Order'])))+\
            "{:<5s}".format(str(row['11']))+\
            "{:<5s}".format(str(row['12']))+\
            "{:<5s}".format(str(row['13']))+\
            "{:<5s}".format(str(row['14']))+\
            "{:<5s}".format(str(row['15']))+\
            "{:<5s}".format(str(row['16']))+\
            "{:<5s}".format(str(row['17']))+\
            "{:<5s}".format(str(row['18']))+\
            "{:<5s}".format(str(row['19']))+\
            "{:<5s}".format(str(row['20']))+\
            "{:<5s}".format(str(row['21']))+'\n'

        fd.write(oneLine)
    fd.close() 

