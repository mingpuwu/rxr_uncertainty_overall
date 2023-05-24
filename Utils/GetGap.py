import numpy as np

def gap(templow, temphigh, interval):
    lowinter = int(templow/interval)
    lowresidue = templow%interval
    highinter = int(temphigh/interval)
    highresidue = temphigh%interval
    num = (temphigh-templow)/interval+1
    gapdata = []
    if(lowresidue == 0 and highresidue == 0): #都可以整除
        gapdata = np.linspace(templow, temphigh, int(num)).tolist()
    elif(lowresidue != 0 and highresidue == 0):
        fisrtinter = (lowinter+1)*interval
        gapdata.append(templow)
        gapdata.extend(np.linspace(fisrtinter, temphigh, int((temphigh-fisrtinter)/interval)+1).tolist())
    elif(lowresidue == 0 and highresidue != 0):
        endinter = (highinter)*interval
        gapdata.extend(np.linspace(templow, endinter, int((endinter-templow)/interval)+1).tolist())
        gapdata.append(temphigh)
    else:
        gapdata.append(templow)
        fisrtinter = (lowinter+1)*interval
        endinter = (highinter)*interval
        gapdata.extend(np.linspace(fisrtinter, endinter, int((endinter-fisrtinter)/interval)+1).tolist())
        gapdata.append(temphigh)
    
    return gapdata

def GetGap(templow, temphigh, interval = 0):
    #如果给定了这个值的话，按照给定的值来
    if(interval != 0):
        #print("--------------interval is :",interval)
        return gap(templow, temphigh, interval)
    
    #如果没有给inteval,自己判定
    if(temphigh - templow < 500):
        interval = 25
    elif((temphigh - templow >= 500) and (temphigh - templow <=1000)):
        interval = 50
    else:
        interval = 100
    #print("--------------interval is :",interval)
    return gap(templow, temphigh, interval)
