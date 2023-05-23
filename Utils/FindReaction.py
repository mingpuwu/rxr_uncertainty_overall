import os

ttfilename = 'Data/TTMechanism.txt'
henrryname = 'Data/C3MechV3.3.MECH'

ttdata = []
henrrydata = []

#needtype:'LOW'/''
def FindReaction(excelname,needtype = '',ttname = '', henryname = ''):
    ttreaction = ''
    excelname = os.path.basename(excelname)
    name = excelname.split('_')
    ttreaction = name[1]
    ttreaction = ttreaction + '<=>'
    ttreaction = ttreaction + name[2].split('.')[0]
    if(len(ttname) != 0 ):
        ttreaction = ttname

    print("want to find ",ttreaction)
    fd = open(ttfilename, 'r',encoding='gbk')
    alldata = fd.readlines()

    for index,i in enumerate(alldata):
        if '=' in i:
            r = i.split()[0]
            if(ttreaction == r):
                print("TT find it ",i)
                if(needtype == 'LOW'): #果如是low的需要使用公式的下一行数据
                    splitdata = alldata[index+1].split()
                    ttdata.append(splitdata[2])
                    ttdata.append(splitdata[3])
                    ttdata.append(splitdata[4])
                else:
                    splitdata = i.split()
                    ttdata.append(splitdata[1])
                    ttdata.append(splitdata[2])
                    ttdata.append(splitdata[3])
    fd.close()

    if(len(ttdata) == 0):
        print("error TT reaction not find ",ttreaction)

    #thrid body的形式带(+M)

    Preequal_oneItem = excelname.split('_')[1].split('+')[0]
    Preequal_twoItem = excelname.split('_')[1].split('+')[1]
    Afterequal_oneItem = excelname.split('_')[2].split('.')[0].split('+')[0]
    Afterequal_twoItem = excelname.split('_')[2].split('.')[0].split('+')[1]
    print(Preequal_oneItem, Preequal_twoItem, Afterequal_oneItem, Afterequal_twoItem)
    PreoneReaction = Preequal_oneItem+"+"+Preequal_twoItem
    PretwoReaction = Preequal_twoItem+"+"+Preequal_oneItem
    AfteroneReaction = Afterequal_oneItem+'+'+Afterequal_twoItem
    AftertwoReaction = Afterequal_twoItem+'+'+Afterequal_oneItem
    #print(oneReaction)
    #print(twoReaction)
    fd = open(henrryname, 'r',encoding='gbk')
    alldata = fd.readlines()

    #使用指定的名字查找
    if(len(henryname) != 0):
        for index,i in enumerate(alldata):
            if('=' in i):
                r = i.split()[0]
                if(henryname == r):
                    print("hen find it ",r)
                    if(needtype == 'LOW'):#需要用下一行的数据
                        splitdata = alldata[index+1].split()
                        henrrydata.append(splitdata[1])
                        henrrydata.append(splitdata[2])
                        henrrydata.append(splitdata[3].replace('/',''))
                    else:
                        splitdata = alldata[index]
                        henrrydata.append(splitdata[1])
                        henrrydata.append(splitdata[2])
                        henrrydata.append(splitdata[3])
    else:
        for i in alldata:
            prefirst = ''
            presecond = ''

            if('=' in i):
                r = i.split()[0]
                prefirst = r.split('=')[0]
                presecond = r.split('=')[1]
            if((PreoneReaction == prefirst ) or (PretwoReaction == prefirst)):#拼出的公式是不是在前面
                #print("first",i)
                if((AfteroneReaction == presecond) or (AftertwoReaction == presecond)):
                    print("hen find it ",i)
                    splitdata = i.split()
                    henrrydata.append(splitdata[1])
                    henrrydata.append(splitdata[2])
                    henrrydata.append(splitdata[3])
            elif ((AfteroneReaction == prefirst) or (AftertwoReaction == prefirst)):
                #print("second",i)
                if((PreoneReaction == presecond) or (PretwoReaction == presecond)):
                    print("hen find it ",i)
                    splitdata = i.split()
                    henrrydata.append(splitdata[1])
                    henrrydata.append(splitdata[2])
                    henrrydata.append(splitdata[3])
    fd.close()

    if(len(henrrydata) == 0):
        print("error hen reaction not find ",ttreaction)

    return ttdata, henrrydata