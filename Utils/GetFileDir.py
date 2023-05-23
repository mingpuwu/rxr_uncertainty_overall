import os

#num反应序列，type:forwardexcel,reverseexcel,forwardtxt,reversetxt,dir
#根据反应序列返回文件名路径，比如返回excel文件名称，GetFileDir(99,'excel')
def GetFileDir(num,type):
    dirlist = os.listdir('Data')
    dir = ''
    for i in dirlist:
        #print(i)
        n = i.split('_')[0]
        #print(n)
        if(n == str(num)):
            dir = i
            break
    #print(dir)
    if(type == 'forwardexcel'):
        return os.path.join('Data',dir,dir+'.xlsx') 
    elif(type == 'reverseexcel'):
        return os.path.join('Data',dir,dir+'_reverse'+'.xlsx') 
    elif(type == 'forwardtxt'):
        return os.path.join('Data',dir,dir+'.txt') 
    elif(type == 'reversetxt'):
        return os.path.join('Data',dir,dir+'_reverse'+'.txt')
    elif(type == 'dir'):
        return os.path.join('Data',dir) 
    
if __name__ == '__main__':
    name = GetFileDir(99,'excel')
    print(name)