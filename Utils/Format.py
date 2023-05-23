def Format(Plot, name, lowt, hight, A, n, EAR,k0):
    oneLine = "{:<20s}".format(Plot) + \
    "{:<30s}".format(name)+\
    "{:<5s}".format('1')+\
    "{:<10s}".format(str(lowt))+\
    "{:<10s}".format(str(hight))+\
    "{:.5e}".format(A)+\
    "    "+\
    "{:<10s}".format(str(round(n,2)))+\
    "{:<20s}".format(str(round(EAR,4)))+\
    "{:.5e}".format(k0)+\
    "    "+\
    "{:<5s}".format(str(2))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(1))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+\
    "{:<5s}".format(str(0))+'\n'

    return oneLine