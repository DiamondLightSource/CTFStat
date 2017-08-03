vers = '0.2'
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import logging

global FigNum # number of figures


FigNum = 0
#------------------------------Inputs -----------------------------------------#
class Inp(object):
    _registry = []
    def __init__(self, flag, value, req):
        self._registry.append(self)
        self.flag = flag
        self.value = value
        self.req = req


def create_inp(flag, value, req):
    errmsg = '''\n    CFStat is a quality assurance program. It provide statistics about quality of micrographs during acquisition.
\nUsage: CFTStat.py  [OPTION]... star_file  (E.g.: CFTStat.py --star <micrographs_all_gctf.star>
\nOptions
\n--star, --verbose               specify input file type (star file)\n'''
    Argument = Inp(flag, value, req)
    
        
    if Argument.value == False:
        if Argument.flag in sys.argv:
            Argument.value = True
        else:
            Argument.value = False
    return Argument.value


#*********Parse the star file  *********#
def starfile_Parser(f):
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    for i in alldata:
        if '#' in i:
            labelsdic[i.split('#')[0]] = int(i.split('#')[1])-1
        if len(i.split()) > 3:
            data.append(i.split())
        if len(i.split()) < 3:
            header.append(i.strip("\n"))
    return(labelsdic,header,data)

#---- read the CTF information from the inputed star file 
def plot_ctfstats(alldata,mresvals,mastigvals):
    global FigNum
    global dpi
    data = []
    for i in alldata:
        if len(i.split()) > 3:
            data.append(i.split())
        if '_rlnDefocusU' in i:
            dfucol = int(i.split('#')[-1])-1
        if '_rlnDefocusV' in i:
            dfvcol = int(i.split('#')[-1])-1
        if '_rlnDefocusAngle' in i:
            dfacol = int(i.split('#')[-1])-1
        if '_rlnMicrographName' in i:
            namecol = int(i.split('#')[-1])-1
        if '_rlnFinalResolution' in i:
            rescol = int(i.split('#')[-1])-1
        if 'rlnCtfMaxResolution' in i:
            rescol = int(i.split('#')[-1])-1
            
    v,u,a,names,res = [],[],[],[],[]
    for i in data:
        u.append(float(i[dfucol]))
        v.append(float(i[dfvcol]))
        a.append(float(i[dfacol]))
        names.append(i[namecol])
        res.append(float(i[rescol]))
    amin = min(a)
    amax = max(a)
    maxfactor = 1/amax
    

    micsdic = {}     #-- images dictionary
    count = 0
    for i in names:
        micsdic[i] = (u[count],v[count],a[count],res[count]) 
        count+=1

    scaleda = []       #---   plot defocus
    for i in a:
        scaleda.append(i*maxfactor)
    

    plt.subplot2grid((3,2), (0, 0), colspan=2)       # --  dfplot
    plt.scatter(u, v, s=10, c=scaleda)        
    plt.xlabel('DefocusU',fontsize=10)
    plt.ylabel('DefocusV',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)
    

    astig = []       #--  astigmatism plots
    count = 0
    for i in u:
        astig.append(abs(i - v[count]))
        count +=1    
    plt.subplot2grid((3,2), (1, 0))
    n, bins, patches = plt.hist(astig, 100, facecolor='blue', alpha=0.75)
    plt.xlabel('Astigmatism',fontsize=10)
    plt.ylabel('Micrographs',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)
    
    dfumean = np.mean(u)
    dfvmean = np.mean(v)

    plt.subplot2grid((3,2), (2, 0))      #--  plot  resolution
    n, bins, patches = plt.hist(res, 100, facecolor='green', alpha=0.75)
    plt.xlabel('Resolution',fontsize=10)
    plt.ylabel('Micrographs',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=8)

    plt.subplot2grid((3,2), (1, 1))
    plt.plot(mresvals)

    plt.tick_params(axis='y', which='major', labelsize=8)
    plt.xticks([])
    
    nums = range(len(mresvals))
    fit = np.polyfit(nums,mresvals,1)
    x = np.arange(min(nums),max(nums))
    eq = '{0}*x+{1}'.format(fit[0],fit[1])
    y = eval(eq)
    plt.plot(x,y,color='Red')
    plt.xlabel('Res/Time {0}%'.format(round(fit[0]*100,7)),fontsize=10)
    plt.ylabel('Resolution % mean',fontsize=10)
    
    plt.subplot2grid((3,2), (2, 1))
    plt.plot(mastigvals)
    plt.tick_params(axis='y', which='major', labelsize=8)
    plt.xticks([])

    fit = np.polyfit(nums,mastigvals,1)
    eq = '{0}*x+{1}'.format(fit[0],fit[1])
    y = eval(eq)
    plt.plot(x,y,color='Red')
    plt.xlabel('Astig/Time {0}%'.format(round(fit[0]*100,7)),fontsize=10)
    plt.ylabel('Astig % mean',fontsize=10)
    
    plt.tight_layout()
    plt.savefig('CTF_QA_{0}.png'.format(FigNum))
    if filter == True:
        plt.show(block=False)
    FigNum +=1 
    return (micsdic)
    
#----  total running -------------------------------
def running_total(dic):
    dickeys = dic.keys()
    dickeys.sort()
    means = []
    running_total = 0
    n = 1
    for dicline in dickeys:
        for val in dic[dicline]:
            running_total = running_total + float(val)
            n+=1
    mean = running_total/n

    for dicline in dickeys:
        for val in dic[dicline]:
            means.append(float(val)/mean)    
    return(means)
#----------------------------------------------------
# this is static code need to be changed. it was written like this to deal with the error message that shows when inputting arguments
thefile = 'micrographs_all_gctf.star'
alldata = open(thefile,'r').readlines()
filter = create_inp('--f',False,False)
print("**** CTFStat v{0} ****".format(vers))


#-----total computations---------------
def read_time_stamp(file):
    (star_labels,star_header,star_data) = starfile_Parser(file)
    astigdatesdic = {}
    resdatesdic = {}
    for i in star_data:
        date = ''.join(i[star_labels['_rlnMicrographName ']].split('_')[-3:-1])
        if date not in resdatesdic:
           resdatesdic[date] = [i[star_labels['_rlnFinalResolution ']]]
        else:
            resdatesdic[date].append(i[star_labels['_rlnFinalResolution ']])
        
        astig = abs(float(i[star_labels['_rlnDefocusU ']])-float(i[star_labels['_rlnDefocusV ']]))
        if date not in astigdatesdic:
           astigdatesdic[date] = [astig]
        else:
            astigdatesdic[date].append(astig)
    return(astigdatesdic,resdatesdic)
#------------------------------------------------------------

(astigdatesdic,resdatesdic) = read_time_stamp(thefile) 

astigmeans = running_total(astigdatesdic)
resmeans = running_total(resdatesdic)

# standard calculations
micrographs = plot_ctfstats(alldata,resmeans,astigmeans)

# filtering if necesary
if filter == False:
    sys.exit('Done:  CTF_QA.png Saved')
