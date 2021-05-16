import scipy.optimize # for fitter
import sys
from math import *

myp=2.7928
mp=0.938272

def loaddata(filename):
    data=[]
    maxnorm=0;
    for l in open(filename):
        try:
            values=l.split()
            energy=float(values[0])
            q2=float(values[3])
            cs=float(values[4])
            dcs=float(values[5])
            norms=[int (i) for i in values[10].split(":")]
            for n in norms:
                if maxnorm<n:
                    maxnorm=n
            data.append([energy/1000,q2,cs,dcs,norms])
        except:
            continue
    return data,maxnorm

def poly(q2,params): # this calculates 1+param[0]*q^2+param[1]*q^4...
    n=len(params)
    value=1;
    q=q2
    for p in params:
        value+=p*q
        q=q*q2
    return value
        
def tau(q2):
    return q2/4/mp**2;

def theta(E,q2):
    return 2*acos(sqrt((-4*E*E*mp+2*E*q2+mp*q2)/(-4*E*E*mp+2*E*q2))
)

def dip(q2):
    return 1/(1+q2/0.71)**2

def cs(E,q2,parge,pargm):
    ge2=poly(q2,parge)**2
    gm2=(myp*poly(q2,pargm))**2
    ta=tau(q2)
    th=theta(E,q2)
    cross= (ge2+ta*gm2)/(1+ta)+2*ta*gm2*tan(th/2)**2
    ge2=dip(q2)**2
    gm2=ge2*myp**2
    dipcross=(ge2+ta*gm2)/(1+ta)+2*ta*gm2*tan(th/2)**2
    return cross/dipcross


    

data,numnorm=loaddata(sys.argv[1]);
iter=0
print "Loaded",len(data),"data points,",numnorm,"normalization parameters"


def fitfunc(params):
    global iter
    chi2=0
    res=[]
    for d in data:
        model=cs(d[0],d[1],params[0:10],params[10:20]) # calculate model value

        norm=1
        for np in d[4]:
            norm*=params[np-1 +20] #calculate norm product
            
        datacs=d[2]*norm
        datacse=d[3]*norm

        err=(datacs-model)/datacse
        chi2+=err*err
        res.append(err)
    if iter %10==0: print iter, chi2 , params
    iter+=1
    return res;

# this a seed to increase the fit convergence speeed
p=[-3.368,14.56,-88.19,453.62,-1638.79,+3980.71,-6312.63,+6222.36,-3443.22,+814.41,-2.59,1.02,23.49,-93.0372,+140.79,-0.365,-305.67,+444.62,-273.66,+64.58]
erg=scipy.optimize.leastsq(fitfunc,p+[1]*numnorm,full_output=1)

x=erg[0]
print "Final parameters: ",x

of=open(sys.argv[2],"w")

for qi in range (100):
    q2=qi*0.01
    print >>of,q2,poly(q2,x[0:10]),poly(q2,x[10:20])
                    
of.close()
