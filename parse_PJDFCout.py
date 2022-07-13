#python3 parser of the output from WISE-steep-MCMC-PJDFC
#run this in the directory made for the object

def jd2ymd(jd):
    """
        Returns a date corresponding to the given Julian day number.
    """
    jd2 = float(jd) + 0.5
    dfrac, dint = modf(jd2)
    x2 = floor(jd - 1721119.5)
    c2 = floor((4*x2 + 3)/146097)
    x1 = x2 - floor(146097*c2/4)
    c1 = floor((100*x1 + 99)/36525)
    x0 = x1 - floor(36525*c1/100)
    year = 100*c2 + c1
    month = floor((5*x0 + 461)/153)
    day = x0 - floor((153*month - 457)/5) + 1
    if month > 12:
        month = month - 12
        year = year + 1    
    return year, month, day + dfrac

import sys
import os
import matplotlib.pyplot as plt
import numpy
from math import modf,floor

if sys.version[0] != '3':
    print("Please run with python3")
    exit()

#cp, not mv, so that there is a record of the actual outputs
    
os.system("/bin/cp fort.21 PJDFC.out")
#fort.22 is the postscript format of the SED
os.system("/bin/cp fort.22 SED_data.out")


#^*^ parse and plot SED_data.out here 
#best fit solution used for SED
inf=open("SED_data.out")
line1=inf.readline().rstrip()
line2=inf.readline().rstrip()

best={}
dat1=line1.split()
best['pole_ra']=numpy.degrees(float(dat1[1]))
best['pole_dec']=numpy.degrees(float(dat1[2]))
best['pV_median']=numpy.e**(float(dat1[3]))
best['period']=numpy.e**(float(dat1[4]))
best['Gamma']=numpy.e**(float(dat1[5]))
best['pjdfc']=(float(dat1[6])) #period * thermal inertia Gamma * diameter * crater frac * color
x=float(dat1[7])
best['crater_frac']=numpy.e**(x)/(1 + numpy.e**x)
best['ir_albedo_ratio']=numpy.e**(float(dat1[8]))
best['chisq']=float(dat1[10])
best['penalties']=float(dat1[9])-best['chisq']

dat2=line2.split()
best['pV']=float(dat2[1])
best['diameter']=float(dat2[2])
best['Theta']=float(dat2[3])

foo=inf.readline()
epoch_condition={}
epoch_condition['Delta']=[]
epoch_condition['R_helio']=[]
epoch_condition['phase']=[]
epoch_condition['subSunLat']=[]
epoch_condition['subEarthLat']=[]

nepoch=0
nexte=inf.readline()
while nexte[0:4]!='/wvl':
    nepoch+=1
    edat=nexte.split()
    epoch_condition['Delta'].append(float(edat[1]))
    epoch_condition['R_helio'].append(float(edat[2]))
    epoch_condition['phase'].append(float(edat[3]))
    epoch_condition['subSunLat'].append(float(edat[4]))
    epoch_condition['subEarthLat'].append(float(edat[5]))
    nexte=inf.readline()

wave=[]
nextl=inf.readline()
while 'def' not in nextl:
    wdat=nextl.strip().split()
    for w in wdat:
        wave.append(float(w))
    nextl=inf.readline()

color=[]
ecount=0
esed=[]
obsdatax=[]
obsdatay=[]
obsdatayerr=[]
while True:
    line=inf.readline()
    if line=="":
        break
    
    if 'SRGB' in line:
        ecount+=1
        (rr,gg,bb,foo)=line.rstrip().split()
        rout=hex(int(float(rr)*255))[2:]
        gout=hex(int(float(gg)*255))[2:]
        bout=hex(int(float(bb)*255))[2:]
        if len(rout)==1:
            rout='0'+rout
        if len(gout)==1:
            gout='0'+gout
        if len(bout)==1:
            bout='0'+bout

        
        color.append('#%2s%2s%2s'%(rout,gout,bout))
        obsdatax.append([])
        obsdatay.append([])
        obsdatayerr.append([])

    if 'setgray' in line:
        ecount+=1
        color.append('#444444')
        obsdatax.append([])
        obsdatay.append([])
        obsdatayerr.append([])

        
    if 'QQ' in line:
        (xx,err,yy,foo)=line.rstrip().split()
        obsdatax[-1].append(float(xx))
        obsdatay[-1].append(float(yy))
        obsdatayerr[-1].append(float(err))

#/QQ {/Flux exch def /sigma exch def /lambda exch def 
#        lambda log Flux log tophys NP 3 0 360 arc fill
#        NP lambda log Flux log sigma 2.5 div sub pmt
#        lambda log Flux log sigma 2.5 div add plt S
#        } def


        
    if '/ft' in line:
        esed.append([])
        line=inf.readline()
        while 'def doit' not in line:
            flux=float(line.strip())
            esed[-1].append(flux)
            line=inf.readline()

datelabel=[]
cshfile=os.popen("ls run*.csh").readline().rstrip()
print("\n***Using "+cshfile+" to get MJDs. Make sure this is right***\n\n")
cshlines=open(cshfile)
readbit=0
for line in cshlines.readlines():
    dat=line.rstrip().split(',')
    if len(dat)!=12:
        continue
    mjd=float(dat[0])
    (yy,mm,dd)=jd2ymd(2400000.5+mjd)
    dateout="%4i-%02i-%02i"%(yy,mm,dd)
    datelabel.append(dateout)


color=["#000000","#ff0000","#0000ff","#dd00dd","#ee7700","#00ee77","#999999","#cccc55","#55cccc","#cc55cc"]
lstyle=["solid","dashed","dotted"]

    
plt.figure(figsize=[6,4])
plt.subplots_adjust(left=0.15,right=0.95,top=0.95,bottom=0.15)
ax=plt.gca()
lowy=1e99
highy=0
for i in range(len(esed)):

    plt.plot(wave,esed[i],color=color[i%10],ls=lstyle[i%3],label=datelabel[i])    

    if min(esed[i])<lowy:
        lowy=min(esed[i])
    if max(esed[i])>highy:
        highy=max(esed[i])

    
    ax.set_xscale('log')
    ax.set_yscale('log')

    erroroutl=[]
    errorouth=[]
    obsdatay_go=[]
    for j in range(len(obsdatay[i])):
        if obsdatayerr[i][j]>0:    
            obsdatay_go.append(obsdatay[i][j])

        else:
            #measured flux is negative, so plot at an artifically low point, and make error bars correct
            obsdatay_go.append(1e-99)
            
        erroroutl.append( obsdatay[i][j] * (1 - 1/(10**(obsdatayerr[i][j]/2.5))))
        errorouth.append( obsdatay[i][j] * (10**(obsdatayerr[i][j]/2.5) -1))  
        
    plt.errorbar(obsdatax[i],obsdatay_go,yerr=[erroroutl,errorouth],ecolor=color[i%10],fmt='o',color=color[i%10])
    
limtest=ax.get_ylim()
if limtest[0]<1e-99:
    plt.ylim(lowy/10.,highy*10)
plt.legend(loc=2)
plt.xlabel("Wavelength (microns)",fontsize=12)
plt.ylabel(r"$\nu$F$_\nu$",fontsize=12)
plt.savefig("bestfit_SED.png")
plt.savefig("bestfit_SED.pdf")
plt.close()


#always rerun, to get analysis parameters
#if os.path.isfile("fort.4"):
#    print("PJDFC parser output files already exist; using those")
#else:

os.system("echo 'PJDFC.out' | ../TPM/read-WISE-rc-MCMC-PJDFC") 

os.system("/bin/cp fort.2 Dhist.dat")
os.system("/bin/cp fort.32 Dhist_fine.dat")
os.system("/bin/cp fort.3 DvsPeriod.dat")
os.system("/bin/cp fort.4 DvsAlb.dat")

#Dhist:
#Ned's code spits out a list of numbers (index 1 below)
#entry 1 = 1 m
#entry 24 = 10 m
#entry 47 = 100 m
#entry 70 = 1 km
#final entry 123 = 200 km
#Dhist_fine has a factor of 10 more bins, each 1/10th the width
#But, I'm just going to use the diameters from DvsAlb to build my own histograms


D=[]
pV=[]
inf=open("DvsAlb.dat")
for line in inf.readlines():
    dat=line.rstrip().split()
    D.append(float(dat[0]))
    pV.append(float(dat[1]))

    
plt.figure(figsize=[6,6])
ax=plt.gca()
#logD=[numpy.log10(d) for d in D]
#logpV=[numpy.log10(p) for p in pV]
#plt.scatter(logD,logpV,s=1.5,c='black',marker='o',edgecolors='')
plt.scatter(D,pV,s=1,marker='o',color='k')
ax.set_xscale('log')
ax.set_yscale('log')

(Hbin,xbin,ybin)=numpy.histogram2d(D,pV,bins=30)
xbino=[(xbin[i+1] - xbin[i])/2. + xbin[i] for i in range(len(xbin)-1)]
ybino=[(ybin[i+1] - ybin[i])/2. + ybin[i] for i in range(len(ybin)-1)]
Hout=numpy.transpose(Hbin)
Hmax=numpy.amax(Hout)

foo=list(Hout.flatten())
foo.sort()
totH=sum(foo)

sumH=0
sig1=-1
sig2=-1
for i in range(len(foo)):
    sumH+=foo[i]
    if sumH>0.05*totH and sig2<0:
        sig2=i
    if sumH>0.32*totH and sig1<0:
        sig1=i

#print(sig2,foo[sig2],sig1,foo[sig1])

#plt.contour(xbino,ybino,Hout,colors='#777777',levels=[0.25*Hmax,0.5*Hmax,0.75*Hmax])
#plt.text(0.2,-0.2,"Contours at 25%,\n50%, 75% of peak")
plt.contour(xbino,ybino,Hout,colors='#cc0000',levels=[foo[sig2],foo[sig1]])
plt.text(max(D)/2.,0.5,"Contours contain 68%\nand 95.5% of all points")
#plt.ylim(-2,0)
plt.ylim(0.01,1)
plt.xlabel("diameter (km)")
plt.ylabel("albedo")
plt.savefig("DvsAlb.png")
plt.savefig("DvsAlb.pdf")
plt.close()

logD=[numpy.log10(d) for d in D]
foo=logD[:]
nD=len(logD)
foo.sort()
sig1low=foo[int(nD*0.16)]
sig1high=foo[int(nD*0.84)]
sig2low=foo[int(nD*0.025)]
sig2high=foo[int(nD*0.975)]

if sig2high-sig2low > 0.2:
    histstep=0.01
elif sig2high-sig2low > 0.02:
    histstep=0.001
else:
    histstep=0.0001

plt.figure(figsize=[6,6])
histlimlow=int(min(logD)*1000)/1000.
histlimhigh=(int(max(logD)*1000)+1)/1000.
plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),histtype='step',color='black')


plt.axvline(sig1low,color='#cc0000',ls='dashed')
plt.axvline(sig1high,color='#cc0000',ls='dashed')
plt.axvline(sig2low,color='#cc0000',ls='dotted')
plt.axvline(sig2high,color='#cc0000',ls='dotted')

plt.xlabel("log diameter (km)")
plt.ylabel("number of Monte Carlo results")
plt.savefig("diam_hist.png")


Dper=[]
period=[]
inf=open("DvsPeriod.dat")
for line in inf.readlines():
    dat=line.rstrip().split()
    Dper.append(float(dat[0]))
    period.append(float(dat[1]))

if max(period)-min(period)==0:
    print("Fixed period used, skipping D vs Period plot")
else:

    plt.figure(figsize=[6,6])
    
    plt.scatter(Dper,period,s=1.5,c='black',marker='o',edgecolors='')
    
    (Hbin,xbin,ybin)=numpy.histogram2d(Dper,period,bins=20)
    xbino=[(xbin[i+1] - xbin[i])/2. + xbin[i] for i in range(len(xbin)-1)]
    ybino=[(ybin[i+1] - ybin[i])/2. + ybin[i] for i in range(len(ybin)-1)]
    Hout=numpy.transpose(Hbin)
    Hmax=numpy.amax(Hout)
    
    foo=list(Hout.flatten())
    foo.sort()
    totH=sum(foo)
    
    sumH=0
    sig1=-1
    sig2=-1
    for i in range(len(foo)):
        sumH+=foo[i]
        if sumH>0.05*totH and sig2<0:
            sig2=i
        if sumH>0.32*totH and sig1<0:
            sig1=i
    
    plt.contour(xbino,ybino,Hout,colors='#cc0000',levels=[foo[sig2],foo[sig1]])
    plt.text(numpy.mean(Dper),0.95*max(period),"Contours contain 68%\nand 95.5% of all points")
    plt.ylim(-2,0)
    plt.xlabel("diameter (km)")
    plt.ylabel("rotation period")
    plt.savefig("DvsPer.png")
    plt.close()
