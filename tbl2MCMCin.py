#this code will take a series of tbl files, and generate the input
#file needed for Ned's MCMC thermophysical model.

#To run, make a separate subdir for each object with all irsa_*.tbl
#output files, then run this in there

#The output of this is a series of STDOUT statements that need to be put into a csh script, then executed
#The output of the csh-called fortran routine will then need to be parsed with a later TBD script

#2018-11-27  JRM

#I will bite the bullet, and write this in python 3

def robust_avg(meas,sig):

    #astropy table reader provides a masked column, so use that to skip bad values
    #meas=[meas_masked[i] for i in range(len(meas_masked)) if not meas_masked.mask[i]]
    #sig=[sig_masked[i] for i in range(len(sig_masked)) if not sig_masked.mask[i]]
    #this has been moved before the call
 
    
    #get the robust average of a list of data with sigmas
    np=len(meas)

    emx=2.0 #for points more than 2 sig away from the mean, use a linear extrapolation for the chisq

    #print(meas,sig)
    
    for i in range(len(sig)):
        if sig[i]==0:
            sig[i]=1e-4    
    
    if np!=len(sig):
        warnings.warn("Error: Number of measurements and sigmas must be equal")
        mu=9.99
        sigma=9.99
        rcs=0
    elif np==0:
        mu=9.99
        sigma=9.99
        rcs=0
    elif np==1:
        mu=meas[0]
        sigma=sig[0]
        rcs=0
    else:
        m_min=min(meas)
        m_max=max(meas)
        xc = (m_min+m_max)/2.
        dmax= m_max-m_min
        if dmax<=0:
            mu=xc
            sigma=(1/sum([1/s**2 for s in sig]))**0.5
            rcs=0
        else:
            dx=[x-xc for x in meas]
            dx.sort()

            #use a len 100 array per element to search for interpolated minimum
            if np<10:
                nx=100
            else:
                nx=30
            n=nx*np
            newdx=numpy.zeros(n)

            for i in range(np-1):
                for jx in range(nx):
                    f=jx/nx  #fraction of the way from 0 to 0.999
                    j=int(jx + i*nx) #location in new array
                    out=f * dx[i+1] + (1-f)*dx[i]
                    #print(f,dx[i],dx[i+1],out)
                    newdx[j]=out

            Y=numpy.zeros(n)
            for i in range(n):
                Y[i]=0
                for j in range(np):
                    E = abs((meas[j] - xc - newdx[i])/sig[j])
                    if E > emx:
                        CS = 2 * emx * E - emx**2
                    else:
                        CS=E**2
                    Y[i]+=CS
            best=min(Y)
            best_index=numpy.where(Y==best)[0][0]

            if best_index==0:
                D = (newdx[2] - newdx[0])/2.
            elif best_index==n-1:
                D = (newdx[n-1] - newdx[n-3])/2.
            else:
                D = (newdx[best_index+1] - newdx[best_index-1])/2.

            D=max(D,dmax/1024.) #relic from Ned's code. I don't know why it's here; or rather, I don't know where 1024 came from
                                #it appears the 1024 limits the iterations below to 10 passes if it doesn't converge

            ddy = 0
            newy=[0,0,0]
            while ddy<=0 and D<dmax:
                D*=2
                for i in range(3):
                    newy[i]=0
                    T = newdx[best_index] + (i-1)*D
                    for j in range(np):
                        E = abs((meas[j]-xc-T)/sig[j])
                        if E > emx:
                            CS=2 * emx * E - emx**2
                        else:
                            CS=E**2
                        newy[i]+=CS
                ddy=newy[2] - 2*newy[1] + newy[0]
                dy = (newy[2] - newy[0])/2.
            mu= xc + newdx[best_index] - D*dy/ddy
            sigma = D * (2/ddy)**0.5
            rcs=0
            for j in range(np):
                E = abs((meas[j]-mu)/sig[j])
                if E > emx:
                    CS=2*emx*E - emx**2
                else:
                    CS=E**2
                rcs+=CS
    
    
    return mu,sigma,rcs 

    

def get_dist(n,mjd):
    #for object name n, get geocentric distance at time mjd

    #example:
    #https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE='OBSERVER'&COMMAND='155140'&CENTER='500'&REF_SYSTEM=J2000&START_TIME=JD2455400.5&STOP_TIME=JD2455400.6&STEP_SIZE=1d&QUANTITIES='20'&CSV_FORMAT='YES'

    jd=mjd+2400000.5
    startt=str(jd)
    stopt=str(jd+0.1)

    #note, single-quotes have to be converted to %27 for horizons to parse this correctly
    call_line="https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=%27OBSERVER%27&COMMAND=%27"+n+"%27&CENTER=%27500@wise%27&REF_SYSTEM=J2000&START_TIME=JD"+startt+"&STOP_TIME=JD"+stopt+"&STEP_SIZE=1d&QUANTITIES=%2720%27&CSV_FORMAT=%27YES%27"
    
    horizons_ret=urlopen(call_line).readlines()
    for i in range(len(horizons_ret)):
        if '$$SOE' in str(horizons_ret[i]):
            break
    dat=horizons_ret[i+1]
    dist=float(str(dat).split(',')[3])
    return dist
    


def unpackMPCname(s):
    nn2c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    if len(s)==5:
        if s[4]=='P':
            #periodic comet
            outn=str(int(s[0:4]))+'P'
        elif s[4]=='S':
            #satellite
            outn=s[0]+'/'+str(int(s[1:4]))+'S'
        else:
            #numbered obj
            outn='('+str(nn2c.index(s[0])*10000 + int(s[1:]))+')'
    elif len(s)==7:
        #PLS object
        if s[0:3] in ['PLS','T1S','T2S','T3S']:
            outn=s[3:]+" "+s[0]+'-'+s[1]
        else:
            #unnumbered asteroid
            year=nn2c.index(s[0])*100 + int(s[1:3])
            if s[4:6]=='00':
                prov=s[3]+s[6]
            else:
                prov=s[3]+s[6]+str(nn2c.index(s[4])*10 + int(s[5]))
            outn=str(year)+" "+prov
    elif len(s)==8:
        #parabolic comet
        year=nn2c.index(s[1])*100 + int(s[2:4])        
        if s[7] not in ['0','1','2','3','4','5','6','7','8','9']:
            prov=s[4]+s[7]+str(nn2c.index(s[5])*10 + int(s[6]))
            outn=s[0]+'/'+str(year)+" "+prov
        else:
            outn=s[0]+'/'+str(year)+" "+s[4]+str(int(s[5:7]))
    else:
        print("Did not recognize input format of:",s)
        outn=s
    return outn


#MAIN 


ZP={'w1':20.752,'w2':19.596,'w3':17.800,'w4':12.945}
nonlin={'w1':4,'w2':6.1,'w3':4,'w4':0}

epoch_break=10 #days

from astropy.table import Table

import numpy
import os
import warnings
import sys

try:
    from urllib.request import urlopen
except ImportError:
    print("You must run this with python3")
    exit()

if len(sys.argv)<2:
    print("Runline: tbl2MCMC.py (packed designation)")
    exit()
    

objname=sys.argv[1]
    
#because horizons is an https, this lets urllib accept unverified ssl certificates
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

#write header

outf=open("run_"+objname+".csh",'w')

print("gfortran -o WISE-steep-MCMC-PJDFC ~/MCMC/TPM/WISE-steep-MCMC-PJDFC.f ~/MCMC/TPM/Nedlib.a")
outf.write("gfortran -o WISE-steep-MCMC-PJDFC ~/MCMC/TPM/WISE-steep-MCMC-PJDFC.f ~/MCMC/TPM/Nedlib.a\n")

#objname=os.popen("pwd").readline().rstrip().split('/')[-1]

if objname=="Test":
    objname="F5140"


print("echo ",objname)
outf.write("echo "+objname+'\n')

mpcorb=os.popen("grep '^"+objname+"' ~/MCMC/mpcorb.s3m").readline()
H=mpcorb.split()[8]
G=mpcorb.split()[9]
Herr=0.2

#get LC period here ^*^

lcdb_name=unpackMPCname(objname).replace('(',"").replace(')',"")
if ' ' in lcdb_name:
    lcdb_back=os.popen("grep ' "+lcdb_name+" ' ~/MCMC/lc_summary.tab").readline()
else:
    lcdb_back=os.popen("grep '^[ ]*"+lcdb_name+" ' ~/MCMC/lc_summary.tab").readline()

try:
    period=float(lcdb_back[85:95])
except:
    period=0

    
print("echo H="+H+" from mpcorb.s3m")
outf.write("echo H="+str(H)+" from mpcorb.s3m\n")


print("./WISE-steep-MCMC-PJDFC << LAST")
outf.write("./WISE-steep-MCMC-PJDFC << LAST\n")

print(H+","+str(Herr)+","+str(period)+","+objname)
outf.write(str(H)+","+str(Herr)+","+str(period)+","+objname+'\n')

#print(H+","+G+","+objname)

infs=os.popen("ls irsa*.tbl").readlines()


for inline in infs:
    epochs=[]
    #inf=atpy.Table(inline.rstrip())
    inf=Table.read(inline.rstrip(),format='ipac')

    allmjd=list(inf["mjd"])
    allmjd.sort()

    mjd_start=allmjd[0]
    for i in range(len(allmjd)-1):
        if (allmjd[i+1]- allmjd[i])>epoch_break:
            epochs.append((mjd_start,allmjd[i]))
            mjd_start=allmjd[i+1]
    epochs.append((mjd_start,allmjd[-1]))


    for e in epochs:
        #e_data=inf.where((inf.mjd>=e[0])&(inf.mjd<=e[1]))
        e_data=inf[ (inf['mjd']>=e[0]) & (inf['mjd']<=e[1])] 

        #we want to use an actual measured mjd, not a mean here, so that we can request the wise-centric distance accurately
        #so, use the int divide, not float
        
        midi=int(len(e_data)/2)
        
        midtime=e_data["mjd"][midi]
        midra=e_data["ra"][midi]
        middec=e_data["dec"][midi]

        horizon_name=unpackMPCname(objname).replace('(',"").replace(')',"").replace(" ",'%27')
        if len(objname)==5:
            horizon_name=horizon_name+'%3B'
        
        #get geo dist from horizons here for midpoint
        obs_dist=get_dist(horizon_name,midtime)

        oline="%9.3f,%7.3f,%+7.3f,%7.5f"%(midtime,midra,middec,obs_dist)

        for band in ['w1','w2','w3','w4']:
            if band+'flux' in e_data.columns and "_3b" not in inline:

                meas_masked=e_data[band+"flux"]
                sig_masked=e_data[band+"sigflux"]

                try:
                    meas=[meas_masked[i] for i in range(len(meas_masked)) if not (meas_masked.mask[i] or sig_masked.mask[i])]
                    sigs=[sig_masked[i] for i in range(len(sig_masked)) if not (meas_masked.mask[i] or sig_masked.mask[i])]
                except AttributeError:
                    #if there are no bad values, a mask object is not created, which is soooooo helpful
                    meas=meas_masked
                    sigs=sig_masked

                (meanflux,meanferr,rchi)=robust_avg(meas,sigs)

                meanmag= ZP[band] - 2.5*numpy.log10(abs(meanflux))
                
                rsig = 2.5*numpy.log10(1 + meanferr/abs(meanflux))


                if len(meas)>1:
                    meansig=rsig*(max(1,rchi/(len(meas)-1)))**0.5
                else:
                    meansig=rsig

                meansig=max(meansig,0.03) #set minimum mag error based on repeatability

                if meanflux<0:
                    meansig*=-1
                elif meanmag<nonlin[band]:
                    if band=='w2':
                        #linear fit to the K-W2 plot in the WISE explanatory supplement, sec6_3c, figure 8b.
                        magoff=(0-0.55)/(6.1-3.5)*meanmag + 1.29
                        meanmag=meanmag+magoff
                    if band=='w3':
                        #from Wright et al. 2018 arXiv posting
                        updated_mag = 0.86*meanmag + 0.49
                        meanmag=updated_mag

                    #dropping this because for MBAs it was forcing the fits to be be dominated by 2-band data
                    #meansig=0.2
                
                oline=oline+",%5.2f,%+5.2f"%(meanmag,meansig)


            elif band+'mpro' in e_data.columns:

                if "_3b" not in inline:
                    #if fluxes are not present for some reason, use mags instead
                    warnings.warn("Can't find flux for "+band+", using mpro instead")

                meas_masked=e_data[band+"mpro"]
                sig_masked=e_data[band+"sigmpro"]

                try:
                    meas=[meas_masked[i] for i in range(len(meas_masked)) if not (meas_masked.mask[i] or sig_masked.mask[i])]
                    sigs=[sig_masked[i] for i in range(len(sig_masked)) if not (meas_masked.mask[i] or sig_masked.mask[i])]
                except AttributeError:
                    #if there are no bad values, a mask object is not created, which is soooooo helpful
                    meas=meas_masked
                    sigs=sig_masked

                if len(sigs)>0:

                        
                    (meanmag,meansig,rchi)=robust_avg(meas,sigs)
                    meansig=max(meansig,0.03) #set minimum mag error based on repeatability
    
                    if meanmag<nonlin[band]:
                        if band=='w2':
                            #linear fit to the K-W2 plot in the WISE explanatory supplement, sec6_3c, figure 8b.
                            magoff=(0-0.55)/(6.1-3.5)*meanmag + 1.29
                            meanmag=meanmag+magoff
                        if band=='w3':
                            #from Wright et al. 2018 arXiv posting
                            updated_mag = 0.86*meanmag + 0.49
                            meanmag=updated_mag
    
                    #    meansig=0.2
    
                    oline=oline+",%5.2f,%+5.2f"%(meanmag,meansig)
                            
                else:
                    #band totally dropped
                    oline=oline+", 9.99,+9.99"
            else:
                #band missing
                oline=oline+", 9.99,+9.99"
        print(oline)
        outf.write(oline+'\n')

print("LAST")
outf.write("LAST\n")

outf.close()
os.system("chmod +x run_"+objname+".csh")

print("\n\nnow run with:")
print("./run_"+objname+".csh")

        

