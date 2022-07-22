#this code will take a series of tbl files, and generate the input
#file needed for Ned's MCMC thermophysical model.

#To run, make a separate subdir for each object with all irsa_*.tbl
#output files, then run this in there

#The output of this is a series of STDOUT statements that need to be put into a csh script, then executed
#The output of the csh-called fortran routine will then need to be parsed with a later TBD script

#2018-11-27  JRM

#I will bite the bullet, and write this in python 3

import numpy as np
import os
import warnings
import sys
import ssl
from astropy.table import Table
try:
    from urllib.request import urlopen
except ImportError:
    print("You must run this with python3")
    exit()


def robust_avg(measurements, sigmas):

    #astropy table reader provides a masked column, so use that to skip bad values
    #meas=[meas_masked[i] for i in range(len(meas_masked)) if not meas_masked.mask[i]]
    #sig=[sig_masked[i] for i in range(len(sig_masked)) if not sig_masked.mask[i]]
    #this has been moved before the call
 
    
    #get the robust average of a list of data with sigmas
    len_measurements = len(measurements)
    emx = 2.0 #for points more than 2 sig away from the mean, use a linear extrapolation for the chisq
    
    for i in range(len(sigmas)):
        if sigmas[i] == 0:
            sigmas[i] = 1e-4    
    
    if len_measurements != len(sigmas):
        warnings.warn("Error: Number of measurements and sigmas must be equal")
        mu = 9.99
        sigma = 9.99
        rcs = 0
    elif len_measurements == 0:
        mu = 9.99
        sigma = 9.99
        rcs = 0
    elif len_measurements == 1:
        mu = measurements[0]
        sigma = sigmas[0]
        rcs = 0
    else:
        m_min = min(measurements)
        m_max = max(measurements)
        xc = (m_min + m_max) / 2.
        d_max = m_max - m_min
        if d_max <= 0:
            mu = xc
            sigma = (1 / sum([1 / s ** 2 for s in sigmas])) ** 0.5
            rcs = 0
        else:
            dx = [x - xc for x in measurements]
            dx.sort()

            #use a len 100 array per element to search for interpolated minimum
            if len_measurements < 10:
                nx = 100
            else:
                nx = 30
            n = nx * len_measurements
            new_dx = np.zeros(n)

            for i in range(len_measurements - 1):
                for jx in range(nx):
                    f = jx / nx  #fraction of the way from 0 to 0.999
                    j = int(jx + i * nx) #location in new array
                    out = f * dx[i + 1] + (1 - f) * dx[i]
                    #print(f,dx[i],dx[i+1],out)
                    new_dx[j] = out

            y = np.zeros(n)
            
            for i in range(n):
                y[i] = 0
                
                for j in range(len_measurements):
                    e = abs((measurements[j] - xc - new_dx[i]) / sigmas[j])
                    
                    if e > emx:
                        cs = 2 * emx * e - emx ** 2
                    
                    else:
                        cs = e ** 2
                    y[i] += cs
            
            best = min(y)
            best_index = np.where(y == best)[0][0]

            if best_index == 0:
                d = (new_dx[2] - new_dx[0]) / 2.
            elif best_index == n - 1:
                d = (new_dx[n - 1] - new_dx[n - 3]) / 2.
            else:
                d = (new_dx[best_index + 1] - new_dx[best_index - 1]) / 2.

            d = max(d, d_max / 1024.) #relic from Ned's code. I don't know why it's here; or rather, I don't know where 1024 came from
                                #it appears the 1024 limits the iterations below to 10 passes if it doesn't converge
            ddy = 0
            new_y = [0, 0, 0]
            
            while ddy <= 0 and d < d_max:
                d *= 2
                
                for i in range(3):
                    new_y[i] = 0
                    t = new_dx[best_index] + (i - 1) * d
                    
                    for j in range(len_measurements):
                        e = abs((measurements[j] - xc - t) / sigmas[j])
                        
                        if e > emx:
                            cs = 2 * emx * e - emx ** 2
                        
                        else:
                            cs = e ** 2
                        new_y[i] += cs
                
                ddy = new_y[2] - 2 * new_y[1] + new_y[0]
                dy = (new_y[2] - new_y[0]) / 2.
            mu = xc + new_dx[best_index] - d * dy / ddy
            sigma = d * (2 / ddy) ** 0.5
            rcs = 0
            
            for j in range(len_measurements):
                e = abs((measurements[j] - mu) / sigmas[j])
                
                if e > emx:
                    cs = 2 * emx * e - emx ** 2
                
                else:
                    cs = e ** 2
                rcs += cs
    
    return mu,sigma,rcs 


def get_distance(name, mjd):
    """
    Retrieves the geocentric distance of an object at a given time.

    Arguments: name (str) -- The MPC unpacked designation of the object
               mjd (float) -- the time in modified julian days
    """
    #example:
    #https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE='OBSERVER'&COMMAND='155140'&CENTER='500'&REF_SYSTEM=J2000&START_TIME=JD2455400.5&STOP_TIME=JD2455400.6&STEP_SIZE=1d&QUANTITIES='20'&CSV_FORMAT='YES'

    jd = mjd + 2400000.5
    start_time = str(jd)
    stop_time = str(jd + 0.1)

    # Note: single-quotes have to be converted to %27 for horizons to 
    # parse this correctly
    call_line = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&T" 
    call_line += "ABLE_TYPE=%27OBSERVER%27&COMMAND=%27" + name
    call_line += "%27&CENTER=%27500@wise%27&REF_SYSTEM=J2000&START_TIM" 
    call_line += "E=JD" + start_time + "&STOP_TIME=JD" + stop_time
    call_line += "&STEP_SIZE=1d&QUANTITIES=%2720%27&CSV_FORMAT=%27YES%27"

    horizons_file = urlopen(call_line).readlines()
    for i in range(len(horizons_file)):
        if "$$SOE" in str(horizons_file[i]):
            break
    datum = horizons_file[i + 1]
    distance = float(str(datum).split(',')[3])
    return distance
    


def unpack_MPC_name(packed_name):
    """
    Generates the unpacked MPC designation for a given object.

    Arguments: packed_name (str) -- The packed name of the object.

    Returns: (str) -- The unpacked name of the object
    """

    char_map = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    char_map += "abcdefghijklmnopqrstuvwxyz"

    if len(packed_name) == 5:

        if packed_name[4] == 'P':
            # Periodic Comet
            unpacked_name = str(int(packed_name[0:4])) + 'P'

        elif packed_name[4] == 'S':
            # Satellite
            unpacked_name = packed_name[0] + '/'
            unpacked_name += str(int(packed_name[1:4])) + 'S'

        else:
            # Numbered Object
            input_1 = char_map.index(packed_name[0]) * 10000
            input_2 = int(packed_name[1:])
            unpacked_name = '(' + str(input_1 + input_2) + ')'

    elif len(packed_name) == 7:

        # PLS Object
        if packed_name[0:3] in ["PLS", "T1S", "T2S", "T3S"]:
            unpacked_name = packed_name[3:] + " " + packed_name[0]
            unpacked_name += '-' + packed_name[1]

        else:
            # Unnumbered Asteroid
            year = char_map.index(packed_name[0]) * 100
            year += int(packed_name[1:3])

            if packed_name[4:6] == "00":
                prov = packed_name[3] + packed_name[6]

            else:
                input_1 = char_map.index(packed_name[4]) * 10
                input_2 = int(packed_name[5])
                prov = packed_name[3] + packed_name[6]
                prov += str(input_1 + input_2)
    
            unpacked_name = str(year) + " " + prov
    
    elif len(packed_name) == 8:
        # Parabolic Comet
        year = char_map.index(packed_name[1]) * 100
        year += int(packed_name[2:4])
        
        if packed_name[7] not in ['0', '1', '2', '3', '4', '5', 
                                  '6', '7', '8', '9']:
            prov = packed_name[4] + packed_name[7]
            input_1 = char_map.index(packed_name[5]) * 10
            input_2 = int(packed_name[6])
            prov += str(input_1 + input_2)
            unpacked_name = packed_name[0] + '/'
            unpacked_name += str(year) + " " + prov
        
        else:
            unpacked_name = packed_name[0] + '/' + str(year)
            unpacked_name += " " + packed_name[4]
            unpacked_name += str(int(packed_name[5:7]))
    
    else:
        print("Did not recognize input format of:", packed_name)
        unpacked_name = packed_name
    
    return unpacked_name


#MAIN 


z_p = {"w1" : 20.752, "w2" : 19.596, "w3" : 17.800, "w4" : 12.945}
non_lin = {"w1" : 4, "w2" : 6.1, "w3" : 4, "w4" : 0}

epoch_break = 10 #days

if len(sys.argv)<2:
    print("Runline: tbl2MCMC.py (packed designation)")
    exit()
    
obj_name = sys.argv[1]
    
# Because horizons is an https, this lets urllib accept unverified 
# ssl certificates
ssl._create_default_https_context = ssl._create_unverified_context

#write header
output_file = open("run_" + obj_name + ".csh", 'w')

print("gfortran -o ../WISE-steep-MCMC-PJDFCS ../WISE-steep-MCMC-PJDFCS.f ../READ_TABLE-steep.f ../Nedlib.a")
#output_file.write("gfortran -o ../WISE-steep-MCMC-PJDFCS ../WISE-steep-MCMC-PJDFCS.f ../READ_TABLE-steep.f ../Nedlib.a\n")

#objname=os.popen("pwd").readline().rstrip().split('/')[-1]

if obj_name == "Test":
    obj_name = "F5140"


print("echo ", obj_name)
output_file.write("echo " + obj_name + '\n')

mpcorb = os.popen("grep '^" + obj_name + "' ./mpcorb.s3m").readline()
print(mpcorb)
h_mag = mpcorb.split()[8]
g_mag = mpcorb.split()[9]
h_error = 0.2

#get LC period here ^*^

lcdb_name = unpack_MPC_name(obj_name).replace('(',"").replace(')',"")
if ' ' in lcdb_name:
    lcdb_back = os.popen("grep ' " + lcdb_name + " ' lc_summary.tab").readline()
else:
    lcdb_back = os.popen("grep '^[ ]*" + lcdb_name + " ' lc_summary.tab").readline()

try:
    period = float(lcdb_back[85:95])
except:
    period = 0

    
print("echo H=" + h_mag + " from mpcorb.s3m")
output_file.write("echo H=" + str(h_mag) + " from mpcorb.s3m\n")
print("../WISE-steep-MCMC-PJDFCS << LAST")
output_file.write("../WISE-steep-MCMC-PJDFCS << LAST\n")
print(h_mag + "," + str(h_error) + "," + str(period) + "," + obj_name)
output_file.write(str(h_mag) + "," + str(h_error) + "," + str(period) + "," + obj_name + '\n\n')
files = os.popen("ls irsa*.tbl").readlines()
print(files)

# Opening all input data files
for inline in files:
    epochs = []
    data_table = Table.read(inline.rstrip(), format="ipac")

    # Reading all mjd
    all_mjd = list(data_table["mjd"])
    all_mjd.sort()

    mjd_start = all_mjd[0]
    
    for i in range(len(all_mjd) - 1):
        if (all_mjd[i + 1] - all_mjd[i]) > epoch_break:
            epochs.append((mjd_start, all_mjd[i]))
            mjd_start = all_mjd[i + 1]
    
    epochs.append((mjd_start, all_mjd[-1]))
    print("Writing to .csh file")

    for epoch in epochs:
        
        # Generates data table containing only epoch relevant data
        epoch_data = data_table[(data_table['mjd'] >= epoch[0]) & (data_table['mjd'] <= epoch[1])] 
        input_table = Table()
        #print(type(input_table))

        # Finds middle index value for accurate middle distance calling
        mid_index = int(len(epoch_data) / 2)
        mid_time = epoch_data["mjd"][mid_index]
        mid_ra = epoch_data["ra"][mid_index]
        mid_dec = epoch_data["dec"][mid_index]

        horizon_name = unpack_MPC_name(obj_name).replace('(', "").replace(')', "").replace(" ", "%27")
        if len(obj_name) == 5:
            horizon_name = horizon_name + "%3B"
        
        mid_distance = get_distance(horizon_name, mid_time)
        output_line = "%9.3f,%7.3f,%+7.3f,%7.5f"%(mid_time, mid_ra, mid_dec, mid_distance)
        for band in ["w1", "w2", "w3", "w4"]:

            # For non three band data and present flux
            if band + "flux" in epoch_data.columns and "_3b" not in inline:

                masked_measurements = epoch_data[band + "flux"]
                masked_sigmas = epoch_data[band + "sigflux"]

                try:
                    measurements = [masked_measurements[i] for i in range(len(masked_measurements)) if not (masked_measurements.mask[i] or masked_sigmas.mask[i])]
                    sigmas = [masked_sigmas[i] for i in range(len(masked_sigmas)) if not (masked_measurements.mask[i] or masked_sigmas.mask[i])]
                
                except AttributeError:
                    # All measurements are good data points
                    measurements = masked_measurements
                    sigmas = masked_sigmas

                mean_flux, mean_ferr, r_chi = robust_avg(measurements, sigmas)
                mean_mag = z_p[band] - 2.5 * np.log10(abs(mean_flux))
                r_sig = 2.5 * np.log10(1 + mean_ferr / abs(mean_flux))

                if len(measurements) > 1:
                    mean_sig = r_sig * (max(1, r_chi / (len(measurements) - 1))) ** 0.5
                else:
                    mean_sig = r_sig

                # set minimum mag error based on repeatability
                mean_sig = max(mean_sig, 0.03) 

                if mean_flux < 0:
                    mean_sig *= -1

                elif mean_mag < non_lin[band]:
                    
                    if band == "w2":
                        #linear fit to the K-W2 plot in the WISE explanatory supplement, sec6_3c, figure 8b.
                        mag_off = (0 - 0.55) / (6.1 - 3.5) * mean_mag + 1.29
                        mean_mag = mean_mag + mag_off
                    
                    if band == "w3":
                        #from Wright et al. 2018 arXiv posting
                        updated_mag = 0.86 * mean_mag + 0.49
                        mean_mag = updated_mag
                
                mean_amp = (0.5 * abs(np.max(epoch_data[band + "flux"]) - np.min(epoch_data[band + "flux"])) / np.mean(epoch_data[band + "flux"])) * 5 * np.log10(np.e)
                output_line = output_line + ",%5.2f,%+5.2f,%5.2f"%(mean_mag, mean_sig, mean_amp)

            # For lacking flux
            elif band + "mpro" in epoch_data.columns:
                
                # For the non three band case
                if "_3b" not in inline:
                    warnings.warn("Can't find flux for " + band + ", using mpro instead")

                masked_measurements = epoch_data[band + "mpro"]
                masked_sigmas = epoch_data[band + "sigmpro"]

                try:
                    measurements = [masked_measurements[i] for i in range(len(masked_measurements)) if not (masked_measurements.mask[i] or masked_sigmas.mask[i])]
                    sigmas = [masked_sigmas[i] for i in range(len(masked_sigmas)) if not (masked_measurements.mask[i] or masked_sigmas.mask[i])]
                
                except AttributeError:
                    # All measurements are good data points
                    measurements = masked_measurements
                    sigmas = masked_sigmas

                if len(sigmas) > 0:
                    mean_mag, mean_sig, r_chi = robust_avg(measurements, sigmas)
                    mean_sig = max(mean_sig, 0.03) # set minimum mag error based on repeatability
    
                    if mean_mag < non_lin[band]:

                        if band == "w2":
                            #linear fit to the K-W2 plot in the WISE explanatory supplement, sec6_3c, figure 8b.
                            mag_off = (0 - 0.55)/(6.1-3.5) * mean_mag + 1.29
                            mean_mag = mean_mag + mag_off
                        
                        if band == 'w3':
                            #from Wright et al. 2018 arXiv posting
                            updated_mag = 0.86 * mean_mag + 0.49
                            mean_mag = updated_mag

                    mean_amp = (0.5 * abs(np.max(epoch_data[band + "flux"]) - np.min(epoch_data[band + "flux"])) / np.mean(epoch_data[band + "flux"])) * 5 * np.log10(np.e)
                    output_line = output_line + ",%5.2f,%+5.2f,%5.2f"%(mean_mag, mean_sig, mean_amp)

                # Band is completely dropped      
                else:
                    output_line = output_line + ", 9.99,+9.99,0."
            
            # Band is not present
            else:

                output_line = output_line + ", 9.99,+9.99,0."
        
        print(output_line)
        output_file.write(output_line + '\n')


output_file.write("LAST\n")
print("LAST")
output_file.close()
os.system("chmod +x run_" + obj_name + ".csh")
print("\n\nnow run with:")
print("./run_" + obj_name + ".csh")

        

