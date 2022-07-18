"""
python3 parser of the output from WISE-steep-MCMC-PJDFC
To be run this in the directory made for the object.
"""
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from math import modf, floor


def julian_days_utc_converter(jd):
    """
    Returns a date in utc format corresponding to the given 
    Julian day number.

    Arguments: jd (float) -- a julian date entry

    Returns: (tup) -- a utc format entry with year[0], month[1], and
                      day[2].
    """
    jd_adjusted = float(jd) + 0.5
    decimal_day = modf(jd_adjusted)[0]
    x_2 = floor(jd - 1721119.5)
    c_2 = floor((4 * x_2 + 3) / 146097)
    x_1 = x_2 - floor(146097 * c_2 / 4)
    c_1 = floor((100 * x_1 + 99) / 36525)
    x_0 = x_1 - floor(36525 * c_1 / 100)
    year = 100 * c_2 + c_1
    month = floor((5 * x_0 + 461) / 153)
    day = x_0 - floor((153 * month - 457) / 5) + 1
    if month > 12:
        month = month - 12
        year = year + 1    
    return year, month, day + decimal_day



def retrieve_MCMC_data():
    #cp, not mv, so that there is a record of the actual outputs  
    os.system("/bin/cp ./fort.21 PJDFC.out")
    #fort.22 is the postscript format of the SED
    os.system("/bin/cp ./fort.22 SED_data.out")

    #^*^ parse and plot SED_data.out here 
    # best fit solution used for SED
    SED_file = open("SED_data.out")
    print("\nReading in best fit solution used for SED")
    datum_minor = SED_file.readline().rstrip().split()
    datum_major = SED_file.readline().rstrip().split()
    best_fit = {}

    best_fit["pole_ra"] = np.degrees(float(datum_minor[1]))
    best_fit["pole_dec"] = np.degrees(float(datum_minor[2]))
    best_fit["pv_median"] = np.e ** (float(datum_minor[3]))
    best_fit["period"] = np.e ** (float(datum_minor[4]))
    best_fit["gamma"] = np.e ** (float(datum_minor[5]))
    best_fit["pjdfc"] = (float(datum_minor[6])) #period * thermal inertia Gamma * diameter * crater frac * color
    x = float(datum_minor[7])
    best_fit["crater_frac"] = np.e ** (x) / (1 + np.e ** x)
    best_fit["ir_albedo_ratio"] = np.e ** (float(datum_minor[8]))
    best_fit["chisq"] = float(datum_minor[10])
    best_fit["penalties"] = float(datum_minor[9]) - best_fit['chisq']
    best_fit["pv"] = float(datum_major[1])
    best_fit["diameter"] = float(datum_major[2])
    best_fit["theta"] = float(datum_major[3])

    # Reads in all epochs and their corresponding conditions
    print("Reading in all epochs")
    SED_file.readline()
    epoch_condition = {}
    epoch_condition["delta"] = []
    epoch_condition["r_helio"] = []
    epoch_condition["phase"] = []
    epoch_condition["sub_sun_lat"] = []
    epoch_condition["sub_earth_lat"] = []

    epoch = SED_file.readline()
    while epoch[0:4] != "/wvl":
        epoch_data = epoch.split()
        epoch_condition["delta"].append(float(epoch_data[1]))
        epoch_condition["r_helio"].append(float(epoch_data[2]))
        epoch_condition["phase"].append(float(epoch_data[3]))
        epoch_condition["sub_sun_lat"].append(float(epoch_data[4]))
        epoch_condition["sub_earth_lat"].append(float(epoch_data[5]))
        epoch = SED_file.readline()

    # Reads in all wavelength data for corresponding epochs
    print("Reading all wavelengths")
    wavelengths = []
    wave_line = SED_file.readline()
    while "def" not in wave_line:
        wavelength_data = wave_line.strip().split()
        for wavelength in wavelength_data:
            wavelengths.append(float(wavelength))
        wave_line = SED_file.readline()
    
    color = []
    esed = []
    data_x = []
    data_y = []
    data_y_error = []

    while True:
        line = SED_file.readline()
        
        if line == "":
            break
        
        if "/ft" in line:
            esed.append([])
            line = SED_file.readline()
            while "def doit" not in line:
                flux = float(line.strip())
                esed[-1].append(flux)
                line = SED_file.readline()

        elif "QQ" in line:
            x_datum, error_datum, y_datum = line.rstrip().split()[0:3]
            data_x[-1].append(float(x_datum))
            data_y[-1].append(float(y_datum))
            data_y_error[-1].append(float(error_datum))

        elif "SRGB" in line:
            rr, gg, bb = line.rstrip().split()[0:3]
            red = hex(int(float(rr) * 255))[2:]
            green = hex(int(float(gg) * 255))[2:]
            blue = hex(int(float(bb) * 255))[2:]
            if len(red) == 1:
                red = '0' + red
            if len(green) == 1:
                green = '0' + green
            if len(blue) == 1:
                blue ='0' + blue
            
            color.append('#%2s%2s%2s'%(red, green, blue))
            data_x.append([])
            data_y.append([])
            data_y_error.append([])

        elif "setgray" in line:
            color.append("#444444")
            data_x.append([])
            data_y.append([])
            data_y_error.append([])

    datelabels = []
    cshfile = os.popen("ls run*.csh").readline().rstrip()
    print("\n***Using "+cshfile+" to get MJDs. Make sure this is right***\n")
    cshlines = open(cshfile)
    for line in cshlines.readlines():
        epoch_data = line.rstrip().split(',')
        if len(epoch_data) != 12:
            continue
        mjd = float(epoch_data[0])
        year, month, day = julian_days_utc_converter(2400000.5 + mjd)
        date = "%4i-%02i-%02i"%(year, month, day)
        datelabels.append(date)

    best_fit_plotters = [esed, wavelengths, datelabels, data_x, data_y, data_y_error]
    return best_fit, best_fit_plotters, epoch_condition, wavelengths


def generate_SED_plot(esed, wavelengths, datelabels, data_x, data_y, data_y_error):
    """
    Generates a plot of best fit for wavelengths across different epochs.
    The outputted graph is in log 10 scale.

    Arguments: esed (list) -- a series of flux lists for each wavelengths
               wavelengths (list) -- a series of wavelengths to plot against
               datelabels (list) -- a set of labels for each epoch
               data_x (list) -- a set of error points for the x axis
               data_y (list) -- a set of error points for the y axis
               data_y_error (list) -- error values for the y axis

    Returns: None, generates two plots titled best_fit_SED.pdf/png
    """

    print("Generating best fit plot for SED")
    
    colors = ["#000000", "#ff0000", "#0000ff", "#dd00dd", "#ee7700", 
        "#00ee77", "#999999", "#cccc55", "#55cccc", "#cc55cc"]
    linestyles = ["solid","dashed","dotted"]

    plt.figure(figsize=[6,4])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    ax = plt.gca()
    y_low = 1e99
    y_high = 0

    for i in range(len(esed)):

        # Self corrects x and y limits
        if min(esed[i]) < y_low:
            y_low = min(esed[i])
        if max(esed[i]) > y_high:
            y_high = max(esed[i])

        plt.plot(wavelengths, esed[i], 
                color=colors[i % 10], 
                ls=linestyles[i % 3], 
                label=datelabels[i])    

        ax.set_xscale("log")
        ax.set_yscale("log")

        lower_errors = []
        higher_errors = []
        data_y_adjusted=[]

        for j in range(len(data_y[i])):
            
            if data_y_error[i][j] > 0:    
                data_y_adjusted.append(data_y[i][j])

            else:
                #measured flux is negative, so plot at an artifically low point, and make error bars correct
                data_y_adjusted.append(1e-99)

            lower_errors.append(data_y[i][j] * (1 - 1 / (10 ** (data_y_error[i][j] / 2.5))))
            higher_errors.append(data_y[i][j] * (10 ** (data_y_error[i][j] / 2.5) - 1))
            
        plt.errorbar(data_x[i], data_y_adjusted, 
                    yerr=[lower_errors, higher_errors],
                    ecolor=colors[i%10],
                    fmt='o',
                    color=colors[i%10])

    plt.legend(loc=2)
    plt.xlabel("Wavelength (microns)", fontsize=12)
    plt.ylabel(r"$\nu$F$_\nu$", fontsize=12)
    plt.savefig("bestfit_SED.png")
    plt.savefig("bestfit_SED.pdf")
    plt.close()


def generate_diameter_vs_albedo_plot():
    """
    Generates a plot of diameter and albedo solutions for the MCMC output.

    Arguments: None

    Returns: None, generates two plots titled diameter_vs_albedo.png/pdf
    """

    diameters = []
    albedos = []
    
    data_file = open("DvsAlb.dat")
    for line in data_file.readlines():
        datum = line.rstrip().split()
        diameters.append(float(datum[0]))
        albedos.append(float(datum[1]))
    
    # Generate base figure
    plt.figure(figsize=[6,6])
    ax = plt.gca()
    plt.scatter(diameters, albedos, s=1,marker='o',color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set up bins for contour levels
    H_bins, x_bins, y_bins = np.histogram2d(diameters, albedos, bins=30)
    x_bins_adjusted = [(x_bins[i+1] - x_bins[i])/2. + x_bins[i] for i in range(len(x_bins)-1)]
    y_bins_adjusted = [(y_bins[i+1] - y_bins[i])/2. + y_bins[i] for i in range(len(y_bins)-1)]
    H_bins_adjusted = np.transpose(H_bins)

    H_bins_flat = list(H_bins_adjusted.flatten())
    H_bins_flat.sort()
    H_total = sum(H_bins_flat)

    H_cumulative = 0
    level_1 =- 1
    level_2 =- 1
    for i in range(len(H_bins_flat)):
        H_cumulative += H_bins_flat[i]
        if H_cumulative > 0.05 * H_total and level_2 < 0:
            level_2 = i
        if H_cumulative > 0.32 * H_total and level_1 < 0:
            level_1 = i

    # Add contour linings to plot
    plt.contour(x_bins_adjusted, y_bins_adjusted, H_bins_adjusted,
                colors="#cc0000",
                levels=[H_bins_flat[level_2], H_bins_flat[level_1]])
    plt.text(max(diameters)/2.,0.5,"Contours contain 68%\nand 95.5% of all points")
    plt.ylim(0.01, 1)
    plt.xlabel("Diameter (km)")
    plt.ylabel("Albedo")
    plt.savefig("diameter_vs_albedo.png")
    plt.savefig("diameter_vs_albedo.pdf")
    plt.close()


def generate_diameter_histogram():
    """
    Generates a histogram of diameter solutions for the MCMC output.

    Arguments: None

    Returns: None, generates two plots titled diameter_histogram.png/pdf
    """

    diameters = []
    data_file = open("DvsAlb.dat")
    for line in data_file.readlines():
        datum = line.rstrip().split()
        diameters.append(float(datum[0]))

    # Set up log scale plot arguments and standard deviation lines
    log_diameters = [np.log10(diameter) for diameter in diameters]
    log_diameters_copy = log_diameters[:]
    diameters_count = len(log_diameters)
    log_diameters_copy.sort()
    sigma_1_low = log_diameters_copy[int(diameters_count * 0.16)]
    sigma_1_high = log_diameters_copy[int(diameters_count * 0.84)]
    sigma_2_low = log_diameters_copy[int(diameters_count * 0.025)]
    sigma_2_high = log_diameters_copy[int(diameters_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    plt.figure(figsize=[6, 6])
    hist_low_limit = int(min(log_diameters) * 1000) / 1000.
    hist_high_limit = (int(max(log_diameters) * 1000) + 1) / 1000.
    plt.hist(log_diameters, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    plt.axvline(sigma_1_low,color='#cc0000',ls='dashed')
    plt.axvline(sigma_1_high,color='#cc0000',ls='dashed')
    plt.axvline(sigma_2_low,color='#cc0000',ls='dotted')
    plt.axvline(sigma_2_high,color='#cc0000',ls='dotted')
    
    # Labeling plot
    plt.xlabel("Log diameter (km)")
    plt.ylabel("Number of Monte Carlo Results")
    plt.savefig("diameter_histogram.png")
    plt.savefig("diameter_histogram.pdf")


def generate_diameter_vs_period_plot():
    """
    Generates a plot of diameter and period solutions for the MCMC output.

    Arguments: None

    Returns: None, generates two plots titled diameter_vs_period.png/pdf
    """

    diameters = []
    periods = []
    data_file = open("DvsPeriod.dat")
    for line in data_file.readlines():
        datum = line.rstrip().split()
        diameters.append(float(datum[0]))
        periods.append(float(datum[1]))

    if max(periods) - min(periods) == 0:
        print("Fixed period used, skipping D vs Period plot")

    else:
        # Set up base figure
        plt.figure(figsize=[6, 6])    
        plt.scatter(diameters, periods, s=1.5, c="black", marker='o', edgecolors='')
        
        # Set up bins for contour levels
        H_bins, x_bins, y_bins = np.histogram2d(diameters, periods, bins=20)
        x_bins_adjusted = [(x_bins[i + 1] - x_bins[i]) / 2. + x_bins[i] for i in range(len(x_bins) - 1)]
        y_bins_adjusted = [(y_bins[i + 1] - y_bins[i]) / 2. + y_bins[i] for i in range(len(y_bins) - 1)]
        H_bins_adjusted  = np.transpose(H_bins)
        H_output = np.amax(H_bins_adjusted)
        
        H_bins_adjusted = list(H_bins_adjusted.flatten())
        H_bins_adjusted.sort()
        H_total = sum(H_bins_adjusted)

        H_cumulative = 0
        level_1 =- 1
        level_2 =- 1
        for i in range(len(H_bins_adjusted)):
            H_cumulative += H_bins_adjusted[i]
            if H_cumulative > 0.05 * H_total and level_2 < 0:
                level_2 = i
            if H_cumulative > 0.32 * H_total and level_1 < 0:
                level_1 = i
        
        # Add contour lines to plot
        plt.contour(x_bins_adjusted, y_bins_adjusted, H_bins_adjusted, 
                    colors="#cc0000", 
                    levels=[H_bins_adjusted[level_2], H_bins_adjusted[level_1]])
        plt.text(np.mean(diameters), 0.95 * max(periods), 
                "Contours contain 68%\nand 95.5% of all points")
        plt.ylim(-2, 0)
        plt.xlabel("Diameter (km)")
        plt.ylabel("Rotation period")
        plt.savefig("diameter_vs_period.png")
        plt.savefig("diameter_vs_period.pdf")
        plt.close()


if sys.version[0] != '3':
    print("Please run with python3")
    exit()

best_fit, best_fit_plotters, epoch_condition, wavelengths = retrieve_MCMC_data()
generate_SED_plot(best_fit_plotters[0], best_fit_plotters[1], 
                  best_fit_plotters[2], best_fit_plotters[3], 
                  best_fit_plotters[4], best_fit_plotters[5])
generate_diameter_histogram()
generate_diameter_vs_albedo_plot()
generate_diameter_vs_period_plot()

print("Echoing relevant files\n")
os.system("echo 'PJDFC.out' | ../read-WISE-rc-MCMC-PJDFC") 
os.system("/bin/cp fort.2 Dhist.dat")
os.system("/bin/cp fort.32 Dhist_fine.dat")
os.system("/bin/cp fort.3 DvsPeriod.dat")
os.system("/bin/cp fort.4 DvsAlb.dat")