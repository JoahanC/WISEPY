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


def retrieve_output_data():
    
    output_file = open("PJDFC.out")
    output_file.readline()

    diameters = []
    chis = []
    gammas = []

    for line in output_file.readlines():
        datum = line.strip().split()
        if len(datum) < 10:
            #capture odd cases where col 6 was runing into col 5.  Only a handful, so just skip 
            continue
    
        diameters.append(np.e ** float(datum[5]))
        chis.append(float(datum[8]))
        gammas.append(np.e ** float(datum[4]))

    return diameters, chis, gammas


def retrieve_MCMC_data():
    #cp, not mv, so that there is a record of the actual outputs  
    os.system("/bin/cp ./fort.21 PJDFC.out")
    #fort.22 is the postscript format of the SED
    os.system("/bin/cp ./fort.22 SED_data.out")

    #^*^ parse and plot SED_data.out here 
    # best fit solution used for SED
    SED_file = open("SED_data.out")
    print("Reading in best fit solution used for SED.")
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
    first_line = True

    while True:
        line = SED_file.readline()
        
        if "SRGB" not in line and first_line:
            color.append("#444444")
            data_x.append([])
            data_y.append([])
            data_y_error.append([])

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
    print("\n*** Using "+cshfile+" to get MJDs. Make sure this is right*** \n")
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

    print("Generating best fit plot for SED.")
    
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

        # Plot flux for each wavelength
        plt.plot(wavelengths, esed[i], 
                color=colors[i % 10], 
                ls=linestyles[i % 3], 
                label=datelabels[i])    

        # Adjust y values and set up error ranges
        lower_errors = []
        higher_errors = []
        data_y_adjusted=[]
        for j in range(len(data_y[i])):

            lower_errors.append(data_y[i][j] * (1 - 1 / (10 ** (data_y_error[i][j] / 2.5))))
            higher_errors.append(data_y[i][j] * (10 ** (data_y_error[i][j] / 2.5) - 1))

            if data_y_error[i][j] > 0:    
                data_y_adjusted.append(data_y[i][j])

            else:
                # Measured flux is negative, so plot at an artifically
                # low point, and make error bars correct
                data_y_adjusted.append(1e-99)

        # Set scale and include error bars
        ax.set_xscale("log")
        ax.set_yscale("log")
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

    print("Generating diameter vs albedo plots.")
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

    print("Generating diameter histogram.")
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
    plt.xlabel("Diameter (km)")
    plt.ylabel("Number of Monte Carlo Results")
    plt.savefig("diameter_histogram.png")
    plt.savefig("diameter_histogram.pdf")


def generate_diameter_vs_period_plot():
    """
    Generates a plot of diameter and period solutions for the MCMC output.

    Arguments: None

    Returns: None, generates two plots titled diameter_vs_period.png/pdf
    """

    print("Generating diameter vs period plots.")
    diameters = []
    periods = []
    data_file = open("DvsPeriod.dat")
    for line in data_file.readlines():
        datum = line.rstrip().split()
        diameters.append(float(datum[0]))
        periods.append(float(datum[1]))

    if max(periods) - min(periods) == 0:
        print("Fixed period used, skipping diameter vs period plot")

    else:
        # Generate base figure
        plt.figure(figsize=[6,6])
        ax = plt.gca()
        plt.scatter(diameters, periods, s=1,marker='o',color='k')
        ax.set_xscale('log')
        ax.set_yscale('log')

        # Set up bins for contour levels
        H_bins, x_bins, y_bins = np.histogram2d(diameters, periods, bins=20)
        x_bins_adjusted = [(x_bins[i + 1] - x_bins[i]) / 2. + x_bins[i] for i in range(len(x_bins) - 1)]
        y_bins_adjusted = [(y_bins[i + 1] - y_bins[i]) / 2. + y_bins[i] for i in range(len(y_bins) - 1)]
        H_bins_adjusted  = np.transpose(H_bins)
        H_output = np.amax(H_bins_adjusted)
        
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

        # Add contour lines to plot
        plt.contour(x_bins_adjusted, y_bins_adjusted, H_bins_adjusted, 
                    colors="#cc0000", 
                    levels=[H_bins_flat[level_2], H_bins_flat[level_1]])
        plt.text(np.mean(diameters), 0.95 * max(periods), 
                "Contours contain 68%\nand 95.5% of all points")

        plt.xlabel("Diameter (km)")
        plt.ylabel("Rotation period")
        plt.savefig("diameter_vs_period.png")
        plt.savefig("diameter_vs_period.pdf")
        plt.close()


def generate_diameter_vs_gamma_plot(diameters, gammas):
    """
    Generates a plot of diameter and gamma solutions for the MCMC output.

    Arguments: diameters (list) -- a series of diameter solutions
               gammas (list) -- a series of gamma solutions

    Returns: None, generates two plots titled diameter_vs_gamma.png/pdf
    """
    print("Generating diameter vs gamma plots.")
    plt.figure(figsize=[6,6])
    ax = plt.gca()
    plt.scatter(diameters, gammas, s=1, marker='o', color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel("Diameter (km)")
    plt.ylabel("Gamma")
    plt.savefig("diameter_vs_gamma.png")
    plt.savefig("diameter_vs_gamma.pdf")
    plt.close()


def generate_diameter_vs_chi_plots(diameters, chis):
    """
    Generates a plot of diameter and chi solutions for the MCMC output.

    Arguments: diameters (list) -- a series of diameter solutions
               chis (list) -- a series of chi solutions

    Returns: None, generates two plots titled diameter_vs_chi.png/pdf
             and two zoomed plots titled diameter_vs_chi_zoom.png/pdf
    """

    print("Generating diameter vs chi plots.")
    plt.figure(figsize=[6, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.loglog(diameters, chis, color='k', marker='o', ms=0.5, ls="None")
    plt.xlabel("Diameter (km)", fontsize=13)
    plt.ylabel(r"fit $\chi^2$", fontsize=13)
    plt.savefig("diameter_vs_chi.png")
    plt.savefig("diameter_vs_chi.pdf")
    plt.close()

    plt.figure(figsize=[6, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.plot(diameters, chis, color='k', marker='o', ms=0.5, ls='None')
    plt.xlabel("Diameter (km)", fontsize=13)
    plt.ylabel(r"fit $\chi^2$", fontsize=13)
    plt.ylim(0.9 * min(chis), 2 * min(chis))  
    plt.savefig("diameter_vs_chi_zoom.png")
    plt.savefig("diameter_vs_chi_zoom.pdf")
    plt.close()

def generate_gamma_vs_chi_plots(gammas, chis):
    """
    Generates a plot of gamma and chi solutions for the MCMC output.

    Arguments: gammas (list) -- a series of gamma solutions
               chis (list) -- a series of chi solutions

    Returns: None, generates two plots titled gamma_vs_chi.png/pdf
    """

    print("Generating gamma vs chi plots.")
    plt.figure(figsize=[6, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.loglog(gammas, chis, color='k', marker='o', ms=0.5, ls='None')
    plt.xlabel(r"Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)", fontsize=13)
    plt.ylabel(r"fit $\chi^2$", fontsize=13)

    plt.savefig("gamma_vs_chi.png")
    plt.savefig("gamma_vs_chi.pdf")
    plt.close()


def display_MCMC_results():
    with open("best_fit.txt", 'r') as file:
        line_1 = file.readline().split()
        line_1[4] = "I=1 ... MNC ="
        output_1 = f"{line_1[0]} {line_1[1]} {line_1[2]} {line_1[3]}"
        output_1 += f" {line_1[4]} {line_1[5]} {line_1[6]}s"
        line_2 = file.readline().split()
        output_2 = f"Patch & Total Weights: {line_2[4]} & {line_2[5]}"
        print(output_2)
        line_3 = file.readline().strip()
        print(line_3)
        line_4 = file.readline().strip()
        print(line_4)
        file.readline()
        print("Out of 1 ... NMC loop\n\n*** Properties ***\n")
        line_6 = file.readline().split()
        output_6 = f"Diameter (Mean): {line_6[1][0:6]} +/- {line_6[2]}"
        print(output_6)
        output_7 = f"Diameter (Median): {line_6[4][0:6]} +{line_6[5][0:4]}%/-{line_6[6]}"
        print(output_7)
        line_8 = file.readline().split()
        output_8 = f"p_V (Mean): {line_8[2]} +/- {line_8[4]}"
        print(output_8)
        output_9 = f"p_V (Median): {line_8[6][0:6]} +/-n{line_8[8][0:4]}"
        print(output_9)
        line_10 = file.readline().split()
        output_10 = f"Theta_1 (Mean): {line_10[1][0:6]} +/- {line_10[2]}"
        print(output_10)
        output_11 = f"Theta_1 (Median): {line_10[4][0:6]} +{line_10[5][0:4]}%/-{line_10[6]}"
        print(output_11)
        line_12 = file.readline().split()
        if len(line_12) == 6:
            output_12 = f"Period (hour): {line_12[3][0:5]} +{line_12[4]}%/-{line_12[5]}"
        if len(line_12) == 5:
            output_12 = f"Period (hour): {line_12[3][0:4]} +{line_12[3][5:]}%/{line_12[4]}"
        print(output_12)
        line_13 = file.readline().split()
        if len(line_13) == 2:
            output_13 = f"Sqrt(Kappa*Rho*C): {line_13[1][0:5]} +{line_13[2]}&/{line_13[3]}"
        print(output_13)
        line_14 = file.readline().split()
        output_14 = f"Crater Fraction: {line_14[2][0:5]} +{line_14[2][6:11]}/-{line_14[2][12:]}"
        print(output_14)
        line_15 = file.readline().split()
        output_15 = f"p_IR/p_V: {line_15[1][0:5]} +{line_15[1][6:10]}%/-{line_15[1][11:]}"
        print(output_15)
        line_16 = file.readline().split()
        output_16 = f"Pole peak at: {line_16[4]} {line_16[5]}"
        print(output_16)
        line_17 = file.readline().split()
        output_17 = f"Mean pole at: {line_17[4]} {line_17[5]}"
        print(output_17)
        line_18 = file.readline().split()
        output_18 = f"Moment eigenvector: {line_18[2]} at RA,DEC: "
        output_18 += f"{line_18[5]} {line_18[6]}"
        print(output_18)
        line_19 = file.readline().split()
        output_19 = f"Moment eigenvector: {line_19[2]} at RA,DEC: "
        output_19 += f"{line_19[5]} {line_19[6]}"
        print(output_19)
        line_20 = file.readline().split()
        output_20 = f"Moment eigenvector: {line_20[2]} at RA,DEC: "
        output_20 += f"{line_20[5]} {line_20[6]}\n"
        print(output_20)


print("Echoing relevant files\n")
os.system("echo 'PJDFC.out' | ../read-WISE-rc-MCMC-PJDFC") 
os.system("/bin/cp fort.2 Dhist.dat")
os.system("/bin/cp fort.32 Dhist_fine.dat")
os.system("/bin/cp fort.3 DvsPeriod.dat")
os.system("/bin/cp fort.4 DvsAlb.dat")
#display_MCMC_results()

best_fit, best_fit_plotters, epoch_condition, wavelengths = retrieve_MCMC_data()
diameters, chis, gammas = retrieve_output_data()
generate_SED_plot(best_fit_plotters[0], best_fit_plotters[1], 
                  best_fit_plotters[2], best_fit_plotters[3], 
                  best_fit_plotters[4], best_fit_plotters[5])
generate_diameter_histogram()
generate_diameter_vs_albedo_plot()
generate_diameter_vs_period_plot()
generate_diameter_vs_gamma_plot(diameters, gammas)
generate_diameter_vs_chi_plots(diameters, chis)
generate_gamma_vs_chi_plots(gammas, chis)