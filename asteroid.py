from helpers import *
import os
import re
import warnings
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning
import matplotlib.pyplot as plt


class Asteroid:

    def __init__(self, unpacked_mpc_code, bands=[2]):
        """
        Constructor for the Asteroid class which represents an asteroid object with all
        associated MPC and WISE observations.
        
        Parameters
        ----------

        unpacked_mpc_code : int
            The unpacked MPC designation of the object.

        bands : list
            A list of the bands to be considered for this object. Each value must be an int.

        """
        self.bands = bands
        self.unpacked_mpc_code = unpacked_mpc_code
        self.packed_name = pack_MPC_name(unpacked_mpc_code)
        self.band_files = {}
        self.get_mpc_observations()
        self.split_wise_observations = {}
        self.total_wise_observations = []
        self.ra_dec_values = {}
        self.split_unique_observations = {}
        self.total_unique_observations = []
        self.source_ids = {}

        for band in bands:
            self.band_files[band] = return_input_files(self.packed_name, band)
            band_wise_observations = self.get_wise_observations(band)
            self.split_wise_observations[band] = band_wise_observations
            for observation in band_wise_observations:
                self.total_wise_observations.append(band_wise_observations[observation])
            self.ra_dec_values[band] = self.generate_ra_dec(self.band_files[band][1])
            unique_band_observations = self.generate_unique_observations(band)
            self.split_unique_observations[band] = unique_band_observations
            for observation in unique_band_observations:
                self.total_unique_observations.append(unique_band_observations[observation])
        
        for band in self.split_unique_observations:
            self.source_ids[band] = []
            for observation in self.split_unique_observations[band]:
                self.source_ids[band].append(self.split_unique_observations[band][observation][0][:9])

    
    def run_comparison(self):
        for band in self.bands:
            self.generate_new_table(band)
            self.generate_full_table(band)
        self.generate_flux_snr_plots()


    def get_mpc_observations(self):
        """
        Populates the Asteroid object with all MPC observations from a .txt in the 
        /input_data/ directory. The file must follow the following naming scheme:

        ./input_data/[[packed_MPC_name]].txt
        
        Parameters
        ----------
        None


        Returns
        -------
        None, populates the self.mpc_observations field with a dictionary containing the
        utc dates of each observation as the key while the observation index[0] and 
        julian dates[1] form a 2 element list for each key.
        """
        dates = []
        ids = []
        obs_ids = []

        with open(return_input_files(self.packed_name)[0], 'r') as mpc_file:
            for idx, row in enumerate(mpc_file.readlines()):
                row_string = row
                id_string = row_string[:15] 
                ids.append(id_string)
                row_string = row_string[15:]
                year = row_string[:4]
                month = row_string[5:7]
                row_string = row_string[8:]
                elements = row_string.split()
                day = elements[0]
                dec_day = day[2:9]
                time_string = decimal_day_converter(dec_day)
                day = day[:2]
                date = year+'-'+month+'-'+day+time_string
                dates.append(date[:23])
                observation_id = row_string[47:]
                obs_ids.append(observation_id)
        
        utc_dates = np.array(dates)
        time_object = Time(utc_dates, format='isot')
        julian_dates = time_object.jd
        observations = {}

        for idx, date in enumerate(utc_dates):
            observations[date] = [obs_ids[idx], julian_dates[idx]]
        self.mpc_observations = observations

    
    def generate_ra_dec(self, irsa_file):
        """
        Returns all the RA DEC values for the WISE observations of this object
        of a given band.

        Parameters
        ----------
        irsa_file : str
            The Astropy table file from IRSA containing all of the observation
            data for a given band.

        Returns
        -------
        A dictionary of source_ids with the ids being the keys, and the values
        being a two element list of RA values[0] and DEC values[1].
        """
        data_object = Table.read(irsa_file, format='ipac')
        source_ids = list(data_object['source_id'])
        ra = list(data_object['ra'])
        dec = list(data_object['dec'])
        info = {}
        for index, source_id in enumerate(source_ids):
            info[source_id[:9]] = [ra[index], dec[index]]
        return info


    def pull_fluxes(self, data_object, bands):
        """
        Pulls all flux values and flux uncertainties for a given object.

        Parameters
        ----------
        data_object : Table
            The data for a given object.
        
        bands : int
            The number of wavelength bands to be looked at.
        
        Returns
        -------
            All flux values and flux uncertainties for the number of bands
            requested as an n-element tuple.
        """
        w1_flux = list(data_object['w1flux'])
        w1_flux_sigma = list(data_object['w1sigflux'])
        w2_flux = list(data_object['w2flux'])
        w2_flux_sigma = list(data_object['w2sigflux'])
        if bands == 2:
            return (w1_flux, w1_flux_sigma,
                w2_flux, w2_flux_sigma)
        if bands == 3:
            w3_flux = list(data_object['w3flux'])
            w3_flux_sigma = list(data_object['w3sigflux'])
            return (w1_flux, w1_flux_sigma,
                w2_flux, w2_flux_sigma,
                w3_flux, w3_flux_sigma)
        if bands == 4:
            w3_flux = list(data_object['w3flux'])
            w3_flux_sigma = list(data_object['w3sigflux'])
            w4_flux = list(data_object['w4flux'])
            w4_flux_sigma = list(data_object['w4sigflux'])
            return (w1_flux, w1_flux_sigma,
                w2_flux, w2_flux_sigma,
                w3_flux, w3_flux_sigma,
                w4_flux, w4_flux_sigma)


    def get_wise_observations(self, bands=2):
        """
        Returns WISE observation data for a given band.
        
        Parameters
        ----------
        band : int 
            The wavelength band being queried.
        
        Returns
        -------
        A dictionary where the keys correspond to the utc dates for each image
        and an n-element list with source ids[0], julian dates[1], RA[2], DEC[3],
        and corresponding flux and flux uncertainty measurements[4+].
        """

        data_object = Table.read(return_input_files(self.packed_name, bands)[1], format='ipac')
        dates = list(data_object['mjd'])
        source_ids = list(data_object['source_id'])
        ra = list(data_object['ra'])
        dec = list(data_object['dec'])
        fluxes = self.pull_fluxes(data_object, bands)
        
        time_object = Time(dates, format='mjd', scale='utc')
        julian_dates = time_object.jd
        utc_dates = time_object.isot

        images = {}
        for idx, date in enumerate(utc_dates):
            if bands == 2:
                images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                    ra[idx], dec[idx], 
                                    fluxes[0][idx], fluxes[1][idx], 
                                    fluxes[2][idx], fluxes[3][idx]]
            if bands == 3:
                images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                    ra[idx], dec[idx], 
                                    fluxes[0][idx], fluxes[1][idx], 
                                    fluxes[2][idx], fluxes[3][idx],
                                    fluxes[4][idx], fluxes[5][idx]]
            if bands == 4:
                images[str(date)] = [source_ids[idx], julian_dates[idx], 
                                    ra[idx], dec[idx], 
                                    fluxes[0][idx], fluxes[1][idx], 
                                    fluxes[2][idx], fluxes[3][idx],
                                    fluxes[4][idx], fluxes[5][idx],
                                    fluxes[6][idx], fluxes[7][idx]]
        return images


    def generate_unique_observations(self, band):
        """
        Compares the observational instances between the MPC and WISE dataset 
        files for a given object and populates the objects 
        total_unique_observations and split_unique_observations fields for
        a given band.
        
        Parameters
        ----------
        band : int
            The band being queried.
        
        Returns
        -------
        None, populates the object's fields with a subset of wise observations.
        """
        # Year Counting
        self.wise_years = []
        self.mpc_years = []
        for datum in self.split_wise_observations[band]:
            year = int(datum[:4])
            if year not in self.wise_years:
                self.wise_years.append(year)
        for datum in self.mpc_observations:
            year = int(datum[:4])
            if year not in self.mpc_years:
                self.mpc_years.append(year)
        
        wise_years_str = ""
        for year in self.wise_years:
            wise_years_str += str(year) + ", "
        
        mpc_years_str = ""
        for year in self.mpc_years:
            mpc_years_str += str(year) + ", "
        
        # Acquire utc dates
        wise_utc = list(self.split_wise_observations[band].keys())
        mpc_utc = list(self.mpc_observations.keys())

        mpc_intervals = {}
        utc_delta = TimeDelta("11.0", format='sec')
        for datum in mpc_utc:
            base_time = Time(datum, format='isot')
            lower_bound = base_time - utc_delta
            upper_bound = base_time + utc_delta
            mpc_intervals[datum] = [lower_bound.jd, upper_bound.jd]

        wise_intervals = {}
        for datum in wise_utc:
            base_time = Time(datum, format='isot')
            wise_intervals[datum] = [base_time.jd]

        new_epochs = {}
        for epoch in wise_intervals:
            recorded = False
            for observation in mpc_intervals:
                if wise_intervals[epoch] >= mpc_intervals[observation][0] and wise_intervals[epoch] <= mpc_intervals[observation][1]:
                    recorded = True
            if not recorded:
                new_epochs[epoch] = self.split_wise_observations[band][epoch]

        return new_epochs


    def generate_new_table(self, band):
        """
        Generates an ipac format table from an existing WISE table for a set of
        new images documented as a series of source ids in a /loader_data/ file.
        The /loader_data/ file for this object must fit the following naming 
        scheme:

        /loader_data/[[packed_mpc_name]][[.tbl/_3band.tbl/_cryo.tbl]]

        Parameters
        ----------
        band : int
            The band being queried.

        Returns
        -------
        None, generates a new .tbl file of ipac format in /observations/new/[[packed_mpc_name]]/
        """
        
        warnings.simplefilter('ignore', category=AstropyUserWarning)

        # Reads in all WISE source ids for a given asteroid and band set
        wise_file = f"database_files/wise/{self.packed_name}_{band}band.tbl"
        data_object = Table.read(wise_file, format='ipac')
        wise_sids = list(data_object['source_id'])
        wise_sids = [sid[0:9] for sid in wise_sids]

        # Reads in unique source ids
        unique_sids = []
        with open(f"loader_data/{self.packed_name}_{band}band.txt", 'r') as file:
            for line in file:
                unique_sids.append(line.rstrip())
        # Generates a mask with unique source ids and writes tbl file
        mask = np.isin(wise_sids, unique_sids)
        t_new = data_object[mask]
        t_new.write(f"observations/new/{self.packed_name}_{band}bands.tbl", 
                format="ipac", overwrite=True)


    def return_new_sids(self, band):
        """
        Returns a list of all the newly recovered source ids located in
        /observations/new/.
        
        Parameters
        ----------
        band : int
            The band being queried.

        Returns
        -------
        A list of all new source ids.
        """
        
        warnings.simplefilter('ignore', category=AstropyUserWarning)

        # Reads in all WISE source ids for a given asteroid and band set
        new_sids = []
        try:
            wise_file = f"observations/new/{self.packed_name}_{band}bands.tbl"
            data_object = Table.read(wise_file, format='ipac')
            sids = list(data_object['source_id'])
            sids = [sid[0:9] for sid in sids]
            new_sids.extend(sids)
        except:
            print('')
        return new_sids


    def generate_full_table(self, band):
        """
        Generates a .tbl file containing all useful WISE observational instances
        which are useful for thermophysical modeling.

        Parameters
        ----------
        band : int
            The band being queried.

        Returns
        -------
        None, generates a new .tbl file of ipac format in /observations/all/[[packed_mpc_name]]/
        """
        warnings.simplefilter('ignore', category=AstropyUserWarning)

        # Reads in all known WISE observations
        wise_obs = {}
        for date in self.split_wise_observations[band]:
            wise_obs[date] = self.split_wise_observations[band][date]

        # Read in all known MPC observations
        reported_obs = {}
        for date in self.mpc_observations:
            if self.mpc_observations[date][0][7:] == "C51":
                reported_obs[date] = self.mpc_observations[date]

        # Acquire utc dates
        wise_utc = list(wise_obs.keys())
        mpc_utc = list(self.mpc_observations.keys())
        
        mpc_intervals = {}
        utc_delta = TimeDelta("11.0", format='sec')
        for datum in mpc_utc:
            base_time = Time(datum, format='isot')
            lower_bound = base_time - utc_delta
            upper_bound = base_time + utc_delta
            mpc_intervals[datum] = [lower_bound.jd, upper_bound.jd]

        wise_intervals = {}
        for datum in wise_utc:
            base_time = Time(datum, format='isot')
            wise_intervals[datum] = [base_time.jd]

        existing_epochs = {}
        for epoch in wise_intervals:
            recorded = False
            for observation in mpc_intervals:
                if wise_intervals[epoch] >= mpc_intervals[observation][0] and wise_intervals[epoch] <= mpc_intervals[observation][1]:
                    recorded = True
            if recorded:
                existing_epochs[epoch] = wise_obs[epoch]

        existing_sids = []
        for epoch in existing_epochs:
            existing_sids.append(existing_epochs[epoch][0][:9])

        new_sids = self.return_new_sids(band)
        all_sids = existing_sids + new_sids

        data_object = Table.read(self.band_files[band][1], format='ipac')
        wise_sids = list(data_object['source_id'])
        wise_sids = [sid[0:9] for sid in wise_sids]
        mask = np.isin(wise_sids, all_sids)
        t_new = data_object[mask]
        t_new.write(f"observations/all/{self.packed_name}_{band}bands.tbl", 
                format="ipac", overwrite=True)


    def display_unique_observations(self):
        """
        Displays the number of epochs cataloged in the Asteroid object as well as 
        the unique years in which data is present.
        """
        wise_years_str = ""
        for year in self.wise_years:
            wise_years_str += str(year) + ", "
        
        mpc_years_str = ""
        for year in self.mpc_years:
            mpc_years_str += str(year) + ", "

        print("Epochs observed in the MPC database:", len(self.mpc_observations))
        print("Epochs observed in the WISE database:", len(self.total_wise_observations))
        print("Years in which MPC data was collected: " + mpc_years_str[:-2])
        print("Years in which WISE data was collected: " + wise_years_str[:-2])
        print(f"New epochs detected in WISE catalog: {len(self.total_unique_observations)}")

    
    def generate_flux_snr_plots(self):
        """
        Generates a series of flux and SNR plots with associated error bars
        for a given MPC object and its targeted bandset.

        Parameters
        ----------
        bands : int
            The wavelength band being plotted.

        Returns
        -------
        A series of plots in /plots/flux_plots/ and /plots/snr_plots/.
        """
        if self.packed_name not in os.listdir("./plots/flux_plots/"):
            os.mkdir(f"./plots/flux_plots/{self.packed_name}")
        if self.packed_name not in os.listdir("./plots/snr_plots/"):
            os.mkdir(f"./plots/snr_plots/{self.packed_name}")
        for band in self.bands:
            if f"{band}_band" not in os.listdir(f"./plots/flux_plots/{self.packed_name}"):
                os.mkdir(f"./plots/flux_plots/{self.packed_name}/{band}_band")
            if f"{band}_band" not in os.listdir(f"./plots/snr_plots/{self.packed_name}"):
                os.mkdir(f"./plots/snr_plots/{self.packed_name}/{band}_band")
            
            new_data = Table.read(f"observations/new/{self.packed_name}_{band}bands.tbl", format="ipac")
            all_data = Table.read(f"observations/all/{self.packed_name}_{band}bands.tbl", format="ipac")
            mjd_new = list(new_data["mjd"])
            mjd_all = list(all_data["mjd"])
            
            new_flux_values = {}
            for i in range(band):
                new_flux_values[i + 1] = [list(new_data[f"w{i + 1}flux"]), list(new_data[f"w{i + 1}sigflux"])]
            for set in new_flux_values:
                template_new_plot(self.packed_name, mjd_new, new_flux_values[set], "Flux", set, band)
            
            new_snr_values = {}
            for i in range(band):
                new_snr_values[i + 1] = [list(new_data[f"w{i + 1}snr"]), np.zeros(len(list(new_data[f"w{i + 1}snr"])))]
            for set in new_snr_values:
                template_new_plot(self.packed_name, mjd_new, new_snr_values[set], "SNR", set, band)
            
            all_flux_values = {}
            for i in range(band):
                all_flux_values[i + 1] = [list(all_data[f"w{i + 1}flux"]), list(all_data[f"w{i + 1}sigflux"])]
            for set in all_flux_values:
                template_composite_plot(self.packed_name, mjd_new, mjd_all, new_flux_values[set], all_flux_values[set], "Flux", set, band)

            all_snr_values = {}
            for i in range(band):
                all_snr_values[i + 1] = [list(all_data[f"w{i + 1}snr"]), np.zeros(len(list(all_data[f"w{i + 1}snr"])))]
            for set in all_snr_values:
                template_composite_plot(self.packed_name, mjd_new, mjd_all, new_snr_values[set], all_snr_values[set], "SNR", set, band)


    def data_sort(self, source_ids, bands=2):
        """
        The sorting algorithm for this object's unique source ids.
        
        Parameters
        ----------
        source_ids : list
            A list of all unique epoch source ids.

        bands: int
            The band set being sorted.
        
        Returns
        A sorted list of .FITS files corresponding to the source ids provided.
        """
        
        # Generate all relevant keys and perform initial sort
        files = os.listdir("wise_images/" + str(self.unpacked_mpc_code) + '/')
        wise_files = []
        for file in files:
            if file[:9] in source_ids:
                wise_files.append(file)


        file_id = {}
        for file in wise_files:
            if file == ".DS_Store" or file == ".DS_S":
                pass
            if int(file[:5]) in file_id.keys():
                file_id[int(file[:5])].append(file)
            else: 
                file_id[int(file[:5])] = [file]

        key_sort = list(file_id.keys())
        key_sort.sort()

        # Sky region sort
        number_sorted = {}
        for key in key_sort:
            for file in wise_files:
                code = int(file[0:5])
                if code == key and code in number_sorted.keys():
                    number_sorted[code].append(file)
                elif code == key and code not in number_sorted.keys():
                    number_sorted[code] = [file]
                else:
                    pass
        
        for number in number_sorted:
            regex = r"[0-9]{5}[a]"
            matches = []
            for test in number_sorted[number]:
                if len(re.findall(regex, test)) == 0:
                    pass
                else:
                    matches.append(test)
            regex = r"[0-9]{5}[b]"
            for test in number_sorted[number]:
                if len(re.findall(regex, test)) == 0:
                    pass
                else:
                    matches.append(test)
            regex = r"[0-9]{5}[c]"
            for test in number_sorted[number]:
                if len(re.findall(regex, test)) == 0:
                    pass
                else:
                    matches.append(test)
            regex = r"[0-9]{5}[r]"
            for test in number_sorted[number]:
                if len(re.findall(regex, test)) == 0:
                    pass
                else:
                    matches.append(test)
            regex = r"[0-9]{5}[s]"
            for test in number_sorted[number]:
                if len(re.findall(regex, test)) == 0:
                    pass
                else:
                    matches.append(test)

            number_sorted[number] = matches

        sorted_run = []
        for number in number_sorted:
            a_temp = {}
            b_temp = {}
            for elem in number_sorted[number]:
                if elem[5] == 'a':
                    if elem[6:9] not in a_temp.keys():
                        a_temp[elem[6:9]] = [elem]
                    else:
                        a_temp[elem[6:9]].append(elem)
                else:
                    if elem[6:9] not in b_temp.keys():
                        b_temp[elem[6:9]] = [elem]
                    else:
                        b_temp[elem[6:9]].append(elem)
            a_keys = list(a_temp.keys())
            b_keys = list(b_temp.keys())
            a_keys.sort()
            b_keys.sort()
            for key in a_keys:
                for elem in a_temp[key]:
                    sorted_run.append(elem)
            for key in b_keys:
                for elem in b_temp[key]:
                    sorted_run.append(elem)

        #for file in sorted_run:
        #    print(file)


        # Test W band sort
        idx = list(range(len(sorted_run)))[::bands]
        w_sorted = []
        interval = {}
        for i in idx:
            for j in range(bands):
                interval[int(sorted_run[i + j][11])] = sorted_run[i + j]
            for k in range(1, bands + 1):
                w_sorted.append(interval[k])
        renamed_sorted_run = []
        for file in w_sorted:
            renamed_sorted_run.append("wise_images/" + str(self.unpacked_mpc_code) + '/' + file)
        return renamed_sorted_run

    
    def make_region(self, fits_file, source_ids):
        """
        Makes a region file for a given .FITS file in the /regions/ folder.
        
        Parameter
        ---------
        file : str
            A .FITS file located in /wise_images/[[packed_mpc_name]].

        source_ids : dict
            A dictionary containing RA[0] and DEC[1] values for each wise
            observation, represented by the keys consiting of source ids. 
        
        Returns
        -------
        None, generates a corresponding .reg file in the /regions/ folder.
        """
        sid = fits_file[13 + len(str(self.unpacked_mpc_code)):22 + len(str(self.unpacked_mpc_code))]
        band = fits_file[23 + len(str(self.unpacked_mpc_code)):25 + len(str(self.unpacked_mpc_code))]
        with open(f"regions/{sid}_{band}.reg", 'w') as file:
            file.write("# Region file format: DS9 version 4.1\n")
            file.write('global color=green dashlist=8 3 width=1 ' 
                    + 'font="helvetica 10 normal roman" select=1 ' 
                    + 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 ' 
                    + 'include=1 source=1\n')
            file.write("fk5\n")
            file.write(f'circle({source_ids[sid][0]},{source_ids[sid][1]},15.000")')

    
    def generate_script(self, band, lower_bound, upper_bound):
        """
        Generates a ds9 terminal script which can be used to view a set of ds9 images 
        for the unique observations of this object.
        
        Parameters
        ----------

        band : int
            The specific band set of images to look at.

        lower_bound : int
            The lower index to begin loading images. This method uses zero-indexing.

        upper_bound : int
            The upper inex to stop loading images.
        
        Returns
        -------
        A string which can be run in the terminal to view a subset of ds9 images.
        """
        print(self.source_ids[band])
        sorted_files = self.data_sort(self.source_ids[band], band)
        print(len(sorted_files))
        region_files = os.listdir("regions")
        for file in region_files:
            os.remove("regions/" + file)

        lookup_table = self.generate_ra_dec(self.band_files[band][1])

        for file in sorted_files:
            sid = file[13 + len(str(self.unpacked_mpc_code)):22 + len(str(self.unpacked_mpc_code))]
            self.make_region(file, lookup_table)

        file_region = {}
        for file in sorted_files:
            sid = file[13 + len(str(self.unpacked_mpc_code)):22 + len(str(self.unpacked_mpc_code))]
            band = file[23 + len(str(self.unpacked_mpc_code)):25 + len(str(self.unpacked_mpc_code))]
            file_region[file] = f"regions/{sid}_{band}.reg"

        run_string = "ds9 -scale log -tile "
        for file in list(file_region.keys())[lower_bound:upper_bound]:
            run_string += file + ' -regions '
            reg_string = file_region[file]
            run_string += reg_string + ' '
        run_string += ' -zmax'
        os.popen(run_string)
        #return run_string


    def load_files(self, load_file, band=2):
        """
        Generates a ds9 terminal script which can be used to view a set of ds9 images 
        for source ids recorded in a .txt in /loader_files/
        
        Parameters
        ----------

        load_file : str
            The name of a .txt file found in the /loader_data/ folder containing source id
            stubs.

        band : int
            The specific band set of images to look at.
        
        Returns
        -------
        A string which can be run in the terminal to view a subset of ds9 images.
        """
        file_stubs = []
        with open(f"loader_data/{load_file}", 'r') as file:
            for line in file:
                file_stubs.append(line.rstrip())

        wise_files = os.listdir(f"wise_images/{self.unpacked_mpc_code}")
        sorted_run = []
        for stub in file_stubs:
            for file in wise_files:
                if file[:9] == stub:
                    sorted_run.append(file)

        if band == 4:
            idx = list(range(len(sorted_run)))[::4]
            w_sorted = []
            for i in idx:
                quartet = {int(sorted_run[i][11]) : sorted_run[i], 
                int(sorted_run[i+1][11]): sorted_run[i+1], 
                int(sorted_run[i+2][11]): sorted_run[i+2], 
                int(sorted_run[i+3][11]): sorted_run[i+3]}
                for i in range(1, 5):
                    w_sorted.append(quartet[i])

        if band == 3:
            idx = list(range(len(sorted_run)))[::3]
            w_sorted = []
            for i in idx:
                triplet = {int(sorted_run[i][11]) : sorted_run[i], 
                int(sorted_run[i+1][11]): sorted_run[i+1], 
                int(sorted_run[i+2][11]): sorted_run[i+2]}
                for i in range(1, 4):
                    w_sorted.append(triplet[i])

        if band == 2:
            idx = list(range(len(sorted_run)))[::2]
            w_sorted = []
            for i in idx:
                pair = {int(sorted_run[i][11]) : sorted_run[i], 
                int(sorted_run[i+1][11]): sorted_run[i+1]}
                for i in range(1, 3):
                    w_sorted.append(pair[i])

        band_file_map = {2: ("input_data/" + str(self.packed_name) + ".txt", 
                        "input_data/" + str(self.packed_name) + ".tbl", 0),
                        3: ("input_data/" + str(self.packed_name) + ".txt",
                        "input_data/" + str(self.packed_name) + "_3band.tbl", 22),
                        4: ("input_data/" + str(self.packed_name) + ".txt",
                        "input_data/" + str(self.packed_name) + "_cryo.tbl", 44)}

        mpc_file, wise_file = band_file_map[band][0], band_file_map[band][1]

        new_epochs = self.split_unique_observations[band]
        good_epochs = {}
        for epoch in new_epochs:
            sid = new_epochs[epoch][0][:9]
            for stub in file_stubs:
                if stub == sid:
                    good_epochs[sid] = new_epochs[epoch]
                    good_epochs[sid][0] = epoch    
    
        #terminal_table(mpc_file, bands, good_epochs)
        
        region_list = []
        for file in w_sorted:
            region_list.append(f"wise_images/{self.unpacked_mpc_code}/" + file)

        region_files = os.listdir("regions")
        for file in region_files:
            os.remove("regions/" + file)
        lookup_table = {}
        for id in self.ra_dec_values[band]:
            print('here', id)
            if id in file_stubs:
                lookup_table[id] = self.ra_dec_values[band][id]
        print(lookup_table)
        print(region_list)
        for file in region_list:
            self.make_region(file, lookup_table)

        ds9_script = "ds9 -tile "
        for file in w_sorted:
            sid, w_band = file[:9], file[10:12]
            ds9_script += f"wise_images/{self.unpacked_mpc_code}/" + file + ' -regions '
            ds9_script += f"regions/{sid}_{w_band}.reg "
        os.popen(ds9_script)