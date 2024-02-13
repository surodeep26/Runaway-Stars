import os
from astropy.table import Table, Column, QTable, join
import yaml
import pandas as pd  # Renamed for clarity
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from astroquery.vizier import Vizier
from scipy.stats import norm

from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_sun, get_body
from IPython.display import Math
from astropy.wcs import WCS
from astroquery.skyview import SkyView
import os
import numpy as np
# from runaway_functions import *
import pandas as pd
import shutil


from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import matplotlib.patheffects as path_effects
import matplotlib.patches as patches
from astropy.visualization.wcsaxes import add_scalebar
from astropy.coordinates import Angle
from psrqpy import QueryATNF

from astroquery.simbad import Simbad
import astropy.coordinates as coord


from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
from astropy import units as u
from IPython.display import display, Math
import yaml
from astropy.stats import sigma_clip



cluster_list = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
def read_yaml_file(file_path):
    '''
    Read the configuration file for the rest of the code. 
    This contains the various parameters for the code to run.
    '''
    with open(file_path, 'r') as yaml_file:
        config = yaml.safe_load(yaml_file)
    return config
config = read_yaml_file('config.yaml')

class Cluster:
    def __init__(self, cluster_name):
        # Assuming 'cluster_list' is a predefined list containing your Astropy table
        self.cluster_table = cluster_list[cluster_list['Cluster'] == cluster_name]
        
        if len(self.cluster_table) == 0:
            raise ValueError(f"No data found for cluster '{cluster_name}' in the cluster list.")
        
        # Extracting values from the table and setting attributes
        self.name = self.cluster_table['Cluster'][0]
        if not os.path.exists(f'{self.name}'):
            os.mkdir(f'{self.name}')

        # self.ra_icrs = self.cluster_table['RA_ICRS'][0]
        # self.de_icrs = self.cluster_table['DE_ICRS'][0]
        ti = (Time('J2000')+1*u.Myr)
        self.coordinates = SkyCoord(ra=self.cluster_table['RA_ICRS'],dec=self.cluster_table['DE_ICRS'],pm_ra_cosdec=self.cluster_table['pmRA'],pm_dec=self.cluster_table['pmDE'],obstime=ti)
        self.diameter_dias = (self.cluster_table['r50'][0]*u.deg*2).to(u.arcmin)
        self.diameter = self.cluster_table['Diameter'][0]*u.arcmin
        self.N = self.cluster_table['N'][0] #no. of members
        self.distance = self.cluster_table['Dist'][0]*u.pc
        self.RV = self.cluster_table['RV'][0]*u.km/u.s
        self.all = self.cluster_table[0]

    def calculate_search_arcmin(self, extra=config['Cluster']['search_extent'], output=False):
        """
        Calculate the search arcminute for the cluster.

        Parameters:
        extra (float, optional): Extent of search around cluster in pc. Default is taken from config file `search_extent`.
        output (bool, optional): Whether to display the output of some useful parameters. Default is False.

        Returns:
        Quantity: The search arcminute value.
        """
        theta = self.diameter / 2  # radius of the cluster in arcmin
        D = self.distance  # Assuming you have a 'distance' attribute in your Cluster class
        r = np.tan(theta) * D
        # display(r)
        search_arcmin = np.arctan((r + extra * u.pc) / D)
        search_arcmin = search_arcmin.to(u.arcminute)
        r, search_arcmin = r.round(3), search_arcmin.round(3)
        
        if output:
            display(Math(r'\mathrm{Dist.}' + f'= {D}' + r'\mathrm{\ Ang. \ Radius}' + f'= {theta}' +
                        r'\mathrm{\ Phy. \ Radius}' + f'= {r}'))
            display(Math(r'\mathrm{Extra \ } '+f'= {extra} pc'))
            display(Math(r'\mathrm{\ search \ arcmin}' + f'= {search_arcmin}'))

        return search_arcmin

    def plot_search_region(self, extra=config['Cluster']['search_extent']):
        """
        Plots and saves the fits file for the region around the given cluster.

        Parameters:
        - extra (float): Additional extent for the search region (default is from config['Cluster']['search_extent']).

        Returns:
        None
        """
        search_arcmin = self.calculate_search_arcmin(extra=extra)
        
        # Define the file path
        fits_file_path = f'{self.name}/{self.name}_extra{extra}pc.fits'
        
        # Check if the file already exists
        if os.path.exists(fits_file_path):
            print(f'fits image exists in {self.name} folder')
            # File exists, no need to download, use the existing file
            images = [fits.open(fits_file_path)]
            # Extract the WCS information
            wcs = WCS(images[0][0].header)

        else:
            # File doesn't exist, get the image data from SkyView
            images = SkyView.get_images(position=self.coordinates,
                                        survey=config['Plots']['skyview_survey'],
                                        radius=2*search_arcmin)
            # Extract the WCS information
            wcs = WCS(images[0][0].header)
            hdu = fits.PrimaryHDU(data=images[0][0].data, header=images[0][0].header)
            hdulist = fits.HDUList([hdu])
            
            # Save the fits file
            hdulist.writeto(fits_file_path, overwrite=True)
        # Plot the image
        print(f'Plotting a {extra} pc region around cluster.')
            
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection': wcs})
        ax.imshow(images[0][0].data, cmap='gray')
        data = images[0][0].data
        # Add labels and title
        ax.set_xlabel('Right Ascension (degrees)')
        ax.set_ylabel('Declination (degrees)')
        ax.set_title(f"{self.name} with {extra}pc search region")
        ax.grid(color='lightgrey',ls='dotted')

        # Add a circle representing the cluster radius
        circle_cluster = plt.Circle((data.shape[1] / 2, data.shape[0] / 2), radius=150*self.diameter/search_arcmin,
                        edgecolor='red', facecolor='none',ls='dashed',label=f'Cluster Diameter = {self.diameter}')
        circle_search_region = plt.Circle((data.shape[1] / 2, data.shape[0] / 2), radius=150*search_arcmin/search_arcmin,
                    edgecolor='green', facecolor='none',ls='dashed',label=f'Search Region = {search_arcmin}')
        ax.add_artist(circle_cluster)
        ax.add_artist(circle_search_region)
        _ = search_arcmin.round(-1)/2 #round to the nearest 5 (by round first to nearest 10 and then divide by 2
        scalebar_length = ((_.to(u.rad))*(self.distance.to(u.m)).to(u.pc)).round(2)
        __ = self.distance.value/1000
        add_scalebar(ax, _,color="yellow",label=f"or {scalebar_length.value}pc (at dist {__:.2f}kpc)",size_vertical=0.5)
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        ax.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontsize=12,fontweight='bold')
        ax.legend()
        # Show the plot
        plt.show()

    def dias_members(self,memb_prob=config['memb_prob'],parallax_quality_threshold=config['parallax_quality_threshold']):
        dias_members = Table.read(f'Clusters/{self.name}.dat', format='ascii.tab')
        dias_members.remove_rows([0,1])
        # Get the column names in the table
        column_names = dias_members.colnames
        column_names.remove('Source')
        # Convert all columns from string to float, replacing non-float values with None
        for column_name in column_names:
            column_data = dias_members[column_name]
            converted_data = []
            for entry in column_data:
                try:
                    converted_data.append(float(entry))
                except ValueError:
                    converted_data.append(None)
            dias_members[column_name] = converted_data
        dias_members.sort('Source')

        mask1 = [value >= memb_prob for value in dias_members['Pmemb']]
        mask2 = [value > parallax_quality_threshold for value in dias_members['Plx']/dias_members['e_Plx']]
        final_mask = (np.array(mask1) & np.array(mask2))
        dias_members = dias_members[final_mask]
        #create a filter to get rid of the stars which don't have a Bp-Rp value in the .dat table
        mask_bprp = [value is not None for value in dias_members["BP-RP"] ]
        plottable_members = dias_members[mask_bprp]
        return plottable_members
    

    def stars_in_region(self,allTables=False):
        """
        Finds the stars around the cluster within a conical frustum.
        Saves the three tables in the cluster folder.
        Qurries VizieR only if tables are not already present in the cluster folder.

        Parameters:
        allTables (bool, optional): Default is false, and returns 1 argument - stars_in_region
        which is the stars_fromDR3 and stars_fromDR3_dis joined.

        Returns:
        stars_in_region (astropy table)

        or (stars_in_region, stars_fromDR3, stars_fromDR3_dis) if allTables is True
               
        """
        if not os.path.exists(f'{self.name}'):
            os.mkdir(f'{self.name}')

        stars_fromDR3_path = f'{self.name}/{self.name}_stars_fromDR3.tsv'
        stars_fromDR3_dis_path = f'{self.name}/{self.name}_stars_fromDR3_dis.tsv'
        stars_in_region_path =  f'{self.name}/{self.name}_stars_in_region.tsv'
        fs_path =  f'{self.name}/{self.name}_fs.tsv'
        nsfs_path =  f'{self.name}/{self.name}_nsfs.tsv'

        if allTables and (os.path.exists(stars_in_region_path)) and os.path.exists(stars_fromDR3_path) and os.path.exists(stars_fromDR3_dis_path) and os.path.exists(fs_path) and os.path.exists(nsfs_path):
            #checks for generate_tables() if all the tables are already present.
            stars_in_region = Table.read(stars_in_region_path,format='ascii.ecsv')
            stars_fromDR3 = Table.read(stars_fromDR3_path,format='ascii.ecsv')
            stars_fromDR3_dis = Table.read(stars_fromDR3_dis_path,format='ascii.ecsv')
            fs = Table.read(fs_path,format='ascii.ecsv')
            nsfs = Table.read(nsfs_path,format='ascii.ecsv')

            print("All tables present (stars_in_region,stars_fromDR3,stars_fromDR3_dis,fs,nsfs)")
            return stars_in_region,stars_fromDR3,stars_fromDR3_dis,fs,nsfs


        if os.path.exists(stars_in_region_path) and (not allTables):
            stars_in_region = Table.read(stars_in_region_path,format='ascii.ecsv')
            print(f'{stars_in_region_path} exists in {self.name} with {len(stars_in_region)} stars')

            return stars_in_region


        if os.path.exists(stars_fromDR3_path):
            stars_fromDR3 = Table.read(stars_fromDR3_path,format='ascii.ecsv')
            print(f'Table stars_fromDR3 exists in {self.name} with {len(stars_fromDR3)} stars, skipping VizieR query')

        else:
            print(f'Querying VizieR {config["cat_stars"]} \n at: {self.coordinates.ra}, {self.coordinates.dec} \n within: {self.calculate_search_arcmin()}')
            print(f'Filters: \n {config["cat_stars_filters"]}')
            stars_fromDR3 = Vizier(columns=["*","+_r"],row_limit = -1).query_region(self.coordinates, 
                                                                                       radius=self.calculate_search_arcmin(), 
                                                                                       catalog=config['cat_stars'],
                                                                                       column_filters=config['cat_stars_filters'])[0]
            print(f'{len(stars_fromDR3)} sources found')
        if os.path.exists(stars_fromDR3_dis_path):
            stars_fromDR3_dis = Table.read(stars_fromDR3_dis_path,format='ascii.ecsv')
            print(f'Table stars_fromDR3_dis exists in {self.name} with {len(stars_fromDR3_dis)} stars, skipping VizieR query')


        else:
            print(f'Querying VizieR {config["cat_star_distances"]} \n at: {self.coordinates.ra}, {self.coordinates.dec} \n within: {self.calculate_search_arcmin()}')

            _l = (1-config['distance_tolerance'])*self.distance.value
            _u = (1+config['distance_tolerance'])*self.distance.value
            _filter = {'rgeo':f'>{_l} && <{_u}'}
            print(f'Filters: \n {_filter}')

            stars_fromDR3_dis = Vizier(columns=["*","+_r"],row_limit = -1).query_region(self.coordinates, 
                                                                                                radius=self.calculate_search_arcmin(), 
                                                                                            catalog=config['cat_star_distances'],
                                                                                            column_filters=_filter)[0]
            print(f'{len(stars_fromDR3_dis)} sources found')
            
        
        
        stars_in_region = join(stars_fromDR3,stars_fromDR3_dis,keys='Source',join_type='inner')
        dias_members=self.dias_members()
        rmRA = stars_in_region['pmRA']-dias_members['pmRA'].mean()
        rmDE = stars_in_region['pmDE']-dias_members['pmDE'].mean()

        e_rmRA = stars_in_region['e_pmRA']+np.sqrt((dias_members['e_pmRA']**2).sum())
        e_rmDE = stars_in_region['e_pmDE']+np.sqrt((dias_members['e_pmDE']**2).sum())

        µ_pec = np.sqrt(rmRA**2+rmDE**2)
        D = stars_in_region['rgeo']/1000

        stars_in_region.add_column(µ_pec, name='µ_pec', index=0)
        stars_in_region.add_column(µ_pec*D*4.74,name='v_pec',index=0)
        stars_in_region.add_column(rmRA,name='rmRA',index=5)
        stars_in_region.add_column(e_rmRA,name='e_rmRA',index=6)
        stars_in_region.add_column(rmDE,name='rmDE',index=7)
        stars_in_region.add_column(e_rmDE,name='e_rmDE',index=8)
        stars_in_region['v_pec'].unit = u.km / u.s
        stars_in_region.sort('Gmag')
        fs = stars_in_region[stars_in_region['v_pec'] >= config['v_runaway']]
        nsfs = stars_in_region[np.array(stars_in_region['v_pec'] >= config['v_walkaway']) & np.array(stars_in_region['v_pec'] < config['v_runaway'])]

        #save the tables
        stars_in_region.write(f'{self.name}/{self.name}_stars_in_region.tsv',format='ascii.ecsv',overwrite=True)
        stars_fromDR3.write(f'{self.name}/{self.name}_stars_fromDR3.tsv',format='ascii.ecsv',overwrite=True)
        stars_fromDR3_dis.write(f'{self.name}/{self.name}_stars_fromDR3_dis.tsv',format='ascii.ecsv',overwrite=True)
        fs.write(f'{self.name}/{self.name}_fs.tsv',format='ascii.ecsv',overwrite=True)
        nsfs.write(f'{self.name}/{self.name}_nsfs.tsv',format='ascii.ecsv',overwrite=True)

        if allTables:
            print(f'{len(stars_in_region)} stars in the region')
            return stars_in_region,stars_fromDR3,stars_fromDR3_dis,fs,nsfs
        else:
            print(f'{len(stars_in_region)} stars in the region')
            return stars_in_region
# def refine_stars_in_region(stars_in_region,v_dispersion=config['v_dispersion'],sigma_clip=config['sigma_clip']):
#     mask_vdisp = [value < v_dispersion for value in stars_in_region['v_pec']]

    def generate_tables(self):
        self.stars_in_region(allTables=True)

    def read_table(self, *tables):
        """
        Read one or more tables associated with the cluster.
        See the particular cluster folder for the tables present.
        Run genetate_tables() to generate the tables.

        Parameters:
        - tables (str): Variable number of table names to be read.

        Returns:
        - table (astropy.table.table.Table or list): If a single table is provided, returns the corresponding table.
        If multiple tables are provided, returns a list of tables.
        If any of the specified tables do not exist, prints a message and skips them.

        Example:
        cluster = Cluster('YourClusterName')
        member_table = cluster.read_table('members')
        stars_in_region,stars_fromDR3,stars_fromDR3_dis,fs,nsfs = cluster.read_table('stars_in_region','stars_fromDR3','stars_fromDR3_dis','fs','nsfs')
        """
        tl = []

        for table in tables:
            file_path = f'{self.name}/{self.name}_{table}.tsv'
            
            # Check if the file exists before trying to read it
            if os.path.exists(file_path):
                tl.append(Table.read(file_path, format='ascii.ecsv'))
            else:
                print(f"File {file_path} not present.")

        return tl[0] if len(tl) == 1 else tl
    def plot_theoretical_isochrone(self):
        fig = plt.figure(figsize=(8,8))
        BP_RP_theo,Gmag_theo = theoretical_isochrone(self)

        ax = fig.add_subplot()
        ax.plot(BP_RP_theo,Gmag_theo,label='Theoretical Isochrone')
        ax.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)",fontsize=14)
        ax.set_ylabel(r"$G$ (mag)",fontsize=14)
        ax.set_title(f"CMD for {self.name}",fontsize=14)
        ax.invert_yaxis()
        ax.legend()


def find_cluster(stars_in_region,refCluster,sigma = config['sigma_clip'],plx_quality = config['plx_quality']):
    """
    Takes in the stars in a region and sigma_clip their pmRA and pmDE separately to keep values within
    a specified sigma. Then further filters stars which have good parallax quality (default>10) and are 
    within the cluster radius (from dias).

    Example Usage:
    cluster = Cluster('Berkeley_97')
    my_cluster = find_cluster(cluster.stars_in_region(),refCluster=cluster.name)
    """

    within_r = Cluster(refCluster).diameter.value/2
    ms = stars_in_region[np.array(~sigma_clip(stars_in_region['pmRA'],sigma=sigma).mask)&np.array(~sigma_clip(stars_in_region['pmDE'],sigma=sigma).mask)]
    my_stars = ms[(np.array(ms['_r_1']<within_r) & np.array(ms['Plx']/ms['e_Plx']>plx_quality))]
    print('mean(v_pec): ',my_stars['v_pec'].mean())
    print('std(v_pec): ',my_stars['v_pec'].std())
    print('v_pec_max: ',my_stars['v_pec'].max())
    return my_stars



def cluster_v_dispersion(cluster_name,memb_prob=config['memb_prob'],parallax_quality_threshold=config['parallax_quality_threshold'],plot=False):
    """
    Calculate the peculiar velocities of cluster members and analyze velocity dispersion.

    Parameters:
    - cluster_name (str): Name of the cluster.
    - memb_prob (float): Membership probability threshold for selecting cluster members.
    - parallax_quality_threshold (float): Parallax quality threshold for selecting cluster members.
    - plot (bool): If True, plot a histogram with a bell curve representing the velocity dispersion.

    Returns:
    - v_pec (numpy.ndarray): Array of peculiar velocities for cluster members.
    - mean_velocity (float): Mean peculiar velocity of the cluster members.
    - std_deviation (float): Standard deviation of peculiar velocities.

    Example:
    v_pec, mean_velocity, std_deviation = cluster_v_dispersion('Berkeley_97', memb_prob=0.8, plot=True)
    """



    cluster = Cluster(cluster_name)
    cluster_members = cluster.dias_members(memb_prob,parallax_quality_threshold)
    µ_pec = np.sqrt(cluster_members['pmRA']**2+cluster_members['pmRA']**2)
    v_pec = µ_pec.value*cluster.distance.value*4.74*(u.km/u.s)/1000
    mean_velocity = np.mean(v_pec)
    std_deviation = np.std(v_pec)
    if plot:
        # Create a range of values for x-axis
        x = np.linspace(min(v_pec), max(v_pec), 100)

        # Generate the bell curve using the mean and standard deviation
        y = norm.pdf(x, mean_velocity, std_deviation)

        # Plot the histogram of velocities
        plt.hist(v_pec, bins=10, density=True, alpha=0.6, color='g')

        # Plot the bell curve
        plt.plot(x, y, 'r--', label=f'Mean: {mean_velocity:.2f}, Std Dev: {std_deviation:.2f}')

        # Add labels and legend
        plt.xlabel('Velocity')
        plt.ylabel('Probability Density')
        plt.legend()
        plt.title(f'{len(v_pec)} members')

        # Show the plot
        plt.show() 
            
    
    
    return v_pec, mean_velocity, std_deviation


def columnsplot(columnx,columny, xlabel=None, ylabel=None):
    from astropy.visualization import quantity_support

    """
    Plot two Astropy columns on the x and y axes.

    Parameters:
        columnx (astropy.table.Column): Astropy column for the x-axis.
        columny (astropy.table.Column): Astropy column for the y-axis.
        xlabel (str, optional): Custom label for the x-axis. Defaults to None.
        ylabel (str, optional): Custom label for the y-axis. Defaults to None.
    """
    with quantity_support():
        fig, ax = plt.subplots()
        ax.plot(columnx, columny, 'o')
        ax.grid(True)
        if xlabel is None:
            xlabel = f"{columnx.name} ({columnx.unit})"
        if ylabel is None:
            ylabel = f"{columny.name} ({columny.unit})"
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.show()



def theoretical_isochrone(cluster,Av=None,logage=None,metallicity=None,output=None, printing = True):
    """
    Retrieves the theoretical isochrone for a given cluster from CMD interface.
    By default takes the parameters of the cluster as input.
    Checks if data file is already present. If not, gets it.

    ## Parameters:
    -   cluster (Cluster class object): For values of extinction, age, metallicity.
    -   Av (float, optional): Extinction value. Defaults to None, which uses cluster.all['Av'].
    -   logage (float, optional): Logarithm of the age. Defaults to None, which uses cluster.all['logage'].
    -   metallicity (float, optional): Metallicity value. Defaults to None, which uses cluster.all['__Fe_H_'].
    -   output (str, optional): 'table' gets the entire table. Default is None, which returns Bp-Rp, Gmag respectively as a tuple of astropy columns.

    ## Returns:
        tuple or astropy.table.Table: If output is None, returns a tuple containing Bp-Rp and Gmag.
                                      If output is 'table', returns the entire table.

    
    ## Example Usage:
    BP_RPtheo, Gmagtheo = theoretical_isochrone(cluster, Av=3.1, logage=6, metallicity=-0.1)

    or, for a complete table,

    theo_isochrone = theoretical_isochrone(cluster, Av=3.1, logage=6, metallicity=-0.1,output='table')

    Also generates a file {cluster.name}\{cluster.name}_compare_data_out_Av3.1_age6.9_FeH-0.1.dat
    
    """
    # Accessing cluster parameters
    if Av is None:
        Av = cluster.all['Av']
    if logage is None:
        logage = cluster.all['logage']
    if metallicity is None:
        metallicity = cluster.all['__Fe_H_']



    #Check if parameters are unchanged
    if (str(Av) == str(cluster.all['Av'])) and (str(logage) == str(cluster.all['logage'])):
        compate_data_out_path = os.path.join(cluster.name,(cluster.name+'_compare_data_out.dat'))
        if printing:
            print(f'{cluster.name}')
            print(f'Av: {str(Av)}')
            print(f'logage: {str(logage)}')
            print(f'[Fe/H]: {str(metallicity)}')
    #otherwise, make a new file to save the new isochrone
    else:
        compate_data_out_path = os.path.join(cluster.name,(cluster.name+f'_compare_data_out_Av{str(Av)}_age{str(logage)}_FeH{str(metallicity)}.dat'))
        if printing:
            print(f'{cluster.name} (Changed)')
            print(f"Av: {str(cluster.all['Av'])} --> {str(Av)}")
            print(f"logage: {str(cluster.all['logage'])} --> {str(logage)}")
            print(f"[Fe/H]: {str(cluster.all['__Fe_H_'])} --> {str(metallicity)}")

    if not os.path.exists(compate_data_out_path):
        if printing:
            print('Getting theoretical isochrone')
        s = Service()
        options = webdriver.ChromeOptions()
        browser = webdriver.Chrome(service=s, options=options)
        browser.get('http://stev.oapd.inaf.it/cgi-bin/cmd')

        #Evolutionary Tracks #from config
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[3]/td[1]/input[1]").click() #PARSEC version 2.0
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[5]/td/input").click() #+ COLIBRI S_37
        #Phtotometric System #from config
        photometricSystem = Select(browser.find_element(By.XPATH,"//select[@name='photsys_file']")) #dropdown list for available photometric systems
        photometricSystem.select_by_value("YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat") # Gaia EDR3 bands
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[2]/table/tbody/tr[5]/td[1]/input").click() #for PARSEC 2.0 #As above, but adopting revised SED for Vega from Bohlin et al. (2020) (namely CALSPEC alpha_lyr_stis_010.fits).
        #Phtotometric System #from config
        #Circumstellar Dust 
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[3]/font/table/tbody/tr[3]/td[1]/input").click() #for M stars: No dust
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[3]/font/table/tbody/tr[3]/td[2]/input").click() #for C stars: No dust
        #Interstellar extinction 
        Av_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[4]/input")
        Av_field.clear()
        Av_field.send_keys(str(Av))
        #Long Period Variability #from config
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[5]/table/tbody/tr[4]/td[1]/input").click() #3. Periods from Trabucchi et al. (2021).
        #Initial Mass Function #from config
        InitialMassFunction = Select(browser.find_element(By.XPATH,"//select[@name='imf_file']"))
        InitialMassFunction.select_by_value("tab_imf/imf_kroupa_orig.dat")
        #ages
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[4]/td[1]/input").click() #click log(age/yr)
        initialAge_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[4]/td[2]/input")
        initialAge_field.clear()
        initialAge_field.send_keys(str(logage))
        finalAge_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[4]/td[3]/input")
        finalAge_field.clear()
        finalAge_field.send_keys(str(logage))
        step_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[4]/td[4]/input")
        step_field.clear()
        step_field.send_keys(0)
        #metallicities
        selectMbyH = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[7]/td[1]/input").click() #click [M/H]
        initialMetallicity_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[7]/td[2]/input")
        initialMetallicity_field.clear()
        initialMetallicity_field.send_keys(str(metallicity))
        finalMetallicity_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[7]/td[3]/input")
        finalMetallicity_field.clear()
        finalMetallicity_field.send_keys(str(metallicity))
        #Submit #dosen't change
        browser.find_element(By.XPATH,"/html/body/form/div/input[4]").click()
        browser.find_element(By.XPATH,"/html/body/form/fieldset[1]/p[1]/a").click()
        data_out = browser.find_element(By.XPATH,'/html/body/pre').text
        if printing:
            print('Obtained isochrone from CMD')
        #save the table
        with open(compate_data_out_path, 'w') as f:
            f.write(data_out)
        if printing:
            print(f"New isochrone file: {compate_data_out_path}")
        #close the browser
        browser.close()
    else:
        if printing:
            print('Theoretical isochrone exists, reading it.')
    # Read the DAT file into an Astropy table
    theoretical_data = Table.read(compate_data_out_path, format='ascii')
    # Define the new column names
    new_column_names = ['Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe', 'logg', 'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2', 'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc', 'Xn', 'Xo', 'Cexcess', 'Z', 'Teff0', 'omega', 'angvel', 'vtaneq', 'angmom', 'Rpol', 'Req', 'mbolmag', 'G_fSBmag', 'G_BP_fSBmag', 'G_RP_fSBmag', 'G_fSB', 'G_f0', 'G_fk', 'G_i00', 'G_i05', 'G_i10', 'G_i15', 'G_i20', 'G_i25', 'G_i30', 'G_i35', 'G_i40', 'G_i45', 'G_i50', 'G_i55', 'G_i60', 'G_i65', 'G_i70', 'G_i75', 'G_i80', 'G_i85', 'G_i90', 'G_BP_fSB', 'G_BP_f0', 'G_BP_fk', 'G_BP_i00', 'G_BP_i05', 'G_BP_i10', 'G_BP_i15', 'G_BP_i20', 'G_BP_i25', 'G_BP_i30', 'G_BP_i35', 'G_BP_i40', 'G_BP_i45', 'G_BP_i50', 'G_BP_i55', 'G_BP_i60', 'G_BP_i65', 'G_BP_i70', 'G_BP_i75', 'G_BP_i80', 'G_BP_i85', 'G_BP_i90', 'G_RP_fSB', 'G_RP_f0', 'G_RP_fk', 'G_RP_i00', 'G_RP_i05', 'G_RP_i10', 'G_RP_i15', 'G_RP_i20', 'G_RP_i25', 'G_RP_i30', 'G_RP_i35', 'G_RP_i40', 'G_RP_i45', 'G_RP_i50', 'G_RP_i55', 'G_RP_i60', 'G_RP_i65', 'G_RP_i70', 'G_RP_i75', 'G_RP_i80', 'G_RP_i85', 'G_RP_i90']
    # Rename the columns
    for old_name, new_name in zip(theoretical_data.colnames, new_column_names):
        theoretical_data.rename_column(old_name, new_name)
    # Calculate BP-RP and add column at the end
    theoretical_data['BP-RP'] = theoretical_data['G_BP_fSBmag']-theoretical_data['G_RP_fSBmag']
    theoretical_data['Gmag'] = theoretical_data['G_fSBmag'] + 5*np.log10(cluster.distance.value)-5
    BP_RPtheo = theoretical_data['BP-RP'] #BP-RP
    Gmagtheo = theoretical_data['Gmag']#-cluster['Av']
    if output == None:
        return BP_RPtheo,Gmagtheo
    if output == 'table':
        return theoretical_data


def get_runaways(cluster,fs,theoretical_data):
    """
    Selects stars from the input astropy table fs which are runaway star candidates by tracing them back to the cluster.
    If the stars happens to be in the cluster in the past 100kyr then it is selected to be a candidate.
    Then the temperature of the star is estimated by comparing its color (Bp-Rp magnitude) with the theoretical isochrone that is provided.

    Parameters:
    -    cluster (Cluster class object): The cluster compared whose runaways are to be found.
    -    fs or nsfs (astropy table): A table of the fast stars (relative to the cluster mean proper motion) in the region (`search_arcmin`) 
        around the cluster centre. Can also take `nsfs` and trace them back with the same algorithm.
    -   theoretical_data (astropy table ascii format downloaded from CMD): theoretical isochrone to estimate the temperatures of the runaways (or walkaways)

    Returns:
    - run (astropy table): A table of the selected stars out of fs which can be traced back within 100kyr or less 
        with their temperature estimates added in a column.

    Example Usage:
    cluster = Cluster('Berkeley_97')
    theoretical_data = theoretical_isochrone(cluster,output="table",printing=False,logage=7)
    fs = cluster.read_table('fs')
    runaways = get_runaways(cluster,fs,theoretical_data)
    display(runaways)
    
    """


    
    if fs['v_pec'].max() < 17.6:
        name = "walkaways"
        print("Stars given in the fs variable have their maximum v_pec <17.6km/s. These are walkaways.")
        
    else:
        name = "runaways"

    file_path = f'{cluster.name}/{cluster.name}_{name}.tsv'     #`name` is either "walkaways" or "runaways"
    
    if os.path.exists(file_path):
        # If the file exists, read it and return its contents
        print(f'{file_path} exists. Skipping traceback and temp estimates.')
        run = Table.read(file_path, format='ascii.ecsv')
        print(f'{len(fs)} stars had been checked, {len(run)} traced back to cluster')
        return run
    
    else:
        #do the tracebacks
        #def the proper motin of cluster
        times = np.linspace(-1e5, 0, 100)*u.year #time in years
        cluster_motion = cluster.coordinates.apply_space_motion(dt=times)

        #read the fast stars and make new table
        fast_stars = QTable()
        fast_stars['coordinates'] = SkyCoord(ra=fs['RA_ICRS_1'],dec=fs['DE_ICRS_1'],pm_ra_cosdec=fs['pmRA'],pm_dec=fs['pmDE'],obstime = (Time('J2000')+1*u.Myr))
        fast_stars['Source'] = fs['Source']
        print(f'Tracing back {len(fast_stars)} stars...')


        # Preallocate array to hold separations for all stars
        sep_array = np.zeros((len(fast_stars), len(times)))

        for i, star in enumerate(fast_stars):
            # Apply space motion for the star for all time steps
            star_motion = fast_stars['coordinates'][i].apply_space_motion(dt=times)

            # Calculate separation for all time steps at once
            sep = cluster_motion.separation(star_motion).to(u.arcmin)
            sep_array[i, :] = sep

        # Convert sep_array to list of lists
        sep_list = sep_array.tolist()

        # Create dictionary with star names as keys and separations as values
        sep_dict = {star['Source']: sep_list[i] for i, star in enumerate(fast_stars)}

        # make a list of stars that were inside the cluster 
        selected_stars_dict = {}
        threshold_separation = cluster.diameter.value/2
        for source_id, separations in sep_dict.items():
            # Check if any separation is less than the threshold
            if any(sep < threshold_separation for sep in np.array(separations)):
                # Add the source id to the list
                selected_stars_dict[source_id] = sep_dict[source_id]

        # Make list with temp
        selected_stars_dict_with_temp = {}

        for source_id, separations in selected_stars_dict.items():
            new_star_bp_rp = fs[fs['Source']==source_id]['BP-RP']

            # Get the BP-RP and Temperature columns from the table
            temperature_column = theoretical_data['logTe']
            bp_rp_column = theoretical_data['BP-RP']
            # Calculate the differences between the new star's BP-RP value and all values in the table
            differences = np.abs(bp_rp_column - new_star_bp_rp)
            # Find the index of the star with the closest BP-RP value
            closest_star_index = np.argmin(differences)
            # Get the temperature of the closest star
            closest_star_temperature = temperature_column[closest_star_index]
            closest_star_temperature = 10 ** closest_star_temperature

            selected_stars_dict_with_temp[source_id] = (selected_stars_dict[source_id],closest_star_temperature)


        sources_to_mask = list(selected_stars_dict_with_temp.keys())
        mask = [source in sources_to_mask for source in fs['Source']]
        run = fs[mask]
        run.add_column(np.array(list(selected_stars_dict_with_temp.values()),dtype=object)[:,1],name='Temp. Est',index =0)
        run['Temp. Est'] = run['Temp. Est']*u.K
        run.sort('Temp. Est',reverse=True)
        # Save
        run.write(file_path,format='ascii.ecsv',overwrite=True)
        run.to_pandas().to_excel(os.path.join(cluster.name,f'{cluster.name}_{name}.xlsx'), index=False)

        return run

#plot the find_members and dias_members on the search region fits
#find out the reason why some clusters were not searched. probably because of the search distance
#3D visualisation
#unique config files for each cluster. add string on top of it indicating last change date