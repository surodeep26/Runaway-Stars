import os
from astropy.table import Table, Column, QTable, join
import yaml
import pandas  # Renamed for clarity
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from IPython.display import IFrame
from matplotlib.widgets import CheckButtons

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
# from runaway_functions import *
# import pandas as pd
import shutil
from regions import CircleSkyRegion


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
    def __init__(self, cluster_name,version='dr3'):
        # Assuming 'cluster_list' is a predefined list containing your Astropy table
        self.cluster_table = cluster_list[cluster_list['Cluster'] == cluster_name]
        self.all = self.cluster_table[0]
        self.name = self.cluster_table['Cluster'][0]


        if len(self.cluster_table) == 0:
            raise ValueError(f"No data found for cluster '{cluster_name}' in the cluster list.")
        if version == 'dr3':
            d = self.dias_members()
            d.sort('Gmag')
            dm = d[d['Gmag']<config['gmag_lim']]['Source','RAdeg','DEdeg','Gmag','e_Gmag','Plx','e_Plx','Pmemb']
            # print(len(dm))
            sources = np.array(dm['Source'])
            sources_int = sources.astype(np.int64)
            sr=self.stars_in_region(no_rmRArmDE=True)
            sr = sr[sr['RUWE']<config['ruwe_lim']]
            mask = np.isin(sr['Source'], sources_int)
            dr3dias = sr[mask]
            # print(dr3dias['rgeo'].mean(),dr3dias['rgeo'].std())
            self.all['Dist'] = dr3dias['rgeo'].mean()
            self.all['e_Dist'] = dr3dias['rgeo'].std()
            self.all['pmRA'] = dr3dias['pmRA'].mean()
            self.all['e_pmRA'] = dr3dias['pmRA'].std()
            self.all['pmDE'] = dr3dias['pmDE'].mean()
            self.all['e_pmDE'] = dr3dias['pmDE'].std()
            self.all['Dist'] = dr3dias['rgeo'].mean()
            self.all['Plx'] = dr3dias['Plx'].mean()
            self.all['e_Plx'] = dr3dias['Plx'].std()
        elif version == 'dr2':
            dr3dias = self.dias_members()
        # Extracting values from the table and setting attributes
        if not os.path.exists(f'{self.name}'):
            os.mkdir(f'{self.name}')
        ti = (Time('J2000')+1*u.Myr)
        self.coordinates = SkyCoord(ra=self.cluster_table['RA_ICRS'],dec=self.cluster_table['DE_ICRS'],pm_ra_cosdec=self.cluster_table['pmRA'],pm_dec=self.cluster_table['pmDE'],obstime=ti)
        self.diameter = self.cluster_table['Diameter'][0]*u.arcmin
        self.N = self.cluster_table['N'][0] #no. of members
        self.diameter_dias = (self.cluster_table['r50'][0]*u.deg*2).to(u.arcmin)
        self.distance = self.cluster_table['Dist'][0]*u.pc
        self.RV = self.cluster_table['RV'][0]*u.km/u.s
        self.Av = self.cluster_table['Av'][0]
        self.logage = self.cluster_table['logage'][0]
        self.age = 10**(self.cluster_table['logage'][0])
        self.FeH = self.cluster_table['__Fe_H_'][0]
        self.members = dr3dias #gets dias_members with the default settings for memb_prob and parallax_threshold_quality
        #to change the memb_prob and parallax_threshold_quality filters for member selection, use something like dias_members(0.6,11)

    
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

    def plot_search_region(self, extra=config['Cluster']['search_extent'],display=True,**kwargs):
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
                                        radius=2*search_arcmin,
                                        **kwargs)
            # Extract the WCS information
            wcs = WCS(images[0][0].header)
            hdu = fits.PrimaryHDU(data=images[0][0].data, header=images[0][0].header)
            hdulist = fits.HDUList([hdu])
            
            # Save the fits file
            hdulist.writeto(fits_file_path, overwrite=True)
        if not display:
            return None
        # Plot the image
        print(f'Plotting a {extra} pc region around cluster.')
        # print(wcs.wcs.ctype) #projection is TAN by default for skyview images

            
        fig, ax = plt.subplots(figsize=(6,6), subplot_kw={'projection': wcs})
        ax.imshow(images[0][0].data, cmap='gray')
        data = images[0][0].data
        # Add labels and title
        ax.set_xlabel('Right Ascension (degrees)')
        ax.set_ylabel('Declination (degrees)')
        ax.set_title(f"{self.name} with {extra}pc search region")
        ax.grid(color='lightgrey',ls='dotted')

        # Add a circle representing the cluster radius
        circle_cluster = plt.Circle((data.shape[1] / 2, data.shape[0] / 2), radius=(data.shape[0] / 2)*self.diameter/(2*search_arcmin),
                        edgecolor='red', facecolor='none',ls='dashed',label=f'Cluster Diameter = {self.diameter}')
        circle_search_region = plt.Circle((data.shape[1] / 2, data.shape[0] / 2), radius=(data.shape[0] / 2)*search_arcmin/search_arcmin,
                    edgecolor='green', facecolor='none',ls='dashed',label=f'Search Region Diameter = {2*search_arcmin}')
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
        images[0].close()  # Close the opened FITS file

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
    

    def stars_in_region(self,allTables=False,no_rmRArmDE=False):
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
            # print(f'{stars_in_region_path} exists in {self.name} with {len(stars_in_region)} stars')

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

        # relative to the cluster stuff
        if not no_rmRArmDE:
            dias_members= self.members
            rmRA = stars_in_region['pmRA']-dias_members['pmRA'].mean()
            rmDE = stars_in_region['pmDE']-dias_members['pmDE'].mean()

            e_rmRA = stars_in_region['e_pmRA']+np.sqrt((dias_members['e_pmRA']**2).mean())
            e_rmDE = stars_in_region['e_pmDE']+np.sqrt((dias_members['e_pmDE']**2).mean())

            e_rmRA = stars_in_region['e_pmRA']+self.all['e_pmRA']
            e_rmDE = stars_in_region['e_pmDE']+self.all['e_pmDE']

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

            stars_in_region = stars_in_region['RA_ICRS_1','DE_ICRS_1','e_RA_ICRS','e_DE_ICRS','_r_1',
                                            'HIP','TYC2','Source','rgeo','Plx','e_Plx',
                                            'v_pec','µ_pec','rmRA','e_rmRA','rmDE','e_rmDE','pmRA','pmDE','e_pmRA','e_pmDE',
                                            'RUWE','Teff','logg','Gmag','BP-RP','BPmag','RPmag','RV','e_RV',
                                            'b_rgeo','B_rgeo','RAVE5','RAVE6']

            fs = stars_in_region[stars_in_region['v_pec'] >= config['v_runaway']]
            nsfs = stars_in_region[np.array(stars_in_region['v_pec'] >= config['v_walkaway']) & np.array(stars_in_region['v_pec'] < config['v_runaway'])]

            #save the tables
            fs.write(f'{self.name}/{self.name}_fs.tsv',format='ascii.ecsv',overwrite=True)
            nsfs.write(f'{self.name}/{self.name}_nsfs.tsv',format='ascii.ecsv',overwrite=True)

        stars_in_region.write(f'{self.name}/{self.name}_stars_in_region.tsv',format='ascii.ecsv',overwrite=True)
        stars_fromDR3.write(f'{self.name}/{self.name}_stars_fromDR3.tsv',format='ascii.ecsv',overwrite=True)
        stars_fromDR3_dis.write(f'{self.name}/{self.name}_stars_fromDR3_dis.tsv',format='ascii.ecsv',overwrite=True)


        if allTables and not no_rmRArmDE:
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
    print(f"{len(my_stars)} kinematic members")
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
    cluster_members = cluster.members
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

    file_path = f'{cluster.name}/{cluster.name}_{name}_all.tsv'     #`name` is either "walkaways" or "runaways"
    
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
        # threshold_separation = cluster.diameter.value #deleteme
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
        run.add_column(np.array(list(selected_stars_dict_with_temp.values()),dtype=object)[:,1],name='Temp. Est',index =1)
        run['Temp. Est'] = run['Temp. Est']*u.K
        run.sort('Temp. Est',reverse=True)
        # Save
        run.write(file_path,format='ascii.ecsv',overwrite=True)
        run.to_pandas().to_excel(os.path.join(cluster.name,f'{cluster.name}_{name}_all.xlsx'), index=False)

        mask = [T > config['runaway_temp'] for T in run['Temp. Est']]
        runaways = run[mask]
        runaways.write(file_path.replace("_all",""),format='ascii.ecsv',overwrite=True)
        runaways.to_pandas().to_excel(os.path.join(cluster.name,f'{cluster.name}_{name}.xlsx'), index=False)


        return run

def get_coord(runaway):
    return coord.SkyCoord(ra=runaway['RA_ICRS_1']*u.deg,dec=runaway['DE_ICRS_1']*u.deg, pm_ra_cosdec=runaway['rmRA']*u.mas/u.year,pm_dec=runaway['rmDE']*u.mas/u.year, frame='icrs')
def plot_cmd(cluster,save=False):
    # Plot CMD
    BP_RP_theo, Gmag_theo = theoretical_isochrone(cluster)

    # fig = plt.figure(figsize=(19, 17))
    # # fig = plt.figure()
    # gs = GridSpec(1, 1, width_ratios=[1, 1])

    # # Line Theorietical curve
    # # ax1 = fig.add_subplot(gs[0, 0])
    # ax1 = fig.subplots()

    fig,ax1 = plt.subplots(figsize=(10, 8.8))

    # plt.rcParams['figure.figsize'] = (10, 8.8)




    ax1.set_ylim(bottom= min(Gmag_theo)-4,top=18)
    ax1.set_xlim(left=min(BP_RP_theo)-0.5, right=max(BP_RP_theo))

    ax1.plot(BP_RP_theo, Gmag_theo, label='Theoretical Isochrone')
    ax1.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)", fontsize=14)
    ax1.set_ylabel(r"$G$ (mag)", fontsize=14)
    ax1.set_title(f"CMD for {cluster.name}", fontsize=14)
    ax1.invert_yaxis()

    #Scatter stars in the region
    sir = cluster.stars_in_region()
    sir_gmag,sir_bp_rp = sir['Gmag'],sir['BP-RP']
    ax1.scatter(sir_bp_rp,sir_gmag,s=2, color='grey',label=f'{len(sir)} stars in the region {cluster.calculate_search_arcmin()}')

    # Scatter dias members
    dias_members = cluster.members
    dias_gmag,dias_bp_rp = dias_members['Gmag'],dias_members['BP-RP']
    ax1.scatter(dias_bp_rp,dias_gmag,s=15, color='black',label=f'{len(dias_members)} Dias Members P > {config["memb_prob"]}')

    # Scatter my members
    my_members = find_cluster(sir,refCluster=cluster.name)
    my_gmag,my_bp_rp = my_members['Gmag'],my_members['BP-RP']
    ax1.scatter(my_bp_rp,my_gmag,s=15,alpha =0.7,marker='x', color='blue',label=rf'{len(my_members)} My Members $\sigma$ < {np.std(my_members["v_pec"]):.2f}')


    # Table for cluster parameters
    cluster_table = [
        ['Members',len(dias_members)],
        [r'Metallicity $[Fe/H]$',cluster.FeH],
        ['Logage',cluster.logage],
        ['Av',cluster.Av],
        ['Distance (pc)',str(round(cluster.distance.value))+"$\pm$"+f'{cluster.all["e_Dist"]}'+'pc']
    ]
    table_bbox = [0.0, 0.84, 0.44, 0.16]  # [left, bottom, width, height]

    table = ax1.table(cellText=cluster_table, cellLoc='right', loc='upper left',bbox=table_bbox)
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for key, cell in table._cells.items():
        cell.set_linewidth(0.5)  # Set the border width
        cell.set_edgecolor('lightgray')  # Set the border color




    # Scatter runaways
    runaways_all = cluster.read_table('runaways_all')
    _runaways_gmag,_runaways_bp_rp = runaways_all['Gmag'],runaways_all['BP-RP']
    mask = [T > config['runaway_temp'] for T in runaways_all['Temp. Est']]
    runaways = runaways_all[mask]
    runaways_gmag,runaways_bp_rp = runaways['Gmag'],runaways['BP-RP']
    ## all runaways
    ax1.scatter(_runaways_bp_rp,_runaways_gmag,s=8,alpha=0.5, 
                c=runaways_all['Temp. Est'],
                cmap='spring_r',norm=plt.Normalize(4000, 23000),
                label=rf'{len(runaways_all)} Runaway(s) $T_{{max}}$ = {max(runaways_all["Temp. Est"]):,.0f} K')

    ## main runaways T > 10,000K
    global annotations
    annotations = []

    def on_release(event):
        global annotations, table2
        for ann in annotations:
            ann.remove()
        annotations = []

        try:
            if table2:
                table2.remove()
        except:
            pass
        plt.draw() # for the release to remove functionality

    def onpick3(event):
        global annotations, table2
        for ann in annotations:
            ann.remove()
        annotations = []
        
        ind = event.ind
        for index in ind:
            index = ind[0]
            temp_est = runaways['Temp. Est'][index]
            ann = ax1.annotate(f'{temp_est:,.0f} K', (runaways_bp_rp[index], runaways_gmag[index]), xytext=(0, -25), textcoords='offset points', fontsize=12, color='black', ha='center', fontweight='bold')
            ann.set_path_effects([pe.withStroke(linewidth=4, foreground='white')])
            annotations.append(ann)
            
            # Add the table below the annotation
            table_data = [
                ['SourceID', f'{runaways[index]["Source"]}'],
                ['Temp. Est.', f'{runaways[index]["Temp. Est"]:,.0f} K'],
                ['Coord.', f'{runaways["RA_ICRS_1"][index]:.4f} {"+" if runaways["DE_ICRS_1"][index] >= 0 else ""}{runaways["DE_ICRS_1"][index]:.4f}'],
                ['Gmag', f'{runaways["Gmag"][index]}'],
                ['Dist.', f'{runaways["rgeo"][index]:.0f}'+'$\pm$'+f'{(runaways["rgeo"][index]-runaways["b_rgeo"][index]):.0f}']

                # Add more rows with actual values
            ]
            table_bbox = [0.55, 0.0, 0.45, 0.2]  # [left, bottom, width, height]
            table2 = ax1.table(cellText=table_data, cellLoc='right', loc='lower center', bbox=table_bbox)
            table2.auto_set_font_size(False)
            table2.auto_set_column_width(col=[0,1])
            table2.set_fontsize(8)
            for key, cell in table2._cells.items():
                cell.set_linewidth(0.5)  # Set the border width
                cell.set_edgecolor('lightgray')  # Set the border color
            
        plt.draw()


    scatter_main = ax1.scatter(runaways_bp_rp,runaways_gmag,picker=True,s=30, 
                            c=runaways['Temp. Est'],cmap='spring_r',norm=plt.Normalize(4000, 23000),
                            label=rf'{len(runaways)} Runaway(s) with T > {config["runaway_temp"]:,} K')
    # Connect the pick_event to the onpick3 function
    fig.canvas.mpl_connect('pick_event', onpick3)
    fig.canvas.mpl_connect('button_release_event', on_release)

    colorbar = fig.colorbar(scatter_main,ax=ax1)
    colorbar.set_label('Temperature (K)')
    ax1.legend(loc='upper right')
    if save:
        plt.savefig(f'{cluster.name}/{cluster.name}_cmd.{save}')



def runCode(cluster,save='png',psr=False):
    if isinstance(cluster,str):
        cluster = Cluster(cluster)

    elif isinstance(cluster,Cluster):
        cluster=cluster
    print(f'{cluster.name:=>50}'+f'{"":=<50}')
    cluster = Cluster(cluster.name)
    cluster.generate_tables()
    theoretical_data = theoretical_isochrone(cluster,output="table",printing=False)
    print(f'{cluster.name+": traceback":->50}'+f'{"":-<50}')
    fs = cluster.read_table('fs')
    runaways_all = get_runaways(cluster,fs,theoretical_data)
    # display(runaways_all)
    mask = [T > config['runaway_temp'] for T in runaways_all['Temp. Est']]
    runaways = runaways_all[mask]
    print(f'{cluster.name+": runaways and isochrone plot":->50}'+f'{"":-<50}')
    plot_cmd(cluster,save=save)
    if psr:
        plot_traceback(cluster,save=save)
    plot_traceback_clean(cluster,save=save)
    display(f"{len(runaways)} Runaway(s) found",runaways)


def plot_traceback_clean(cluster,save=False):
    search_arcmin = cluster.calculate_search_arcmin()
    cluster.plot_search_region(display=None,pixels='800')
    fits_path = f'{cluster.name}/{cluster.name}_extra10pc.fits'
    fits_file = fits.open(fits_path)
    image = fits_file[0]
    wcs = WCS(image.header)

    fig, ax2 = plt.subplots(subplot_kw={'projection': wcs}, figsize=(10, 8.8))
    ax2.imshow(image.data, cmap='gray')

    # Add a circle representing the cluster radius
    circle_cluster = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), radius=(image.shape[0]/2)*cluster.diameter/(2*search_arcmin),
                    edgecolor='red', facecolor='none', ls='dashed', label=f'Cluster Diameter = {cluster.diameter}',linewidth=1)
    circle_search_region = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), radius=(image.shape[0]/2)*search_arcmin/search_arcmin,
                edgecolor='green', facecolor='none', ls='dashed', label=f'Search Region Diameter = {2*search_arcmin}',linewidth=1.5)
    ax2.add_artist(circle_cluster)
    ax2.add_artist(circle_search_region)
    # Read runaways
    runaways_all = cluster.read_table('runaways') #changing this to runaways_all plots all the runaways, not just the >10000K ones
    runaways_all['Source'] = runaways_all['Source'].astype(str)
    runaways_all['Gmag'] = (runaways_all['Gmag'].round(2)).astype(str)
    runaways_all_coords = SkyCoord(ra=runaways_all['RA_ICRS_1'], dec=runaways_all['DE_ICRS_1'])
    runaways_all_pix_coords = wcs.world_to_pixel_values(runaways_all_coords.ra, runaways_all_coords.dec)
    _ = search_arcmin.round(-1)/2 #round to the nearest 5 (by round first to nearest 10 and then divide by 2
    scalebar_length = ((_.to(u.rad))*(cluster.distance.to(u.m)).to(u.pc)).round(2)
    __ = cluster.distance.value/1000
    add_scalebar(ax2, _, color="yellow", label=f"or {scalebar_length.value}pc (at dist {__:.2f}kpc)", size_vertical=0.5)
    x_min, x_max = ax2.get_xlim()
    y_min, y_max = ax2.get_ylim()
    ax2.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontsize=12,fontweight='bold')

    ax2.set_xlabel('Right Ascension (hms)')
    ax2.set_ylabel('Declination (degrees)')
    nam = f"{(cluster.name).replace('_',' ')}"
    ax2.set_title(f"{nam} with a {config['Cluster']['search_extent']} pc search region")
    ax2.grid(color='lightgrey', ls='dotted')
    #######################################
    if len(runaways_all)>0:
        runaway_00, runaway_apdp,runaway_apdm,runaway_amdp,runaway_amdm = [coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=runaways_all['rmRA'],pm_dec=runaways_all['rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']+runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']+runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']+runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']-runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']-runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']+runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']-runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']-runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr))]

        earlier_runaway_00 = runaway_00.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_apdp = runaway_apdp.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_apdm = runaway_apdm.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_amdp = runaway_amdp.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_amdm = runaway_amdm.apply_space_motion(dt = -100_000*u.year)

        earlier_runaway_00_px = wcs.world_to_pixel_values(earlier_runaway_00.ra, earlier_runaway_00.dec)
        earlier_runaway_apdp_px = wcs.world_to_pixel_values(earlier_runaway_apdp.ra, earlier_runaway_apdp.dec)
        earlier_runaway_apdm_px = wcs.world_to_pixel_values(earlier_runaway_apdm.ra, earlier_runaway_apdm.dec)
        earlier_runaway_amdp_px = wcs.world_to_pixel_values(earlier_runaway_amdp.ra, earlier_runaway_amdp.dec)
        earlier_runaway_amdm_px = wcs.world_to_pixel_values(earlier_runaway_amdm.ra, earlier_runaway_amdm.dec)
        runaway_00_px = wcs.world_to_pixel_values(runaway_00.ra, runaway_00.dec)
        # Calculate the displacement vectors
        delta_x = earlier_runaway_apdp_px[0]-earlier_runaway_amdp_px[0]
        delta_y = earlier_runaway_apdp_px[1]-earlier_runaway_apdm_px[1]

        for pxx, pxy, dx, dy in zip(earlier_runaway_00_px[0], earlier_runaway_00_px[1], delta_x, delta_y):
            ellipse = patches.Ellipse(
                (pxx, pxy),
                width=1.5*dx,
                height=1.5*dy,
                fill=True,
                color='g',
                alpha=0.2
            )
            ax2.add_patch(ellipse)

        for c in [earlier_runaway_apdp_px,earlier_runaway_apdm_px,earlier_runaway_amdp_px,earlier_runaway_amdm_px]:
            delta_x = c[0] - runaway_00_px[0]
            delta_y = c[1] - runaway_00_px[1]
            # Draw the vectors
            ax2.quiver(runaway_00_px[0], runaway_00_px[1], delta_x, delta_y, angles='xy', scale_units='xy', scale=1, color='limegreen', width=0.001)
    ############################
    else:
        print("No runaways found")
    # Scatter runaways
    scatter_runaways_all = ax2.scatter(runaways_all_pix_coords[0], runaways_all_pix_coords[1], s=30,picker=True,pickradius=20,
                                    c=runaways_all['Temp. Est'], cmap='spring_r', norm=plt.Normalize(4000, 23000), 
                                    label=f'{len(runaways_all)} Runaway(s)')
    legend = ax2.legend(loc='upper right')

    global annotations
    annotations = []

    def on_release(event):
        global annotations, table2
        try:
            if table2:
                table2.remove()
        except:
            pass
        fig.canvas.draw() # for the release to remove functionality

    def onpick3(event):
        global annotations, table2
        
        ind = event.ind
        
        # print(ind,type(ind),len(ind),event.artist)
        # print(result)
        import matplotlib
        if isinstance(event.artist, matplotlib.collections.PathCollection) and len(ind) == 1:
            for index in ind:
                index = ind[0]
                temp_est = runaways_all['Temp. Est'][index]
                # ann = ax2.annotate(f'{temp_est:,.0f} K', (0,0), xytext=(0, 0), textcoords='offset points', fontsize=12, color='black', ha='center', fontweight='bold')
                # ann.set_path_effects([pe.withStroke(linewidth=4, foreground='white')])
                # annotations.append(ann)
                # Add the table below the annotation
                if (runaways_all[index]["Source"][:3] == 'PSR'): #then this is a psr, table has different format
                    table_data = [
                        ['JNAME', f'{runaways_all[index]["Source"]}'],
                        ['Separation', f'{runaways_all[index]["RAVE6"]} arcmin'],
                        ['Coord.', f'{runaways_all["RA_ICRS_1"][index]:.4f} {"+" if runaways_all["DE_ICRS_1"][index] >= 0 else ""}{runaways_all["DE_ICRS_1"][index]:.4f}'],
                        ['Age', f'{runaways_all["Gmag"][index]}'], #this contains the age
                        ['Dist.', f'{runaways_all["rgeo"][index]:.0f}'+' pc'], #this contains the distance
                        ['Dist_DM', f'{runaways_all["RAVE5"][index]}'] #this contains the dist_dm #dispersion measured distance
                        # Add more rows with actual values
                    ]
                else: #it is a runaway star
                    # Add the table below the annotation
                    table_data = [
                        ['SourceID', f'{runaways_all[index]["Source"]}'],
                        ['Temp. Est.', f'{runaways_all[index]["Temp. Est"]:,.0f} K'],
                        ['Coord.', f'{runaways_all["RA_ICRS_1"][index]:.4f} {"+" if runaways_all["DE_ICRS_1"][index] >= 0 else ""}{runaways_all["DE_ICRS_1"][index]:.4f}'],
                        ['Gmag', f'{runaways_all["Gmag"][index]}'+' mag'],
                        ['Dist.', f'{runaways_all["rgeo"][index]:.0f}'+'$\pm$'+f'{(runaways_all["rgeo"][index]-runaways_all["b_rgeo"][index]):.0f}'+' pc']
                        # Add more rows with actual values
                    ]

                table_bbox = [0.0, 0.0, 0.45, 0.2]  # [left, bottom, width, height]
                table2 = ax2.table(cellText=table_data, cellLoc='right', loc='lower center', bbox=table_bbox)
                table2.auto_set_column_width(col=[0,1])
                for key, cell in table2._cells.items():
                    cell.set_linewidth(0.5)  # Set the border width
                    cell.set_edgecolor('lightgray')  # Set the border color
                    # cell.get_text().set_weight('bold')
            
            fig.canvas.draw()
    # Connect the pick_event to the onpick3 function
    fig.canvas.mpl_connect('pick_event', onpick3)
    fig.canvas.mpl_connect('button_release_event', on_release)



    # colorbar = fig.colorbar(scatter_runaways_all,ax=ax2,label='Temperature (K)')
    plt.tight_layout()

    if save:
        plt.savefig(f'{cluster.name}/{cluster.name}_traceback.{save}')


def plot_traceback(cluster,save=False,psr_circles=True):
    search_arcmin = cluster.calculate_search_arcmin()
    cluster.plot_search_region(display=None,pixels='800')
    fits_path = f'{cluster.name}/{cluster.name}_extra10pc.fits'

    fits_file = fits.open(fits_path)
    image = fits_file[0]

    wcs = WCS(image.header)
    stars_in_region = cluster.stars_in_region()
    fig, ax2 = plt.subplots(subplot_kw={'projection': wcs}, figsize=(10, 8))
    ax2.imshow(image.data, cmap='gray')
    stars_in_region_coords = SkyCoord(ra=stars_in_region['RA_ICRS_1'], dec=stars_in_region['DE_ICRS_1'])
    stars_in_region_pix_coords = wcs.world_to_pixel_values(stars_in_region_coords.ra, stars_in_region_coords.dec)
    members = cluster.members
    dias_members_coords = SkyCoord(ra=members['RA_ICRS_1'], dec=members['DE_ICRS_1'])
    dias_members_pix_coords = wcs.world_to_pixel_values(dias_members_coords.ra, dias_members_coords.dec)
    # Add a circle representing the cluster radius
    circle_cluster = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), radius=(image.shape[0]/2)*cluster.diameter/(2*search_arcmin),
                    edgecolor='red', facecolor='none', ls='dashed', label=f'Cluster Diameter = {cluster.diameter}',linewidth=1)
    circle_search_region = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), radius=(image.shape[0]/2)*search_arcmin/search_arcmin,
                edgecolor='green', facecolor='none', ls='dashed', label=f'Search Region Diameter = {2*search_arcmin}',linewidth=1.5)
    ax2.add_artist(circle_cluster)
    ax2.add_artist(circle_search_region)
    # Scatter Stars in region
    scatter_stars_in_region, = ax2.plot(stars_in_region_pix_coords[0], stars_in_region_pix_coords[1], 
                                        '.',markersize=2, color='grey',label=f'{len(stars_in_region)} Stars in the region')
    # Scatter dias members
    scatter_dias_members, = ax2.plot(dias_members_pix_coords[0], dias_members_pix_coords[1], '.',markersize=10, color='grey',
                                    markeredgecolor='lightgrey',label=f'{len(cluster.members)} Dias Members')
    # ax2.legend()
    # Scatter my members
    my_members =find_cluster(stars_in_region,refCluster=cluster.name)
    my_members_coords = SkyCoord(ra=my_members['RA_ICRS_1'], dec=my_members['DE_ICRS_1'])
    my_members_pix_coords = wcs.world_to_pixel_values(my_members_coords.ra, my_members_coords.dec)
    scatter_my_members, = ax2.plot(my_members_pix_coords[0],my_members_pix_coords[1],
                                    '.',markersize=5, color='blue', label=f'{len(my_members)} My Members')



    # Read runaways
    runaways_all = cluster.read_table('runaways') #changing this to runaways_all plots all the runaways, not just the >10000K ones
    runaways_all['Source'] = runaways_all['Source'].astype(str)
    runaways_all['Gmag'] = (runaways_all['Gmag'].round(2)).astype(str)
    # runaways_all['RAVE5'] = runaways_all['RAVE5'].astype(str)

    psrs = search_psr(cluster.name)
    n_psr = len(psrs)

    display(psrs)
    for i in range(len(psrs)):
        new_row = {col_name: '0' for col_name in runaways_all.colnames}
        new_row['Source'] = 'PSR'+f'{psrs["JNAME"][i]}'
        new_row['RA_ICRS_1'] = psrs['RAJD'][i]
        new_row['DE_ICRS_1'] = psrs['DECJD'][i]
        new_row['rgeo'] = psrs['DIST'][i]*1000
        new_row['RAVE5'] = f"{psrs['DIST_DM'][i]*1000:.2f}" ###################this takes the dist_dm
        new_row['RAVE6'] = f"{psrs['Separation'][i]:.2f}" ###################this takes the dist_dm

        new_row['pmRA'] = psrs['PMRA'][i]
        new_row['pmDE'] = psrs['PMDEC'][i]
        new_row['rmRA'] = psrs['PMRA'][i] - cluster.all['pmRA']
        new_row['rmDE'] = psrs['PMDEC'][i] - cluster.all['pmDE']

        new_row['Gmag'] = f"{psrs['AGE'][i]:.2f}"+" kyr" ###################this takes the age


        µ_pec = np.sqrt(new_row['rmRA']**2+new_row['rmDE']**2)
        new_row['v_pec'] = µ_pec*psrs['DIST'][i]*4.74

        runaways_all.add_row(new_row)
    runaways_all_coords = SkyCoord(ra=runaways_all['RA_ICRS_1'], dec=runaways_all['DE_ICRS_1'])
    runaways_all_pix_coords = wcs.world_to_pixel_values(runaways_all_coords.ra, runaways_all_coords.dec)



    _ = search_arcmin.round(-1)/2 #round to the nearest 5 (by round first to nearest 10 and then divide by 2
    scalebar_length = ((_.to(u.rad))*(cluster.distance.to(u.m)).to(u.pc)).round(2)
    __ = cluster.distance.value/1000
    add_scalebar(ax2, _, color="yellow", label=f"or {scalebar_length.value}pc (at dist {__:.2f}kpc)", size_vertical=0.5)
    x_min, x_max = ax2.get_xlim()
    y_min, y_max = ax2.get_ylim()
    ax2.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontsize=12,fontweight='bold')

    ax2.set_xlabel('Right Ascension (degrees)')
    ax2.set_ylabel('Declination (degrees)')
    ax2.set_title(f"{cluster.name} with {config['Cluster']['search_extent']}pc search region")
    ax2.grid(color='lightgrey', ls='dotted')


    #######################################
    if len(runaways_all)>0:
        runaway_00, runaway_apdp,runaway_apdm,runaway_amdp,runaway_amdm = [coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=runaways_all['rmRA'],pm_dec=runaways_all['rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']+runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']+runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']+runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']-runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']-runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']+runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                            coord.SkyCoord(ra=runaways_all['RA_ICRS_1'],dec=runaways_all['DE_ICRS_1'], pm_ra_cosdec=(runaways_all['rmRA']-runaways_all['e_rmRA']),pm_dec=runaways_all['rmDE']-runaways_all['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr))]

        earlier_runaway_00 = runaway_00.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_apdp = runaway_apdp.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_apdm = runaway_apdm.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_amdp = runaway_amdp.apply_space_motion(dt = -100_000*u.year)
        earlier_runaway_amdm = runaway_amdm.apply_space_motion(dt = -100_000*u.year)

        earlier_runaway_00_px = wcs.world_to_pixel_values(earlier_runaway_00.ra, earlier_runaway_00.dec)
        earlier_runaway_apdp_px = wcs.world_to_pixel_values(earlier_runaway_apdp.ra, earlier_runaway_apdp.dec)
        earlier_runaway_apdm_px = wcs.world_to_pixel_values(earlier_runaway_apdm.ra, earlier_runaway_apdm.dec)
        earlier_runaway_amdp_px = wcs.world_to_pixel_values(earlier_runaway_amdp.ra, earlier_runaway_amdp.dec)
        earlier_runaway_amdm_px = wcs.world_to_pixel_values(earlier_runaway_amdm.ra, earlier_runaway_amdm.dec)
        runaway_00_px = wcs.world_to_pixel_values(runaway_00.ra, runaway_00.dec)
        # Calculate the displacement vectors
        delta_x = earlier_runaway_apdp_px[0]-earlier_runaway_amdp_px[0]
        delta_y = earlier_runaway_apdp_px[1]-earlier_runaway_apdm_px[1]

        for pxx, pxy, dx, dy in zip(earlier_runaway_00_px[0], earlier_runaway_00_px[1], delta_x, delta_y):
            ellipse = patches.Ellipse(
                (pxx, pxy),
                width=1.5*dx,
                height=1.5*dy,
                fill=True,
                color='g',
                alpha=0.2
            )
            ax2.add_patch(ellipse)

        for c in [earlier_runaway_apdp_px,earlier_runaway_apdm_px,earlier_runaway_amdp_px,earlier_runaway_amdm_px]:
            delta_x = c[0] - runaway_00_px[0]
            delta_y = c[1] - runaway_00_px[1]
            # Draw the vectors
            ax2.quiver(runaway_00_px[0], runaway_00_px[1], delta_x, delta_y, angles='xy', scale_units='xy', scale=1, color='limegreen', width=0.001)
        ############################
    else:
        print("No runaways found")
    # Scatter runaways
    scatter_runaways_all = ax2.scatter(runaways_all_pix_coords[0], runaways_all_pix_coords[1], s=30,picker=True,pickradius=20,
                                    c=runaways_all['Temp. Est'], cmap='spring_r', norm=plt.Normalize(4000, 23000), 
                                    label=f'{len(runaways_all)-n_psr} Runaway(s)')
    if n_psr>0:
        scatter_psrs = ax2.scatter(runaways_all_pix_coords[0][-n_psr:], runaways_all_pix_coords[1][-n_psr:], s=20,
                                        c='cyan', 
                                        label=f'{n_psr} Pulsar(s)')
        if psr_circles:
            for ra_pix_psr,dec_pix_psr,ra_psr,dec_psr in zip(runaways_all_pix_coords[0][-n_psr:], runaways_all_pix_coords[1][-n_psr:],
                                            runaways_all_coords[-n_psr:].ra, runaways_all_coords[-n_psr:].dec):
                print(ra_pix_psr,dec_pix_psr)
                center = SkyCoord(ra_psr, dec_psr, unit='deg')
                t = 100*u.kyr
                vs_psr = [300,400,500,600]*u.km/u.s
                for v_psr in vs_psr:
                    v_psr_arcmin = np.arctan(v_psr.to(u.pc/u.kyr)*t/cluster.distance)
                    radius = v_psr_arcmin.to(u.arcmin)
                    sky_reg = CircleSkyRegion(center, radius)
                    pix_reg = sky_reg.to_pixel(wcs)
                    circle_v_spr = patches.Circle((ra_pix_psr,dec_pix_psr),radius=pix_reg.radius,edgecolor='azure', facecolor='aqua',alpha = 0.1)
                    ax2.add_patch(circle_v_spr)

    global annotations
    annotations = []

    def on_release(event):
        global annotations, table2
        try:
            if table2:
                table2.remove()
        except:
            pass
        fig.canvas.draw() # for the release to remove functionality

    def onpick3(event):
        global annotations, table2
        
        ind = event.ind
        
        # print(ind,type(ind),len(ind),event.artist)
        # print(result)
        import matplotlib
        if isinstance(event.artist, matplotlib.collections.PathCollection) and len(ind) == 1:
            for index in ind:
                index = ind[0]
                temp_est = runaways_all['Temp. Est'][index]
                # ann = ax2.annotate(f'{temp_est:,.0f} K', (0,0), xytext=(0, 0), textcoords='offset points', fontsize=12, color='black', ha='center', fontweight='bold')
                # ann.set_path_effects([pe.withStroke(linewidth=4, foreground='white')])
                # annotations.append(ann)
                # Add the table below the annotation
                if (runaways_all[index]["Source"][:3] == 'PSR'): #then this is a psr, table has different format
                    table_data = [
                        ['JNAME', f'{runaways_all[index]["Source"]}'],
                        ['Separation', f'{runaways_all[index]["RAVE6"]} arcmin'],
                        ['Coord.', f'{runaways_all["RA_ICRS_1"][index]:.4f} {"+" if runaways_all["DE_ICRS_1"][index] >= 0 else ""}{runaways_all["DE_ICRS_1"][index]:.4f}'],
                        ['Age', f'{runaways_all["Gmag"][index]}'], #this contains the age
                        ['Dist.', f'{runaways_all["rgeo"][index]:.0f}'+' pc'], #this contains the distance
                        ['Dist_DM', f'{runaways_all["RAVE5"][index]}'] #this contains the dist_dm #dispersion measured distance
                        # Add more rows with actual values
                    ]
                else: #it is a runaway star
                    # Add the table below the annotation
                    table_data = [
                        ['SourceID', f'{runaways_all[index]["Source"]}'],
                        ['Temp. Est.', f'{runaways_all[index]["Temp. Est"]:,.0f} K'],
                        ['Coord.', f'{runaways_all["RA_ICRS_1"][index]:.4f} {"+" if runaways_all["DE_ICRS_1"][index] >= 0 else ""}{runaways_all["DE_ICRS_1"][index]:.4f}'],
                        ['Gmag', f'{runaways_all["Gmag"][index]}'+' mag'],
                        ['Dist.', f'{runaways_all["rgeo"][index]:.0f}'+'$\pm$'+f'{(runaways_all["rgeo"][index]-runaways_all["b_rgeo"][index]):.0f}'+' pc']
                        # Add more rows with actual values
                    ]

                table_bbox = [0.0, 0.0, 0.45, 0.2]  # [left, bottom, width, height]
                table2 = ax2.table(cellText=table_data, cellLoc='right', loc='lower center', bbox=table_bbox)
                table2.auto_set_column_width(col=[0,1])
                for key, cell in table2._cells.items():
                    cell.set_linewidth(0.5)  # Set the border width
                    cell.set_edgecolor('lightgray')  # Set the border color
                    # cell.get_text().set_weight('bold')
            
            fig.canvas.draw()

    legend = ax2.legend(loc='upper right')
    scatter_stars_in_region_legend, scatter_dias_members_legend,scatter_my_members_legend = legend.get_lines()
    scatter_stars_in_region_legend.set_picker(True)
    scatter_stars_in_region_legend.set_pickradius(10)
    scatter_dias_members_legend.set_picker(True)
    scatter_dias_members_legend.set_pickradius(10)
    scatter_my_members_legend.set_picker(True)
    scatter_my_members_legend.set_pickradius(10)

    graphs = {}
    graphs[scatter_stars_in_region_legend] = scatter_stars_in_region
    graphs[scatter_dias_members_legend] = scatter_dias_members
    graphs[scatter_my_members_legend] = scatter_my_members


    colorbar = fig.colorbar(scatter_runaways_all,ax=ax2)
    # Connect the pick_event to the onpick3 function
    fig.canvas.mpl_connect('pick_event', onpick3)
    fig.canvas.mpl_connect('button_release_event', on_release)

    def on_pick(event):
        legend = event.artist
        isVisible = legend.get_visible()
        try:
            graphs[legend].set_visible(not isVisible)
            legend.set_visible(not isVisible)

        except:
            pass
        # fig.canvas.draw() # for the release to remove functionality

    plt.connect('pick_event', on_pick)
    ax2.set_facecolor('black')

    fits_file.close()
    if save:
        plt.savefig(f'{cluster.name}/{cluster.name}_traceback.{save}')
    


def search_psr(cluster,extra=config['psr_extra'],radial_tolerance=config['psr_radial_tolerance']):
    if isinstance(cluster,str):
        cluster = Cluster(cluster)
        print(f'{"Searching pulsars near":->50}'+f' {cluster.name:-<50}')

    elif isinstance(cluster,Cluster):
        cluster=cluster
        print(f'{"Searching pulsars near":->50}'+f' {cluster.name:-<50}')


    print(f'({extra} pc from cluster edge, {radial_tolerance*100:.0f}% radial tolerance around distance {cluster.distance.to(u.kpc)})')

    cluster_coord_ra_str = cluster.coordinates.ra.to_string(unit='hourangle', sep=':', precision=3, pad=True)[0]
    cluster_coord_dec_str = cluster.coordinates.dec.to_string(unit='degree', sep=':', precision=3, pad=True)[0]

    psr_search_deg = cluster.calculate_search_arcmin(extra=extra).to(u.deg).value
    c = [cluster_coord_ra_str,cluster_coord_dec_str,psr_search_deg]
    query = QueryATNF(params=['NAME','JNAME','RAJD','DECJD','DIST','DIST_DM','AGE','PMRA','PMDEC','S400','ASSOC','AGE_I','PX'], circular_boundary=c)
    print(f'{len(query.table)} Pulsar(s) within {psr_search_deg:.2f} degrees ({extra} pc)')

    r_close = (1-radial_tolerance)*cluster.distance.to(u.kpc).value
    r_far = (1+radial_tolerance)*cluster.distance.to(u.kpc).value

    dist_filtered_psr_table = query.table[(query.table['DIST']>r_close) & (query.table['DIST']<r_far)]
    psr_coords = SkyCoord(ra = dist_filtered_psr_table['RAJD'], dec = dist_filtered_psr_table['DECJD'],pm_ra_cosdec = dist_filtered_psr_table['PMRA'],pm_dec = dist_filtered_psr_table['PMDEC'])
    psr_sep = psr_coords.separation(cluster.coordinates).to(u.arcmin)
    dist_filtered_psr_table.add_column(psr_sep,name='Separation')
    dist_filtered_psr_table.sort('Separation')
    dist_filtered_psr_table['AGE'] = dist_filtered_psr_table['AGE'].to(u.kyr)
    dist_filtered_psr_table['AGE_I'] = dist_filtered_psr_table['AGE_I'].to(u.kyr)
    print(f'Of which {len(dist_filtered_psr_table)} Pulsar(s) from {r_close:.2f} kpc to {r_far:.2f} kpc')
    # display(dist_filtered_psr_table['Separation','NAME','JNAME','RAJD','DECJD','DIST','DIST_DM','AGE','PMRA','PMDEC','S400','ASSOC','AGE_I','PX'])

    return dist_filtered_psr_table['Separation','NAME','JNAME','RAJD','DECJD','DIST','DIST_DM','AGE','PMRA','PMDEC','S400','ASSOC','AGE_I','PX']
    #add distance from cluster column




#plot the find_members and dias_members on the search region fits
#find out the reason why some clusters were not searched. probably because of the search distance
#3D visualisation
#unique config files for each cluster. add string on top of it indicating last change date
