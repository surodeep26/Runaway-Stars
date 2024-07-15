import os
import astropy
from astropy.table import Table, Column, QTable, join,vstack
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

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import fits
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
matplotlib.rcParams.update({'font.size': 12})
def simbad(keyword:str):
    return IFrame(f"https://simbad.cds.unistra.fr/simbad/sim-basic?Ident={keyword}+&submit=SIMBAD+search", 1400,500)



cluster_list = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
cluster_list_mod = Table.read('Clusters_from_dias_and_a99_mod', format='ascii.ecsv')
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
    from typing import List

    def __init__(self, cluster_name,version='dr3', add_members: List[int]=[], Pmemb = config['memb_prob'], PlxQuality = config['parallax_quality_threshold']):
        '''
        add_members: Gaia sourceIDs of new members to be added to the cluster, only works with the default `version='dr3'`
        '''
        # Assuming 'cluster_list' is a predefined list containing your Astropy table
        cluster_list = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
        cluster_list_mod = Table.read('Clusters_from_dias_and_a99_mod', format='ascii.ecsv')
        if cluster_list[cluster_list['Cluster'] == cluster_name][0] == cluster_list_mod[cluster_list_mod['Cluster'] == cluster_name][0]:
            cluster_list = cluster_list
            print('Original Data')
        else:
            cluster_list = cluster_list_mod
            print('Modified Data')

        self.cluster_table = cluster_list[cluster_list['Cluster'] == cluster_name]
        self.all = self.cluster_table[0]
        self.name = self.cluster_table['Cluster'][0]
        ti = (Time('J2000')+1*u.Myr)
        self.coordinates = SkyCoord(ra=self.cluster_table['RA_ICRS'],dec=self.cluster_table['DE_ICRS'],pm_ra_cosdec=self.cluster_table['pmRA'],pm_dec=self.cluster_table['pmDE'],obstime=ti)[0]
        self.diameter = self.cluster_table['Diameter'][0]*u.arcmin
        self.diameter_dias = (self.cluster_table['r50'][0]*u.deg*2).to(u.arcmin)
        self.dias_distance = self.cluster_table['Dist'][0]*u.pc
        self.Pmemb = Pmemb
        self.PlxQuality = PlxQuality

        if len(self.cluster_table) == 0:
            raise ValueError(f"No data found for cluster '{cluster_name}' in the cluster list.")
        if version == 'dr3':
            d = self.dias_members(memb_prob=Pmemb)
            
            dm = d[d['Gmag']<config['gmag_lim']]['Source','RAdeg','DEdeg','Gmag','e_Gmag','Plx','e_Plx','Pmemb']
            sources = np.array(dm['Source'])
            sources = np.insert(sources,0,np.array(add_members))

            sources_int_list = sources.astype(np.int64)
            sr=self.stars_in_region(no_rmRArmDE=True)
            mruwe = (sr['RUWE']<config['ruwe_lim'])
            madd = np.isin(sr['Source'],np.array(add_members))
            sr = sr[mruwe | madd]
            mask = np.isin(sr['Source'], sources_int_list)
            # print(sources_int_list)
            dr3dias = sr[mask]
            print(f"{len(dr3dias)}/{len(d)} members matched in dr3 search region,{len(d)-len(dm)} have > {config['gmag_lim']} Gmag")
            # dr3dias.add_row()
            dr3dias.sort('Gmag')

            mask_plx = [value >= PlxQuality for value in (dr3dias['Plx']/dr3dias['e_Plx'])]
            dr3dias = dr3dias[mask_plx]
            print(f"{len(dr3dias)}/{len(dr3dias)+(len(mask_plx)-np.count_nonzero(mask_plx))} members have good parallax")

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
            self.all['N'] = len(dr3dias)
        elif version == 'dr2':
            dr3dias = self.dias_members(memb_prob=Pmemb)
            mask_plx = [value >= PlxQuality for value in (dr3dias['Plx']/dr3dias['e_Plx'])]
            dr3dias = dr3dias[mask_plx]
            self.N = self.cluster_table['N'][0] #no. of members
        elif version == 'dr2raw':
            dr3dias = self.dias_members(memb_prob=0)
            self.N = self.cluster_table['N'][0] #no. of members

        # Extracting values from the table and setting attributes
        if not os.path.exists(f'{self.name}'):
            os.mkdir(f'{self.name}')
        self.distance = self.all['Dist']*u.pc
        self.RV = self.cluster_table['RV'][0]*u.km/u.s
        self.Av = self.cluster_table['Av'][0]
        self.logage = self.cluster_table['logage'][0]
        self.age = 10**(self.cluster_table['logage'][0])
        self.FeH = self.cluster_table['__Fe_H_'][0]
        self.members = dr3dias #gets dias_members with the default settings for memb_prob and parallax_threshold_quality
        #to change the memb_prob and parallax_threshold_quality filters for member selection, use something like dias_members(0.6,11)
        self.N = len(dr3dias)
        self.runaways = self.read_table('runaways')
        self.runaways_all = self.read_table('runaways_all')
    
    def changeParam(self,change: tuple):
        param,new_value = change

        _old_value  = cluster_list[cluster_list['Cluster'] == self.name][0][param]
        cluster_list_mod[param][cluster_list_mod['Cluster'] == self.name] = new_value
        cluster_list_mod.write('Clusters_from_dias_and_a99_mod', format='ascii.ecsv', overwrite=True)
        print(f'Changed {param} from {_old_value} --> {new_value}')
    

    def restoreParam(self, param: str):
        # Read the original and modified tables
        cluster_list_original = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
        cluster_list_modified = Table.read('Clusters_from_dias_and_a99_mod', format='ascii.ecsv')
        # Get the name of the cluster
        cluster_name = self.name
        if param == 'all':
            # Iterate over all columns in the original table
            for column in cluster_list_original.colnames:
                # Get the original value for the current column
                original_value = cluster_list_original[cluster_list_original['Cluster'] == cluster_name][0][column]
                # Get the current value for the current column in the modified table (for logging)
                current_value = cluster_list_modified[cluster_list_modified['Cluster'] == cluster_name][0][column]
                # Restore the original value in the modified table
                cluster_list_modified[column][cluster_list_modified['Cluster'] == cluster_name] = original_value
                if current_value != original_value:
                    print(f'Restored {column}: current {current_value} --> {original_value}')
        else:
            # Handle the case for a single parameter as before
            original_value = cluster_list_original[cluster_list_original['Cluster'] == cluster_name][0][param]
            current_value = cluster_list_modified[cluster_list_modified['Cluster'] == cluster_name][0][param]
            
            cluster_list_modified[param][cluster_list_modified['Cluster'] == cluster_name] = original_value
            print(f'Restored {param}: current {current_value} --> {original_value}')
        # Save the modified table
        cluster_list_modified.write('Clusters_from_dias_and_a99_mod', format='ascii.ecsv', overwrite=True)

    
    def dias_members(self,memb_prob=0):
        dias_members = Table.read(f'./Clusters_Dias/{self.name}.dat', format='ascii.tab')
        # dias_members = Table.read(f'Clusters/{self.name}.dat', format='ascii.tab')
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
        # mask2 = [value >= parallax_quality_threshold for value in (dias_members['Plx']/dias_members['e_Plx'])]
        # final_mask = (np.array(mask1) & np.array(mask2))
        final_mask = (np.array(mask1))
        dias_members = dias_members[final_mask]
        #create a filter to get rid of the stars which don't have a Bp-Rp value in the .dat table
        mask_bprp = [value is not None for value in dias_members["BP-RP"] ]
        plottable_members = dias_members[mask_bprp]
        print(len(dias_members)-len(plottable_members),"memebrs do not have BP-RP from Dias.")

        return plottable_members

    def clean(self,what='except_downloads'):
        print(f"Deleting: {what}")
        folder_path = f'{self.name}'

        # List all files in the folder
        files = os.listdir(folder_path)
        
        # Iterate over each file and delete it
        for file_name in files:
            file_path = os.path.join(folder_path, file_name)
            if what=='everything' and os.path.isfile(file_path):
                os.remove(file_path)
                print(f'{file_path} removed')

            if what=='runaways' and os.path.isfile(file_path) and 'runaway' in file_path:
                print(f'{file_path} removed')
                os.remove(file_path)
            if what=='except_downloads' and os.path.isfile(file_path) and not (('compare_data' in file_path) or ('stars_from' in file_path) or ('extra' in file_path)):
                print(f'{file_path} removed')
                os.remove(file_path)


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
        D = self.dias_distance  # Assuming you have a 'distance' attribute in your Cluster class
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
        ax.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontweight='bold')
        ax.legend()
        images[0].close()  # Close the opened FITS file

        # Show the plot
        plt.show()
    

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

            _l = (1-config['distance_tolerance'])*self.dias_distance.value
            _u = (1+config['distance_tolerance'])*self.dias_distance.value
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
                                            'b_rgeo','B_rgeo','FG','e_FG','FBP','e_FBP','FRP','e_FRP','RAVE5','RAVE6']
            # Adding row for uncertainties in Gmag and BPmag and RPmag
            # values for Gaia G, G_BP, G_RP zero point uncertainties
            sigmaG_0 = 0.0027553202
            sigmaGBP_0 = 0.0027901700
            sigmaGRP_0 = 0.0037793818
            stars_in_region['e_Gmag'] = np.sqrt((-2.5/np.log(10)*stars_in_region['e_FG']/stars_in_region['FG'])**2 + sigmaG_0**2)
            stars_in_region['e_BPmag'] = np.sqrt((-2.5/np.log(10)*stars_in_region['e_FBP']/stars_in_region['FBP'])**2 + sigmaGBP_0**2)
            stars_in_region['e_RPmag'] = np.sqrt((-2.5/np.log(10)*stars_in_region['e_FRP']/stars_in_region['FRP'])**2 + sigmaGRP_0**2)
            stars_in_region['e_BP-RP'] = stars_in_region['e_BPmag']+stars_in_region['e_RPmag']
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
        - tables (str): Variable number of table names to be read. The name of the table is the text after the last "_". 
        For eg., "NGC_4103_psrs.tsv" is called "psrs"

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
        ax.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)")
        ax.set_ylabel(r"$G$ (mag)")
        ax.set_title(f"CMD for {self.name}")
        ax.invert_yaxis()
        ax.legend()

    def latex_text(self):
        ra, dec = self.coordinates.ra, self.coordinates.dec
        ra_str, dec_str = ra.to_string(format='latex')[1:-1], dec.to_string(format='latex')[1:-1]
        dist_str = str(self.all['Dist'])+r"\pm"+str(self.all['e_Dist'])
        print(rf"{self.name.replace('_', ' ')} is located at $\alpha = {ra_str}, \delta = {dec_str}$ at a distance of ${dist_str}$ pc.")



    def latex_table_kinematics(self):
        t_dict = {
            # 'Name':[self.name.replace('_',' ')],
            r'$\alpha (^\circ)$': [f"{(self.coordinates.ra):.2f}".replace(' deg','')],
            r'$\delta (^\circ)$': [f"{(self.coordinates.dec):.2f}".replace(' deg','')],
            "Diam. (')":[self.diameter.value],
            "N":[len(self.members)],
            r"$\mu_{\alpha}^*$ (mas/yr)":[f"{self.all['pmRA']:.2f}"+"$\pm$"+f"{self.all['e_pmRA']:.2f}"],
            r"$\mu_{\delta}$ (mas/yr)":[f"{self.all['pmDE']:.2f}"+"$\pm$"+f"{self.all['e_pmDE']:.2f}"],
            r"$v_R$ (km/s)":[f"{self.all['RV']:.2f}"+"$\pm$"+f"{self.all['e_RV']:.2f}"] if not isinstance(self.all['RV'],np.ma.core.MaskedConstant) else ['N/A'],
            r"$d$ (pc)":[f"{self.all['Dist']:.0f}"+"$\pm$"+f"{self.all['e_Dist']:.0f}"],

        }
        latexdict={'tabletype':'table*','tablealign':'h',
                   'header_start':r'\hline','data_start':r'\hline','data_end': r'\hline', 
                   'caption':f'Kinematic parameters of {self.name.replace("_"," ")}',
                   'preamble':'\label{tab:'+f'{self.name}-kinematics'+'}'}
        latexdict['tablehead'] = r'head'
        # Convert the dictionary to an Astropy table
        astropy_table = Table(t_dict)
        astropy.io.ascii.write(astropy_table, format='latex',output='text.txt',overwrite=True,
                    latexdict=latexdict)
        with open('text.txt', 'r+') as f:
            lines = f.readlines()
            lines[1], lines[2] = lines[2], lines[1]
            f.seek(0)
            f.writelines(lines)
            print(''.join(lines))

        os.remove('text.txt')
    
    def latex_table_members(self,n=None):
        """
        Prints a LaTeX table for the cluster members of the cluster object. The members are the same as what you 
        get from `cluster.members`. Table is sorted according to the Gmag of the stars.
        Parameters:
        - n (int): number of cluster members to be returned
        Example usage:
        cluster = Cluster('IC_2395')
        latex_table_members(cluster,n=10)
        """
        table = self.members
        table.sort('Gmag')
        table = table[:n]
        # Create an empty list to store the formatted distances
        formatted_distances = []
        formatted_parallaxes = []
        formatted_gmags = []
        formatted_bprps = []
        # Iterate over each row in the table
        for row in table:
            distance = row['rgeo']
            plx = row['Plx']
            e_plx = row['e_Plx']
            upper_error = row['B_rgeo']-row['rgeo']
            lower_error = -row['b_rgeo']+row['rgeo']
            gmag = row['Gmag']
            e_gmag = row['e_Gmag']
            bprp = row['BP-RP']
            e_bprp = row['e_BP-RP']
            # Format strings
            formatted_distance = f"${distance:.0f}^{{+{upper_error:.0f}}}_{{-{lower_error:.0f}}}$"
            formatted_parallax = f"${plx:.4f}\pm{e_plx:.4f}$"
            formatted_gmag = f"${gmag:.3f}\pm{e_gmag:.3f}$"
            formatted_bprp = f"${bprp:.3f}\pm{e_bprp:.3f}$"
            # Append the formatted distance to the list
            formatted_distances.append(formatted_distance)
            formatted_parallaxes.append(formatted_parallax)
            formatted_gmags.append(formatted_gmag)
            formatted_bprps.append(formatted_bprp)

        # Add a new column to the table with the formatted distances
        table['formatted_distance'] = formatted_distances
        table['formatted_parallax'] = formatted_parallaxes
        table['formatted_gmag'] = formatted_gmags
        table['formatted_bprp'] = formatted_bprps
        t_dict = {
            'Gaia DR3 Source':table['Source'],
            r"$r_{\text{geo}}$ (pc)":table['formatted_distance'],
            r"$\pi$ (mas)":table['formatted_parallax'],
            r"$G$ (mag)":table['formatted_gmag'],
            r"$G_{\text{BP}}-G_{\text{RP}}$ (mag)":table['formatted_bprp']
            # r"$v_R$ (km/s)":[f"{self.all['RV']:.2f}"+"$\pm$"+f"{self.all['e_RV']:.2f}"] if not isinstance(self.all['RV'],np.ma.core.MaskedConstant) else ['N/A'],
            # r"$d$ (pc)":[f"{self.all['Dist']:.0f}"+"$\pm$"+f"{self.all['e_Dist']:.0f}"],

        }
        latexdict={'tabletype':'table*','tablealign':'h',
                    'header_start':r'\hline','data_start':r'\hline','data_end': r'\hline', 
                    'caption':f'Selected members of {self.name.replace("_"," ")}',
                    'preamble':'\label{tab:'+f'{self.name}-members'+'}'}
        latexdict['tablehead'] = r'head'
        # Convert the dictionary to an Astropy table
        astropy_table = Table(t_dict)
        astropy.io.ascii.write(astropy_table, format='latex',output='text.txt',overwrite=True,
                    latexdict=latexdict)
        with open('text.txt', 'r+') as f:
            lines = f.readlines()
            lines[1], lines[2] = lines[2], lines[1]
            f.seek(0)
            f.writelines(lines)
            print(''.join(lines))

        os.remove('text.txt')



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



def theoretical_isochrone(cluster,Av=None,logage=None,FeH=None,output=None, printing = True):
    metallicity = FeH 
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
    else:
        Av = float(Av)
    if logage is None:
        logage = cluster.all['logage']
    else:
        logage = float(logage)
    if metallicity is None:
        metallicity = cluster.all['__Fe_H_']
    else:
        metallicity = float(metallicity)

    #Check if parameters are unchanged
    if (str(Av) == str(cluster.all['Av'])) and (str(logage) == str(cluster.all['logage'])) and (str(metallicity) == str(cluster.all['__Fe_H_'])):
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
        print(f"getting isochrone form cmd3.7 with Av:{Av:.2f}, logage:{logage:.2f}, metallicity:{metallicity:.2f}")
        options = webdriver.ChromeOptions()
        browser = webdriver.Chrome(service=s, options=options)
        browser.get('http://stev.oapd.inaf.it/cgi-bin/cmd')

        #Evolutionary Tracks #from config
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[3]/td[1]/input[1]").click() #PARSEC version 2.0
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[5]/td/input").click() #+ COLIBRI S_37
        #Phtotometric System #from config
        photometricSystem = Select(browser.find_element(By.XPATH,"//select[@name='photsys_file']")) #dropdown list for available photometric systems
        photometricSystem.select_by_value("YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat") # Gaia EDR3 bands
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[2]/table/tbody/tr[5]/td[1]/input").click() # VBC +new Vega for PARSEC 2.0 #As above, but adopting revised SED for Vega from Bohlin et al. (2020) (namely CALSPEC alpha_lyr_stis_010.fits).
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


def get_runaways(cluster,fs,theoretical_data,separation_factor=2,dist_filter_factor=1):
    """
    Selects stars from the input astropy table fs which are runaway star candidates by tracing them back to the cluster.
    If the stars happens to be in the cluster in the past 100kyr then it is selected to be a candidate.
    Then the temperature of the star is estimated by comparing its color (Bp-Rp magnitude) with the theoretical isochrone that is provided.

    Parameters:
    -    cluster (Cluster class object): The cluster compared whose runaways are to be found.
    -    fs or nsfs (astropy table): A table of the fast stars (relative to the cluster mean proper motion) in the region (`search_arcmin`) 
        around the cluster centre. Can also take `nsfs` and trace them back with the same algorithm.
    -   theoretical_data (astropy table ascii format downloaded from CMD): theoretical isochrone to estimate the temperatures of the runaways (or walkaways)
    - separation_factor: default 2 (since separation<cluster diameter/2 implies within cluster in the past). indicates what is the minimum separation factor to be considered inside cluster.
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
        threshold_separation = cluster.diameter.value/separation_factor
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
            new_star_gmag = fs[fs['Source']==source_id]['Gmag']
            # Get the BP-RP and Temperature columns from the table
            temperature_column = theoretical_data['Teff0']
            bp_rp_column = theoretical_data['BP-RP']
            gmag_column = theoretical_data['Gmag']
            # Calculate the differences between the new star's BP-RP value and all values in the table
            differences_bp_rp = abs(bp_rp_column - new_star_bp_rp)
            # Calculate the differences between the new star's Gmag value and all values in the table
            differences_gmag = abs(gmag_column - new_star_gmag)
            # iso_dist = differences_bp_rp**2+differences_gmag**2 #method 1
            iso_dist = differences_bp_rp #method 2
            # Find the index of the star with the closest BP-RP value
            closest_star_index = np.argmin(iso_dist)
            # Get the temperature of the closest star
            closest_star_temperature = temperature_column[closest_star_index]
            selected_stars_dict_with_temp[source_id] = (selected_stars_dict[source_id],closest_star_temperature)


        sources_to_mask = list(selected_stars_dict_with_temp.keys())
        mask = [source in sources_to_mask for source in fs['Source']]
        run = fs[mask]
        run.add_column(np.array(list(selected_stars_dict_with_temp.values()),dtype=object)[:,1],name='Temp. Est',index =0)
        run['Temp. Est'] = run['Temp. Est']*u.K
        run.sort('Temp. Est',reverse=True)
        # Save
        def check_intersection(range1, range2):
            return max(range1[0], range2[0]) <= min(range1[1], range2[1])

        # Given range
        if dist_filter_factor==None:
            range2 = (0, 
                      15000)
        else:
            range2 = (cluster.all['Dist']-dist_filter_factor*cluster.all['e_Dist'], 
                      cluster.all['Dist']+dist_filter_factor*cluster.all['e_Dist'])

        mask = [check_intersection((row['b_rgeo'], row['B_rgeo']), range2) for row in run]

        mask = np.array(mask)
        # Select stars which are not likele foreground or background
        run= run[mask]


        run.write(file_path,format='ascii.ecsv',overwrite=True)
        run.to_pandas().to_excel(os.path.join(cluster.name,f'{cluster.name}_{name}_all.xlsx'), index=False)

        mask = [T > config['runaway_temp'] for T in run['Temp. Est']]
        runaways = run[mask]
        runaways.write(file_path.replace("_all",""),format='ascii.ecsv',overwrite=True)
        runaways.to_pandas().to_excel(os.path.join(cluster.name,f'{cluster.name}_{name}.xlsx'), index=False)


        return run

def get_coord(runaway):
    return coord.SkyCoord(ra=runaway['RA_ICRS_1']*u.deg,dec=runaway['DE_ICRS_1']*u.deg, pm_ra_cosdec=runaway['rmRA']*u.mas/u.year,pm_dec=runaway['rmDE']*u.mas/u.year, frame='icrs')


def test_isochrones(cluster):
    folder_path = f'{cluster.name}'
    # List all files in the folder
    files = os.listdir(folder_path)
    # Initialize list to store dictionaries
    result = []
    # Iterate over each file
    for file_name in files:
        file_path = os.path.join(folder_path, file_name)
        if "compare_data_out_" in file_path:
            st = file_path.replace('.dat', '').split("_")[-3:]
            # Initialize dictionary to store extracted values
            extracted_values = {}
            # Iterate through the strings
            for string in st:
                if 'Av' in string:
                    extracted_values['Av'] = float(string[2:])
                elif 'age' in string:
                    extracted_values['logage'] = float(string[3:])
                elif 'FeH' in string:
                    extracted_values['FeH'] = float(string[3:])
            # Append dictionary to result list
            result.append(extracted_values)
    return result



def plot_cmd(cluster,save=False,multiple=False,**kwargs):

    # Plot CMD
    BP_RP_theo, Gmag_theo = theoretical_isochrone(cluster,printing=False)
    if not multiple:
        BP_RP_theo, Gmag_theo = theoretical_isochrone(cluster,printing=False,**kwargs)


    fig,ax1 = plt.subplots(figsize=(10, 8.8))

    lines = []
    line, = ax1.plot(BP_RP_theo, Gmag_theo, label='Dias Theoretical Isochrone (for Teff.)', alpha=0.99)
    lines.append(line)

    iso_4_temp0 = f'Av:{str(cluster.Av)}, logage:{str(cluster.logage)}, FeH:{str(cluster.FeH)}'
    # Split the input string into key-value pairs
    pairs = iso_4_temp0.split(', ')
    # Update the values using the dictionary
    for i, pair in enumerate(pairs):
        key, value = pair.split(':')
        if key in kwargs:
            pairs[i] = f"{key}:{kwargs[key]}"
    # Reconstruct the string
    iso_4_temp = ', '.join(pairs)
    print('iso_4_temp0', iso_4_temp0)
    print('iso_4_temp', iso_4_temp)
    if ((not multiple) and (iso_4_temp0!=iso_4_temp)) or (cluster_list[cluster_list['Cluster'] == cluster.name][0] != cluster_list_mod[cluster_list_mod['Cluster'] == cluster.name][0]):
        lines[0].set_label('Dias Theoretical Isochrone (modified)')
    if (not multiple) and (iso_4_temp0==iso_4_temp):
        lines[0].set_label('Dias Theoretical Isochrone')
    ax1.set_ylim(bottom= min(Gmag_theo)-4,top=18)
    ax1.set_xlim(left=min(BP_RP_theo)-0.5, right=max(BP_RP_theo))
    if multiple:
        all_isochrones = test_isochrones(cluster)

        for isoc in all_isochrones:
            BP_RP_theo, Gmag_theo = theoretical_isochrone(cluster,**isoc, printing=False)
            _lbl =  f'{isoc}'.replace("{'Av': ", "Av:").replace(", 'logage': ", ", logage:").replace(", 'FeH': ", ", FeH:").replace("}", "")
            print(_lbl)
            if _lbl == iso_4_temp:
                print(_lbl)
                _lbl += ' (for Teff.)'
                lines[0].set_label('Dias Theoretical Isochrone')
            line, = ax1.plot(BP_RP_theo, Gmag_theo, label=_lbl,alpha=1)
            lines.append(line)

    ax1.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)")
    ax1.set_ylabel(r"$G$ (mag)")
    ax1.set_title(f"CMD for {cluster.name}")
    ax1.invert_yaxis()
    #Scatter stars in the region
    sir = cluster.stars_in_region()
    sir_gmag,sir_bp_rp = sir['Gmag'],sir['BP-RP']
    ax1.scatter(sir_bp_rp,sir_gmag,s=2, color='grey',zorder=1,label=f'{len(sir)} stars in the region {cluster.calculate_search_arcmin()}')
    # Scatter dias members
    dias_members = cluster.members
    dias_gmag,dias_bp_rp = dias_members['Gmag'],dias_members['BP-RP']
    ax1.errorbar(dias_bp_rp,dias_gmag, color='black',zorder=2,fmt='o',xerr=dias_members['e_BP-RP']+0.02,yerr=dias_members['e_Gmag'],
                 label=rf'{len(dias_members)} Dias Members P $\geq$ {cluster.Pmemb}, Q $\geq$ {cluster.PlxQuality}')
    # Table for cluster parameters
    cluster_table = [
        ['Members',len(dias_members)],
        [r'Metallicity $[Fe/H]$',cluster.FeH],
        ['Logage',cluster.logage],
        ['Av',cluster.Av],
        ['Distance (pc)',str(round(cluster.distance.value))+"$\pm$"+f'{cluster.all["e_Dist"]}'+'pc']
    ]

    # Update Metallicity if changed
    if 'FeH' in kwargs and kwargs['FeH'] != cluster.FeH:
        cluster_table[1][1] = f'{cluster.FeH:.2f} --> {kwargs["FeH"]}'

    # Update Logage if changed
    if 'logage' in kwargs and kwargs['logage'] != cluster.logage:
        cluster_table[2][1] = f'{cluster.logage:.2f} --> {kwargs["logage"]}'

    # Update Av if changed
    if 'Av' in kwargs and kwargs['Av'] != cluster.Av:
        cluster_table[3][1] = f'{cluster.Av:.2f} --> {kwargs["Av"]}'


    table_bbox = [0.0, 0.84, 0.44, 0.16]  # [left, bottom, width, height]
    table = ax1.table(cellText=cluster_table, cellLoc='right', loc='upper left',bbox=table_bbox)

    for key, cell in table._cells.items():
        cell.set_linewidth(0.5)  # Set the border width
        cell.set_edgecolor('lightgray')  # Set the border color
    # Scatter runaways
    runaways_all = cluster.read_table('runaways_all') #plot all runaways
    # runaways_all = cluster.read_table('runaways') #plot main runaways
    #recalculate temperatures for different theoretical isochrone
    td = theoretical_isochrone(cluster,output='table',**kwargs,printing=False)
    print(f"calculating temperatures based on: {kwargs}")
    Ttheo = td['Teff0']
    bprptheo = td['BP-RP']
    gmagtheo = td['Gmag']
    for star in runaways_all:
        differences_bprp = abs(bprptheo - star['BP-RP'])
        differences_gmag = abs(gmagtheo - star['Gmag'])
        # differences = differences_bprp**2+differences_gmag**2 #method 1
        differences = differences_bprp #method 2
        closest_star_index = np.argmin(differences)
        new_closest_star_temperature = Ttheo[closest_star_index]
        star['Temp. Est']=new_closest_star_temperature

    _runaways_gmag,_runaways_bp_rp = runaways_all['Gmag'],runaways_all['BP-RP']
    mask = [T > config['runaway_temp'] for T in runaways_all['Temp. Est']]
    runaways = runaways_all[mask]
    runaways_gmag,runaways_bp_rp = runaways['Gmag'],runaways['BP-RP']
    ## all runaways
    try:
        label_all_run = rf'{len(runaways_all)} Runaway(s) $T_{{max}}$ = {max(runaways_all["Temp. Est"]):,.0f} K'
    except:
        label_all_run = rf'{len(runaways_all)} Runaway(s)'
    ax1.scatter(_runaways_bp_rp,_runaways_gmag,s=8,alpha=0.5, zorder=3,
                c=runaways_all['Temp. Est'],
                cmap='spring_r',norm=plt.Normalize(4000, 23000),
                label=label_all_run)
    # main runaways T > 10,000K
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
        if isinstance(event.artist,matplotlib.lines.Line2D):
            return
        global annotations, table2
        for ann in annotations:
            ann.remove()
        annotations = []
        
        ind = event.ind
        for index in ind:
            index = ind[0]
            temp_est = runaways['Temp. Est'][index]
            ann = ax1.annotate(f'{temp_est:,.0f} K', (runaways_bp_rp[index], runaways_gmag[index]), xytext=(0, -25), textcoords='offset points', color='black', ha='center', fontweight='bold')
            ann.set_path_effects([pe.withStroke(linewidth=4, foreground='white')])
            annotations.append(ann)
            
            # Add the table below the annotation
            table_data = [
                ['ID', f'{runaways[index]["Source"]}'],
                ['Temp. Est.', f'{runaways[index]["Temp. Est"]:,.0f} K'],
                ['Coord.', f'{runaways["RA_ICRS_1"][index]:.4f} {"+" if runaways["DE_ICRS_1"][index] >= 0 else ""}{runaways["DE_ICRS_1"][index]:.4f}'],
                ['Gmag', f'{runaways["Gmag"][index]}'],
                ['Dist.', f'{runaways["rgeo"][index]:.0f}'+'$\pm$'+f'{(runaways["rgeo"][index]-runaways["b_rgeo"][index]):.0f}']

                # Add more rows with actual values
            ]
            table_bbox = [0.50, 0.0, 0.50, 0.2]  # [left, bottom, width, height]
            table2 = ax1.table(cellText=table_data, cellLoc='right', loc='lower center', bbox=table_bbox)
            # table2.auto_set_column_width(col=[0,1])

            for key, cell in table2._cells.items():
                cell.set_linewidth(0.5)  # Set the border width
                cell.set_edgecolor('lightgray')  # Set the border color
            
        plt.draw()
    errorbar_runaways = ax1.errorbar(runaways['BP-RP'],runaways['Gmag'],linestyle="None",xerr=runaways['e_BP-RP']+0.02,yerr=runaways['e_Gmag'],color='red')
    scatter_main = ax1.scatter(runaways_bp_rp,runaways_gmag,picker=True,s=30, 
                            c=runaways['Temp. Est'],cmap='spring_r',zorder=4,norm=plt.Normalize(4000, 23000),
                            label=rf'{len(runaways)} Runaway(s) with T > {config["runaway_temp"]:,} K')
    
    
    # Create interactive legend    
    map_legend_to_ax = {} 
    pickradius = 5

    leg = ax1.legend(loc='upper right', fancybox=True, shadow=True)
    for legend_line, ax_line in zip(leg.get_lines(), lines):
        legend_line.set_picker(pickradius)  # Enable picking on the legend line.
        map_legend_to_ax[legend_line] = ax_line
        if ax_line.get_alpha()==1:
            ax_line.set_visible(False)
            legend_line.set_alpha(0.1)

    def on_pick(event):
        # On the pick event, find the original line corresponding to the legend
        # proxy line, and toggle its visibility.
        legend_line = event.artist
        # Do nothing if the source of the event is not a legend line.
        if legend_line not in map_legend_to_ax:
            return

        ax_line = map_legend_to_ax[legend_line]
        visible = not ax_line.get_visible()
        ax_line.set_visible(visible)
        # Change the alpha on the line in the legend, so we can see what lines
        # have been toggled.
        legend_line.set_alpha(1.0 if visible else 0.2)
        fig.canvas.draw()


    # Connect the pick_event to the onpick3 function
    fig.canvas.mpl_connect('pick_event', onpick3)
    fig.canvas.mpl_connect('pick_event', on_pick)
    fig.canvas.mpl_connect('button_release_event', on_release)


    colorbar = fig.colorbar(scatter_main,ax=ax1)
    colorbar.set_label('Temperature (K)')
    # ax1.legend(loc='upper right')
    if save:
        plt.savefig(f'{cluster.name}/{cluster.name}_cmd.{save}')
    return ax1

def runCode(cluster,save='png',psr=False,separation_factor=2,dist_filter_factor=1,**kwargs):
    if isinstance(cluster,str):
        cluster = Cluster(cluster)

    elif isinstance(cluster,Cluster):
        cluster=cluster
    print(f'{cluster.name:=>50}'+f'{"":=<50}')
    cluster.generate_tables()
    theoretical_data = theoretical_isochrone(cluster,output="table",printing=False,**kwargs)
    print(f'{cluster.name+": traceback":->50}'+f'{"":-<50}')
    fs = cluster.read_table('fs')
    runaways_all = get_runaways(cluster,fs,theoretical_data,separation_factor,dist_filter_factor)
    # display(runaways_all)
    mask = [T > config['runaway_temp'] for T in runaways_all['Temp. Est']]
    runaways = runaways_all[mask]
    print(f'{cluster.name+": runaways and isochrone plot":->50}'+f'{"":-<50}')
    cluster = Cluster(cluster.name)
    plot_cmd(cluster,save=save,**kwargs)
    # print(kwargs)
    if psr:
        plot_traceback(cluster,save=save)
    plot_traceback_clean(cluster,save=save, separation_factor=separation_factor)
    print(f"{len(runaways)} Runaway(s) found",runaways)

def reRunCode(cluster,separation_factor=2,dist_filter_factor=1, **kwargs):
    if isinstance(cluster, str):
        cluster = Cluster(cluster)
    elif isinstance(cluster, Cluster):
        cluster = cluster
    cluster.clean()
    runCode(cluster,separation_factor=separation_factor,dist_filter_factor=dist_filter_factor, **kwargs)


def plot_traceback_clean(cluster,save=False, separation_factor=False):
    search_arcmin = cluster.calculate_search_arcmin()
    cluster.plot_search_region(display=None,pixels='800')
    fits_path = f"{cluster.name}/{cluster.name}_extra{config['Cluster']['search_extent']}pc.fits"
    fits_file = fits.open(fits_path)
    image = fits_file[0]
    wcs = WCS(image.header)

    fig, ax2 = plt.subplots(subplot_kw={'projection': wcs}, figsize=(10, 8.8))
    ax2.imshow(image.data, cmap='gray')
    # Add a circle representing the cluster radius
    if separation_factor!=2 and separation_factor!=False:
        test_circle = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), #center
                                 radius=(image.shape[0]/2)*cluster.diameter/(separation_factor*search_arcmin), #radius
                                 edgecolor='red', facecolor='none', ls='dashed', label=f'Test Diameter = {cluster.diameter*2/separation_factor}',linewidth=1, alpha=0.6)
        
        ax2.add_artist(test_circle)

    # Add a circle representing the cluster radius
    circle_cluster = plt.Circle((image.data.shape[1] / 2, image.data.shape[0] / 2), radius=(image.shape[0]/2)*cluster.diameter/(2*search_arcmin),
                    edgecolor='red', facecolor='none', ls='dashed', label=f'Cluster Diameter (r50) = {cluster.diameter}',linewidth=1)
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
    ax2.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontweight='bold')

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
    fits_file.close()

    if save:
        plt.savefig(f'{cluster.name}/{cluster.name}_traceback.{save}')
    
    return ax2,wcs


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
                    edgecolor='red', facecolor='none', ls='dashed', label=f'Cluster Diameter (r50) = {cluster.diameter}',linewidth=1)
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
    # scatter_my_members, = ax2.plot(my_members_pix_coords[0],my_members_pix_coords[1],'.',markersize=5, color='blue', label=f'{len(my_members)} My Members')



    # Read runaways
    runaways_all = cluster.read_table('runaways') #changing this to runaways_all plots all the runaways, not just the >10000K ones
    runaways_all['Source'] = runaways_all['Source'].astype(str)
    runaways_all['Gmag'] = (runaways_all['Gmag'].round(2)).astype(str)
    # runaways_all['RAVE5'] = runaways_all['RAVE5'].astype(str)

    psrs = search_psr(cluster.name)
    n_psr = len(psrs)

    display(psrs)
    if len(runaways_all)<0:
        return None

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
    # _ = str(round(_.value,2))+"'"

    __ = cluster.distance.value/1000
    add_scalebar(ax2, _, color="yellow", label=f"{_.value}' or {scalebar_length.value}pc (at {__:.2f}kpc)", size_vertical=0.5)
    x_min, x_max = ax2.get_xlim()
    y_min, y_max = ax2.get_ylim()
    # ax2.annotate(f"{_:.2f}", xy=(0.83*x_max,0.08*y_max), color='yellow', ha='center',fontweight='bold')

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
                vs_psr = [340]*u.km/u.s #deleteme
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
    # scatter_stars_in_region_legend, scatter_dias_members_legend,scatter_my_members_legend = legend.get_lines()
    scatter_stars_in_region_legend, scatter_dias_members_legend = legend.get_lines()
    scatter_stars_in_region_legend.set_picker(True)
    scatter_stars_in_region_legend.set_pickradius(10)
    scatter_dias_members_legend.set_picker(True)
    scatter_dias_members_legend.set_pickradius(10)
    # scatter_my_members_legend.set_picker(True)
    # scatter_my_members_legend.set_pickradius(10)s

    graphs = {}
    graphs[scatter_stars_in_region_legend] = scatter_stars_in_region
    graphs[scatter_dias_members_legend] = scatter_dias_members
    # graphs[scatter_my_members_legend] = scatter_my_members


    # colorbar = fig.colorbar(scatter_runaways_all,ax=ax2)
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
    
def plot_pm(cluster):
    fig, ax3 = plt.subplots(figsize=(8, 8))
    #Plot 3: proper motions
    ax3.grid('True')
    #Scatter stars in the region
    sir = cluster.stars_in_region()
    sir_pmRA,sir_pmDE = sir['pmRA'],sir['pmDE']
    ax3.scatter(sir_pmRA,sir_pmDE,s=2, color='grey',label=f'{len(sir)} stars in the region {cluster.calculate_search_arcmin()}')
    # Scatter dias members
    dias_members = cluster.members
    dias_pmRA,dias_pmDE = dias_members['pmRA'],dias_members['pmDE']
    ax3.scatter(dias_pmRA,dias_pmDE,s=15, color='black',label=rf'{len(dias_members)} Dias Members P $\geq$ {config["memb_prob"]}, Q $\geq$ {config["parallax_quality_threshold"]}')

    #Scatter runaways
    runaways_all = cluster.read_table('runaways_all')
    runaways = cluster.read_table('runaways')
    try:
        label_all_run = rf'{len(runaways_all)} Runaway(s) $T_{{max}}$ = {max(runaways_all["Temp. Est"]):,.0f} K'
    except:
        label_all_run = rf'{len(runaways_all)} Runaway(s)'
    ax3.scatter(runaways_all['pmRA'],runaways_all['pmDE'],s=10,alpha=0.7, 
            c=runaways_all['Temp. Est'],
            cmap='spring_r',norm=plt.Normalize(4000, 23000),
            label=label_all_run)
    scatter_main = ax3.scatter(runaways['pmRA'],runaways['pmDE'],picker=True,s=40, 
                            c=runaways['Temp. Est'],cmap='spring_r',norm=plt.Normalize(4000, 23000),
                            label=rf'{len(runaways)} Runaway(s) with T > {config["runaway_temp"]:,} K')
    
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
        if isinstance(event.artist,matplotlib.lines.Line2D):
            return
        global annotations, table2
        for ann in annotations:
            ann.remove()
        annotations = []
        
        ind = event.ind
        for index in ind:
            index = ind[0]
            temp_est = runaways['Temp. Est'][index]
            ann = ax3.annotate(f'{temp_est:,.0f} K', (runaways['pmRA'][index], runaways['pmDE'][index]), xytext=(0, -25), textcoords='offset points', color='black', ha='center', fontweight='bold')
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
            table2 = ax3.table(cellText=table_data, cellLoc='right', loc='lower center', bbox=table_bbox)
            table2.auto_set_column_width(col=[0,1])
            for key, cell in table2._cells.items():
                cell.set_linewidth(0.5)  # Set the border width
                cell.set_edgecolor('lightgray')  # Set the border color
            
        plt.draw()
    
    fig.canvas.mpl_connect('pick_event', onpick3)
    fig.canvas.mpl_connect('button_release_event', on_release)

    
    
    ax3.set_xlabel(r"$\mu_{\alpha}^*$ (mas/yr)")
    ax3.set_ylabel(r"$\mu_{\delta}$ (mas/yr)")
    ax3.set_title("Proper Motions")
    ax3.legend(loc="upper left")

def search_psr(cluster, extra=config['psr_extra'], radial_tolerance=config['psr_radial_tolerance'], save=True):
    # Check if cluster is a string
    if isinstance(cluster, str):
        cluster = Cluster(cluster)
        print(f'{"Searching pulsars near":->50}' + f' {cluster.name:-<50} RA:{cluster.coordinates.ra} DE:{cluster.coordinates.dec}')
    
    # Check if cluster is an instance of Cluster
    elif isinstance(cluster, Cluster):
        cluster = cluster
        print(f'{"Searching pulsars near":->50}' + f' {cluster.name:-<50} RA:{cluster_coord.ra} DE:{cluster_coord.dec}')
    
    # Check if cluster is an instance of SkyCoord
    elif isinstance(cluster, SkyCoord):
        cluster_coord = cluster
        print(f'{"Searching pulsars near coordinates":->50}RA:{cluster_coord.ra} DE:{cluster_coord.dec}')
    
    # Invalid type for cluster
    else:
        raise TypeError("cluster must be a string, Cluster, or SkyCoord object")

    # Handle case when cluster is a Cluster object
    if isinstance(cluster, Cluster):
        cluster_coord = cluster.coordinates
        print(f'({extra} pc from cluster edge, {radial_tolerance * 100:.0f}% radial tolerance around distance {cluster.distance.to(u.kpc)})')

    cluster_coord_ra_str = cluster_coord.ra.to_string(unit='hourangle', sep=':', precision=3, pad=True)
    cluster_coord_dec_str = cluster_coord.dec.to_string(unit='degree', sep=':', precision=3, pad=True)

    psr_search_deg = cluster.calculate_search_arcmin(extra=extra).to(u.deg).value if isinstance(cluster, Cluster) else (extra / 60.0)
    c = [cluster_coord_ra_str, cluster_coord_dec_str, psr_search_deg]
    query = QueryATNF(params=['NAME', 'JNAME', 'RAJD', 'DECJD', 'DIST', 'DIST_DM', 'AGE', 'PMRA', 'PMDEC', 'S400', 'ASSOC', 'AGE_I', 'PX', 'P0', 'P1'], circular_boundary=c)
    print(f'{len(query.table)} Pulsar(s) within {psr_search_deg:.2f} degrees ({extra} pc)')

    r_close = (1 - radial_tolerance) * cluster.distance.to(u.kpc).value if isinstance(cluster.distance, astropy.coordinates.distances.Distance) else 0
    r_far = (1 + radial_tolerance) * cluster.distance.to(u.kpc).value if isinstance(cluster.distance, astropy.coordinates.distances.Distance) else float('inf')

    dist_filtered_psr_table = query.table[(query.table['DIST'] > r_close) & (query.table['DIST'] < r_far)]
    psr_coords = SkyCoord(ra=dist_filtered_psr_table['RAJD'], dec=dist_filtered_psr_table['DECJD'], pm_ra_cosdec=dist_filtered_psr_table['PMRA'], pm_dec=dist_filtered_psr_table['PMDEC'])
    psr_sep = psr_coords.separation(cluster_coord).to(u.arcmin)
    dist_filtered_psr_table.add_column(psr_sep, name='Separation')
    dist_filtered_psr_table.sort('Separation')
    dist_filtered_psr_table['AGE'] = dist_filtered_psr_table['AGE'].to(u.kyr)
    dist_filtered_psr_table['AGE_I'] = dist_filtered_psr_table['AGE_I'].to(u.kyr)
    print(f'Of which {len(dist_filtered_psr_table)} Pulsar(s) from {r_close:.2f} kpc to {r_far:.2f} kpc')
    
    maintable = dist_filtered_psr_table['Separation', 'NAME', 'JNAME', 'RAJD', 'DECJD', 'DIST', 'DIST_DM', 'AGE', 'PMRA', 'PMDEC', 'S400', 'ASSOC', 'AGE_I', 'PX']
    
    if save:
        dir_name = cluster.name if isinstance(cluster, Cluster) else 'SkyCoord_Search'
        os.makedirs(dir_name, exist_ok=True)
        maintable.write(f'{dir_name}/{dir_name}_psrs.tsv', format='ascii.ecsv', overwrite=True)
        maintable.to_pandas().to_excel(os.path.join(dir_name, f'{dir_name}_psrs.xlsx'), index=False)
    
    return maintable
    #add distance from cluster column
def nearest_cluster(objectname, output=False):
    result_table = Simbad.query_object(objectname)
    cluster_table = cluster_list['Cluster','RA_ICRS','DE_ICRS']
    ra = result_table['RA'].value.data[0]
    dec = result_table['DEC'].value.data[0]

    # Create a SkyCoord object
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))

    # Input coordinate
    input_coord = coord

    # Convert the input coordinate to a SkyCoord object
    target_coord = coord

    # Calculate the angular separation for each cluster in the table
    separations = []
    for row in cluster_table:
        cluster_coord = SkyCoord(ra=row['RA_ICRS'], dec=row['DE_ICRS'], distance=row['Dist'], unit=(u.deg, u.deg), frame='icrs')
        separation = target_coord.separation(cluster_coord).deg
        separations.append(separation)

    # Find the index of the minimum separation
    min_index = separations.index(min(separations))

    # Get the nearest cluster
    nearest_cluster = cluster_table[min_index]['Cluster']
    if output:
        print(f"The nearest cluster to the {objectname} is: {nearest_cluster}"+f' at a separation of {min(separations)*u.deg.to(u.arcmin):.2f} arcmin')
    return nearest_cluster

def what_line_is_it(line):
    illss = Table.read('spectra/illss_stellar.tsv', format='ascii.ecsv')
    wavelengths = illss['lambda']

    closest_index = np.abs(wavelengths - line).argmin()
    closest_wa = illss[closest_index]  # Corrected line
    return closest_wa

def observe_from_gsh(clustername,obstime=None):
    import astropy
    print('='*50)
    clusters = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
    mask = (clusters['Cluster']==f'{clustername}')
    plt.style.use(astropy_mpl_style)
    quantity_support()
    # Define the Earth location
    gsh = EarthLocation(lon=11.482463*u.deg, lat=50.928449*u.deg, height=370*u.m)
    utcoffset = 1*u.hour  # Eastern Daylight Time
    # Get the current time
    if obstime is not None:
        time = Time(f'{obstime}')- utcoffset
        print('Observation for given time (utc):', time)
    elif obstime is None:
        time = Time.now()
        print('Observation for current time (utc):', time)
    # Define M33 coordinates
    if type(clustername) == astropy.coordinates.sky_coordinate.SkyCoord:
        print('given coordinates')
        m33 = clustername
        label = (round(coord.ra.value,2),round(coord.dec.value,2))
    elif type(clustername) == str:
        print('given cluster name')
        m33 = SkyCoord(clusters[mask][0]['RA_ICRS']*u.deg,clusters[mask][0]['DE_ICRS']*u.deg,frame='icrs')
        label = f'{clustername}'
    print('Coordinates in deg:',m33.to_string('decimal'))
    print('Coordinates in hms:',m33.to_string('hmsdms'))
    # Transform M33 to AltAz coordinates
    m33altaz = m33.transform_to(AltAz(obstime=time, location=gsh))
    print(f"{clustername}'s Altitude = {m33altaz.alt:.2}")
    # Plot the first subplot (Airmass vs. Hours from EDT Midnight)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    # Plot Airmass vs. Hours from EDT Midnight
    midnight = Time(str((time+1*u.day).datetime.date()))
    print('Midnight:', midnight)
    midnight = midnight - utcoffset
    delta_midnight = np.linspace(-2, 10, 100)*u.hour
    frame_July13night = AltAz(obstime=midnight+delta_midnight, location=gsh)
    m33altazs_July13night = m33.transform_to(frame_July13night)
    m33airmasss_July13night = m33altazs_July13night.secz
    axs[0].plot(delta_midnight, m33airmasss_July13night)
    axs[0].set_xlim(-2, 10)
    axs[0].set_ylim(1, 4)
    axs[0].set_xlabel('Hours from gsh Midnight')
    axs[0].set_ylabel('Airmass [Sec(z)]')
    axs[0].set_title('Airmass vs. Hours from Midnight')
    # Plot the second subplot (Altitude vs. Hours from EDT Midnight)
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=gsh)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    moon_July12_to_13 = get_body("moon", times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)
    m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)
    observable = (max(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -12*u.deg)]) > 30*u.deg)
    print('Max elevation:',max(m33altazs_July12_to_13.alt))
    print('Observable:',observable)

    axs[1].plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    axs[1].plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
    sc = axs[1].scatter(delta_midnight, m33altazs_July12_to_13.alt,
                    c=m33altazs_July12_to_13.az.value, label=label, lw=0, s=8,
                    cmap='twilight')
    axs[1].fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0) #sunset
    axs[1].fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs_July12_to_13.alt < -18*u.deg, color='k', zorder=0) #astronomical twilight
    cbar = plt.colorbar(sc, ax=axs[1], label='Azimuth [deg]')
    axs[1].legend(loc='upper left')
    axs[1].set_xlim(-12*u.hour, 12*u.hour)
    axs[1].set_xticks((np.arange(13)*2-12)*u.hour)
    axs[1].set_ylim(0*u.deg, 90*u.deg)
    axs[1].set_xlabel('Hours from gsh Midnight')
    axs[1].set_ylabel('Altitude [deg]')
    axs[1].set_title('Altitude vs. Hours from Midnight')
    # Adjust layout
    plt.tight_layout()
    # Show the plot
    plt.show()
    # print(sunaltazs_July12_to_13.alt < -0*u.deg,len(sunaltazs_July12_to_13.alt < -0*u.deg))
    # print(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)],len(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    # print(delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)],len(delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    # print((delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)])[(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)])>30*u.deg],len((delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)])[(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)])>30*u.deg]))
    sunset = (sunaltazs_July12_to_13.alt < -0*u.deg)
    astro_twilight = (sunaltazs_July12_to_13.alt < -18*u.deg)

    # print(len(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg])),len(astro_twilight))
    try:
        print('Total Observation time (elevation >30 and sunset):',(delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg][-1]-(delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg][0])
    except:
        print('Not Observable during sunset')
    try:    
        print('Good Observation time (elevation >30 and astronomical twilight):',(delta_midnight[astro_twilight])[(m33altazs_July12_to_13.alt[astro_twilight])>30*u.deg][-1]-(delta_midnight[astro_twilight])[(m33altazs_July12_to_13.alt[astro_twilight])>30*u.deg][0])
    except:
        print('Not Observable during astro twilight')

    # print('Good Observation time (elevation >30 and astronomical twilight):',(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg]) and ((delta_midnight[sunset])[astro_twilight]))[-1]-(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg]) and ((delta_midnight[sunset])[astro_twilight]))[0])
    # print(max(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    return observable



def observe_from_gsh_clean(clustername,obstime=None):
    import astropy
    # print('='*50)
    clusters = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
    mask = (clusters['Cluster']==f'{clustername}')
    # plt.style.use(astropy_mpl_style)
    quantity_support()
    # Define the Earth location
    gsh = EarthLocation(lon=11.482463*u.deg, lat=50.928449*u.deg, height=370*u.m)
    utcoffset = 1*u.hour  # Eastern Daylight Time
    # Get the current time
    if obstime is not None:
        time = Time(f'{obstime}')- utcoffset
        # print('Observation for given time (utc):', time)
    elif obstime is None:
        time = Time.now()
        # print('Observation for current time (utc):', time)
    # Define M33 coordinates
    if type(clustername) == astropy.coordinates.sky_coordinate.SkyCoord:
        # print('given coordinates')
        m33 = clustername
        label = (round(coord.ra.value,2),round(coord.dec.value,2))
    elif type(clustername) == str:
        # print('given cluster name')
        m33 = SkyCoord(clusters[mask][0]['RA_ICRS']*u.deg,clusters[mask][0]['DE_ICRS']*u.deg,frame='icrs')
        label = f'{clustername}'
    # print('Coordinates in deg:',m33.to_string('decimal'))
    # print('Coordinates in hms:',m33.to_string('hmsdms'))
    # Transform M33 to AltAz coordinates
    m33altaz = m33.transform_to(AltAz(obstime=time, location=gsh))
    print(f"{clustername}'s Altitude = {m33altaz.alt:.2}")
    # Plot the first subplot (Airmass vs. Hours from EDT Midnight)
    # fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    # Plot Airmass vs. Hours from EDT Midnight
    midnight = Time(str((time+1*u.day).datetime.date()))
    # print('Midnight:', midnight)
    midnight = midnight - utcoffset
    delta_midnight = np.linspace(-2, 10, 100)*u.hour
    frame_July13night = AltAz(obstime=midnight+delta_midnight, location=gsh)
    m33altazs_July13night = m33.transform_to(frame_July13night)
    m33airmasss_July13night = m33altazs_July13night.secz
    # axs[0].plot(delta_midnight, m33airmasss_July13night)
    # axs[0].set_xlim(-2, 10)
    # axs[0].set_ylim(1, 4)
    # axs[0].set_xlabel('Hours from gsh Midnight')
    # axs[0].set_ylabel('Airmass [Sec(z)]')
    # axs[0].set_title('Airmass vs. Hours from Midnight')
    # Plot the second subplot (Altitude vs. Hours from EDT Midnight)
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=gsh)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    moon_July12_to_13 = get_body("moon", times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)
    m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)
    observable = (max(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -12*u.deg)]) > 30*u.deg)
    # print('Max elevation:',max(m33altazs_July12_to_13.alt))
    print('Observable:',observable)

    # axs[1].plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    # axs[1].plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
    # sc = axs[1].scatter(delta_midnight, m33altazs_July12_to_13.alt,
    #                 c=m33altazs_July12_to_13.az.value, label=label, lw=0, s=8,
    #                 cmap='twilight')
    # axs[1].fill_between(delta_midnight, 0*u.deg, 90*u.deg,
    #                 sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0) #sunset
    # axs[1].fill_between(delta_midnight, 0*u.deg, 90*u.deg,
    #                 sunaltazs_July12_to_13.alt < -18*u.deg, color='k', zorder=0) #astronomical twilight
    # cbar = plt.colorbar(sc, ax=axs[1], label='Azimuth [deg]')
    # axs[1].legend(loc='upper left')
    # axs[1].set_xlim(-12*u.hour, 12*u.hour)
    # axs[1].set_xticks((np.arange(13)*2-12)*u.hour)
    # axs[1].set_ylim(0*u.deg, 90*u.deg)
    # axs[1].set_xlabel('Hours from gsh Midnight')
    # axs[1].set_ylabel('Altitude [deg]')
    # axs[1].set_title('Altitude vs. Hours from Midnight')
    # # Adjust layout
    # plt.tight_layout()
    # # Show the plot
    # plt.show()
    # print(sunaltazs_July12_to_13.alt < -0*u.deg,len(sunaltazs_July12_to_13.alt < -0*u.deg))
    # print(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)],len(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    # print(delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)],len(delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    # print((delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)])[(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)])>30*u.deg],len((delta_midnight[(sunaltazs_July12_to_13.alt < -0*u.deg)])[(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)])>30*u.deg]))
    sunset = (sunaltazs_July12_to_13.alt < -0*u.deg)
    astro_twilight = (sunaltazs_July12_to_13.alt < -18*u.deg)

    # # print(len(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg])),len(astro_twilight))
    # try:
    #     print('Total Observation time (elevation >30 and sunset):',(delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg][-1]-(delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg][0])
    # except:
    #     print('Not Observable during sunset')
    # try:    
    #     print('Good Observation time (elevation >30 and astronomical twilight):',(delta_midnight[astro_twilight])[(m33altazs_July12_to_13.alt[astro_twilight])>30*u.deg][-1]-(delta_midnight[astro_twilight])[(m33altazs_July12_to_13.alt[astro_twilight])>30*u.deg][0])
    # except:
    #     print('Not Observable during astro twilight')

    # print('Good Observation time (elevation >30 and astronomical twilight):',(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg]) and ((delta_midnight[sunset])[astro_twilight]))[-1]-(((delta_midnight[sunset])[(m33altazs_July12_to_13.alt[sunset])>30*u.deg]) and ((delta_midnight[sunset])[astro_twilight]))[0])
    # print(max(m33altazs_July12_to_13.alt[(sunaltazs_July12_to_13.alt < -0*u.deg)]))
    return observable
#plot the find_members and dias_members on the search region fits
#find out the reason why some clusters were not searched. probably because of the search distance
#3D visualisation
#unique config files for each cluster. add string on top of it indicating last change date
