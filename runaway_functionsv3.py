import os
import astropy
import numpy as np
from astropy import units as u
from astroquery.vizier import Vizier
from astropy.table import Table, join,vstack
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
import warnings
from astropy.utils.metadata import MergeConflictWarning
from astropy.utils.exceptions import ErfaWarning

import time
# import logging
from typing import List
import yaml
from astropy.stats import sigma_clip
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from matplotlib import pyplot as plt
from astroquery.skyview import SkyView
from regions import CircleSkyRegion, PointSkyRegion, LineSkyRegion
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.simbad import Simbad
import matplotlib.patches as patches

dias2021 = Table.read("dias2021.tsv", format="ascii.ecsv")
maskplx = dias2021['Plx'] > 0.3
maskage = dias2021['logage'] < 7.7
workclusters = []
for clustername in dias2021[maskplx & maskage][:]['Cluster']:
    if clustername not in ['ASCC_79','BH_164','BH_23','Collinder_135','Collinder_140','Gulliver_9','IC_2391','IC_2602','Mamajek_1','Platais_8','UPK_535','UPK_606','UPK_640','Berkeley_59','COIN-Gaia_37','Ivanov_4','LP_1937','Sigma_Ori','UBC_632']:
        workclusters.append(clustername)

def read_yaml_file(file_path):
    '''
    Read the configuration file for the rest of the code. 
    This contains the various parameters for the code to run.
    '''
    with open(file_path, 'r') as yaml_file:
        config = yaml.safe_load(yaml_file)
    return config
config = read_yaml_file('config2.yaml')


class ClusterDias:
    def __init__(self, name:str) -> None:
        #download the catalog if not already in the cwds
        if os.path.exists('./dias2021.tsv'):
            dias2021 = Table.read('dias2021.tsv', format='ascii.ecsv')
        else:
            Vizier.ROW_LIMIT = -1
            dias2021 = Vizier.get_catalogs(catalog="J/MNRAS/504/356/table12")[0]
            mask_logage = dias2021['logage'] < 7.7
            mask_plx = dias2021['Plx'] > 0.3
            dias2021 = dias2021[mask_logage & mask_plx]
        #select the row from Dias catalog
        cluster_row = dias2021[dias2021['Cluster'] == name][0]
        #cluster parameters from the row
        self.name = name
        self.all = cluster_row
        self.r50 = cluster_row['r50']*u.deg
        self.N = cluster_row['N']
        self.skycoord = SkyCoord(ra=cluster_row['RA_ICRS']*u.deg,
                                 dec=cluster_row['DE_ICRS']*u.deg,
                                 distance=cluster_row['Dist']*u.pc,
                                 pm_ra_cosdec=cluster_row['pmRA']*u.mas/u.yr,
                                 pm_dec=cluster_row['pmDE']*u.mas/u.yr,
                                 obstime=(Time('J2000')+1*u.Myr))
        self.ra = self.skycoord.ra
        self.dec = self.skycoord.dec
        self.distance,self.e_distance = self.skycoord.distance,cluster_row['e_Dist']*u.pc
        self.pm_ra_cosdec, self.e_pm_ra_cosdec = self.skycoord.pm_ra_cosdec,cluster_row['e_pmRA']*u.mas/u.yr
        self.pm_dec, self.e_pm_dec = self.skycoord.pm_dec, cluster_row['e_pmDE']*u.mas/u.yr
        self.Av, self.e_Av = cluster_row['Av']*u.mag,cluster_row['e_Av']*u.mag
        self.logage, self.e_logage = cluster_row['logage'],cluster_row['e_logage']
        self.FeH, self.e_FeH = cluster_row['__Fe_H_'],cluster_row['e__Fe_H_']
        self.RV, self.e_RV = cluster_row['RV'],cluster_row['e_RV']
        self.NRV = cluster_row['NRV']
        self.mymembers = self.members()
        

    def members(self) -> None:
        members = (Table.read(f'./Clusters_Dias/{self.name}.dat', format='ascii.tab'))[2:] #first two rows removed
        mask_BPRP_exists = members['BP-RP'] != " ---   "
        members = members[mask_BPRP_exists]
        members['Source'] = members['Source'].astype(np.int64)
        members['Pmemb'] = members['Pmemb'].astype(float)
        members['Plx'] = members['Plx'].astype(float)*u.mas
        members['e_Plx'] = members['e_Plx'].astype(float)*u.mas
        members['RAdeg'] = members['RAdeg'].astype(float)*u.deg
        members['DEdeg'] = members['DEdeg'].astype(float)*u.deg
        members['pmRA'] = members['pmRA'].astype(float)*u.mas/u.yr
        members['pmDE'] = members['pmDE'].astype(float)*u.mas/u.yr
        members['e_pmRA'] = members['e_pmRA'].astype(float)*u.mas/u.yr
        members['e_pmDE'] = members['e_pmDE'].astype(float)*u.mas/u.yr
        members['Gmag'] = members['Gmag'].astype(float)*u.mag
        members['e_Gmag'] = members['e_Gmag'].astype(float)*u.mag
        members['BPmag'] = members['BPmag'].astype(float)*u.mag
        members['e_BPmag'] = members['e_BPmag'].astype(float)*u.mag
        members['RPmag'] = members['RPmag'].astype(float)*u.mag
        members['e_RPmag'] = members['e_RPmag'].astype(float)*u.mag
        members['BP-RP'] = members['BP-RP'].astype(float)*u.mag
        members['e_BP-RP'] = members['e_BPmag']+members['e_RPmag']
        members = members["RAdeg","DEdeg", "Source","Pmemb","Plx","e_Plx","pmRA","e_pmRA","pmDE","e_pmDE","Gmag","e_Gmag","BP-RP","e_BP-RP"]
        return members

    def get_full_catalog() -> Table:
        dias2021 = Vizier.get_catalogs(catalog="J/MNRAS/504/356/table12")[0]
        dias2021.write('dias2021.tsv',format='ascii.ecsv', overwrite=True)
        return dias2021
    
    def theoretical_isochrone(self, params=None, returnparams=False, parsec_version=2):
        # self.members()
        params = params or {}
        Av = float(params.get('Av', None)) if params.get('Av') is not None else round(float(self.Av.value), 1)
        logage = float(params.get('logage', None)) if params.get('logage') is not None else round(float(self.logage), 1)
        FeH = float(params.get('FeH', None)) if params.get('FeH') is not None else round(float(self.FeH), 1)
        
        theo_iso_path = f"./Clusters/{self.name}/{self.name}_compare_data_out_Av{str(Av)}_logage{str(logage)}_FeH{str(FeH)}.isochrone{parsec_version}"
        # print(Av, logage, FeH)
        if os.path.exists(theo_iso_path):
            theo_iso = Table.read(theo_iso_path, format="ascii")
        else:
            theo_iso = get_theoretical_isochrone(Av=Av, logage=logage, FeH=FeH, parsec_version=parsec_version)
            theo_iso['Gmag'] = theo_iso['Gmag'] + 5 * np.log10(self.distance.value) - 5
            theo_iso['G_BP'] = theo_iso["G_BP"]+ 5 * np.log10(self.distance.value) - 5
            theo_iso['G_RP'] = theo_iso["G_RP"]+ 5 * np.log10(self.distance.value) - 5
            if parsec_version==1.2:
                theo_iso['Teff0'] = 10**theo_iso['logTe']
            theo_iso = theo_iso["Mass", "Teff0", "BP-RP", "Gmag", "G_BP", "G_RP", "logg", "logAge", "logL", "logTe", "Mini"]
            theo_iso.write(theo_iso_path, format="ascii", overwrite=True)
        if returnparams:
            return theo_iso, (Av, logage, FeH)
        else:
            return theo_iso

class ClusterCG:
    def __init__(self, name:str) -> None:
        self.name = name
        #download the catalog if not already in the cwds
        if os.path.exists('./CG2020.tsv'):
            CG2020 = Table.read('CG2020.tsv', format='ascii.ecsv')
        else:
            Vizier.ROW_LIMIT = -1
            CG2020 = Vizier.get_catalogs(catalog="J/A+A/633/A99/table1")[0]
        
        cluster_row = CG2020[CG2020['Cluster'] == name][0]
        self.r50 = cluster_row['r50']*u.deg

    def get_full_catalog() -> Table:
        CG2020 = Vizier.get_catalogs(catalog="J/A+A/633/A99/table1")[0]
        CG2020.write('CG2020.tsv',format='ascii.ecsv')
        return CG2020
        
class Cluster:
    def __init__(self, name:str) -> None:
        if not os.path.exists(f'Clusters/{name}'):
            os.mkdir(f'Clusters/{name}')
        clusterDias = ClusterDias(name=name)
        try:
            clusterCG = ClusterCG(name=name)
            self.r50 = clusterCG.r50
        except:
            self.r50 = clusterDias.r50
        self.name = name
        self.ra = clusterDias.skycoord.ra
        self.dec = clusterDias.skycoord.dec
        self.distance = clusterDias.skycoord.distance #don't delete, need for calculating and searching stars in the region
        self.rPhy = np.tan(self.r50) * self.distance
        self.search_arcmin = search_arcmin(self.distance, self.r50)
        # self.r50_phy = np.tan(self.r50) * self.distance #to be changed DR3

        suro2024 = Table.read('suro2024.tsv', format='ascii.ecsv')
        cluster_row = suro2024[suro2024['Cluster'] == name][0]
        self.skycoord = SkyCoord(ra=cluster_row['RA_ICRS']*u.deg,
                                dec=cluster_row['DE_ICRS']*u.deg,
                                distance=cluster_row['Dist']*u.pc,
                                pm_ra_cosdec=cluster_row['pmRA']*u.mas/u.yr,
                                pm_dec=cluster_row['pmDE']*u.mas/u.yr,
                                obstime=(Time('J2000')+1*u.Myr))

        self.all = cluster_row
        self.distance,self.e_distance = cluster_row['Dist']*u.pc,cluster_row['e_Dist']*u.pc
        self.pm_ra_cosdec, self.e_pm_ra_cosdec = cluster_row['pmRA']*u.mas/u.yr,cluster_row['e_pmRA']*u.mas/u.yr
        self.pm_dec, self.e_pm_dec = cluster_row['pmDE']*u.mas/u.yr, cluster_row['e_pmDE']*u.mas/u.yr
        self.Av, self.e_Av = round(cluster_row['Av'],2)*u.mag,cluster_row['e_Av']*u.mag
        self.logage, self.e_logage = round(cluster_row['logage'],2),cluster_row['e_logage']
        self.FeH, self.e_FeH = round(cluster_row['__Fe_H_'],2),cluster_row['e__Fe_H_']
        self.RV, self.e_RV = cluster_row['RV'],cluster_row['e_RV']
        self.NRV = cluster_row['NRV']
        if os.path.exists(f"Clusters/{self.name}/{self.name}_members.tsv"):
            self.mymembers = Table.read(f"Clusters/{self.name}/{self.name}_members.tsv", format="ascii.ecsv")
        else:
            self.mymembers = self.members()
        self.kinematic_cluster = find_cluster(self.stars_in_region())
        #functions
        # self.members = (Table.read(f'./Clusters_Dias/{self.name}.dat', format='ascii.tab'))[2:] 
        # self.members_list = list(self.members['Source'].astype(np.int64))

    def members(self, pmemb:float = config['Cluster']['pmemb'], plxquality:float = config['Cluster']['plxquality'], add_members:List = []):
        sr = self.stars_in_region()
        diasmembers = (ClusterDias(name=self.name).members())['Source','Pmemb']
        #add the add_members in this with Pmemb=1
        configmembers = config.get('added_members',{}).get(self.name)
        print(f"{configmembers} found from config file") if configmembers is not None else None
        add_members = add_members + configmembers if configmembers is not None else []
        #add_members += configmembers if configmembers is not None else []
        #print(add_members)

        for mem_source in add_members:
            #print(1)
            diasmembers.add_row([mem_source,1])
        #display(diasmembers)
        members = join(sr, diasmembers, keys='Source', join_type='inner') #select the ones that dias says is a member
        print(f'{len(members)} out of {len(diasmembers)-len(add_members)} dias members found in search region')
        mask_pmemb = members['Pmemb'] >= pmemb
        mask_plxquality = members['Plx']/members['e_Plx'] >= plxquality
        members = members[mask_pmemb & mask_plxquality]
        members.sort("Gmag")
        members.write(f"Clusters/{self.name}/{self.name}_members.tsv", format="ascii.ecsv", overwrite=True)
        #update kinematic parameters
        self.changeParam(("N", len(members)))
        self.changeParam(("Plx", members['Plx'].mean()))
        self.changeParam(("e_Plx", members['Plx'].std()))
        self.changeParam(("Dist", members['rgeo'].mean()))
        self.changeParam(("e_Dist", members['rgeo'].std()))
        self.changeParam(("pmRA", members['pmRA'].mean()))
        self.changeParam(("e_pmRA", members['pmRA'].std()))
        self.changeParam(("pmDE", members['pmDE'].mean()))
        self.changeParam(("e_pmDE", members['pmRA'].std()))
        # finding RV from the RV mean of the cluster found
        kinematic_cluster = self.kinematic_cluster
        if (np.count_nonzero(~(kinematic_cluster['RV'].mask)))>5:
            self.changeParam(("RV", np.mean(kinematic_cluster['RV'])))
            self.changeParam(("e_RV", np.sqrt(np.sum(kinematic_cluster['e_RV'])**2/((np.count_nonzero(~(kinematic_cluster['RV'].mask))))**2)))
            self.changeParam(("NRV", np.count_nonzero(~(kinematic_cluster['RV'].mask))))
            return members
        elif (np.count_nonzero(~(members['RV'].mask)))>1:
            # self.restoreParam("RV")
            # self.restoreParam("e_RV")
            # self.restoreParam("NRV")
            self.changeParam(("RV", members['RV'].mean()))
            self.changeParam(("e_RV", np.sqrt(np.sum(members['e_RV'])**2/(np.count_nonzero(~(members['RV'].mask)))**2)))
            self.changeParam(("NRV", (np.count_nonzero(~(members['RV'].mask)))))
        else:
            self.restoreParam("RV")
            self.restoreParam("e_RV")
            self.restoreParam("NRV")
        return members


    
    def stars_in_region(self, star=None):
        stars_in_region_path =  f'Clusters/{self.name}/{self.name}_stars_in_region.tsv'
        
        if os.path.exists(stars_in_region_path):
            stars_in_region = Table.read(stars_in_region_path, format='ascii.ecsv')
            # stars_in_region['rmRA'] = stars_in_region['pmRA']-self.pm_ra_cosdec
            # stars_in_region['rmDE'] = stars_in_region['pmDE']-self.pm_dec
            # stars_in_region['e_rmRA'] = stars_in_region['e_pmRA']+self.e_pm_ra_cosdec
            # stars_in_region['e_rmDE'] = stars_in_region['e_pmDE']+self.e_pm_dec
            # stars_in_region['rRV'] = stars_in_region['RV']-self.RV
            # stars_in_region['e_rRV'] = stars_in_region['e_RV']+self.e_RV
            # #include members with high rRV as fast stars
            # stars_in_region['v_pec'] = 4.74*stars_in_region['rgeo'].value/1000*np.sqrt(((stars_in_region['rmRA'].value)**2+(stars_in_region['rmDE'].value)**2))*u.km/u.s
            # stars_in_region['v_pec3d'] = np.sqrt(stars_in_region['v_pec']**2+stars_in_region['rRV']**2)
            # #check
            # stars_in_region['e_vpec'] = 4.74 * stars_in_region['rgeo'].value/1000 * np.sqrt(((stars_in_region['rmRA'].value * stars_in_region['e_rmRA'].value)**2 + (stars_in_region['rmDE'].value * stars_in_region['e_rmDE'].value)**2)) * u.km/u.s
            # stars_in_region['e_vpec3d'] = np.sqrt((stars_in_region['v_pec'] / stars_in_region['v_pec3d'] * stars_in_region['e_vpec'])**2 + (stars_in_region['rRV'] / stars_in_region['v_pec3d'] * stars_in_region['e_rRV'])**2)

            # stars_in_region.write(stars_in_region_path, format='ascii.ecsv', overwrite=True)

        else:
            print("downloading sir")
            stars_in_region = self.get_stars_in_region()
            stars_in_region['rmRA'] = stars_in_region['pmRA']-self.pm_ra_cosdec
            stars_in_region['rmDE'] = stars_in_region['pmDE']-self.pm_dec
            stars_in_region['e_rmRA'] = stars_in_region['e_pmRA']+self.e_pm_ra_cosdec
            stars_in_region['e_rmDE'] = stars_in_region['e_pmDE']+self.e_pm_dec
            stars_in_region['rRV'] = stars_in_region['RV']-self.RV
            stars_in_region['e_rRV'] = stars_in_region['e_RV']+self.e_RV
            #include members with high rRV as fast stars
            stars_in_region['v_pec'] = 4.74*stars_in_region['rgeo'].value/1000*np.sqrt(((stars_in_region['rmRA'].value)**2+(stars_in_region['rmDE'].value)**2))*u.km/u.s
            stars_in_region['v_pec3d'] = np.sqrt(stars_in_region['v_pec']**2+stars_in_region['rRV']**2)
            #check
            stars_in_region['e_vpec'] = 4.74 * stars_in_region['rgeo'].value/1000 * np.sqrt(((stars_in_region['rmRA'].value * stars_in_region['e_rmRA'].value)**2 + (stars_in_region['rmDE'].value * stars_in_region['e_rmDE'].value)**2)) * u.km/u.s
            stars_in_region['e_vpec3d'] = np.sqrt((stars_in_region['v_pec'] / stars_in_region['v_pec3d'] * stars_in_region['e_vpec'])**2 + (stars_in_region['rRV'] / stars_in_region['v_pec3d'] * stars_in_region['e_rRV'])**2)

            stars_in_region.write(stars_in_region_path, format='ascii.ecsv')

        if star is None:
            return stars_in_region
        else:
            return stars_in_region[stars_in_region['Source']==star]
    
    def fast_stars_in_region(self):
        sir = self.stars_in_region()
        mask_fast2d = sir['v_pec'] > 17.6*u.km/u.s
        mask_vpec3ddosentexist = sir['v_pec3d'].mask
        mask_fast3d = sir['v_pec3d'] > 25*u.km/u.s
        
        return sir[(mask_fast2d & mask_vpec3ddosentexist) | mask_fast3d]
    def fs4giesler(self,outlocation=None):
        table = self.fast_stars_in_region()
        g = Table()
        g['TypeInput'] = np.ones_like(table['e_Plx'].value).astype(int)
        g['RA'] = table['SkyCoord'].ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True)
        g['DE'] = table['SkyCoord'].dec.to_string(unit='degree', sep=' ', precision=3, pad=True)
        g['Plx'] = table['rgeo'].to(u.mas, u.parallax())
        g['e_Plx'] = table['e_Plx']
        g['RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        # g['RV'] = table['RV']
        g['e_RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        # g['e_RV'] = table['e_RV']
        g['RVdist'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['pmRA'] = table['pmRA']
        g['e_pmRA'] = table['e_pmRA']
        g['pmDE'] = table['pmDE']
        g['e_pmDE'] = table['e_pmDE']
        g['Source'] = table['Source'].astype(str)
        # g['RV'] = g['RV'].filled(0) #fill masked values with 0
        # g['e_RV'] = g['e_RV'].filled(0) #fill masked values with 0
        #this is the first row for the cluster's motion
        new_row = [1,
           self.ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True), #ra
           self.dec.to_string(unit='degree', sep=' ', precision=3, pad=True), #dec
           ((self.all['Dist']*u.pc).to(u.mas, u.parallax())).value, #plx
           self.all['e_Plx'], #e_plx
        #    self.all['RV'], #RV
           0,   #RV
        #    self.all['e_RV'], #e_RV
           0,   #e_RV
           0, #RV_distribution
           self.all['pmRA'], #pmRA
           self.all['e_pmRA'], #e_pmRA
           self.all['pmDE'], #pmDE
           self.all['e_pmDE'], #e_pmDE
           self.name #ID
           ]
        new_table = Table(names=g.colnames, dtype=[col.dtype for col in g.columns.values()])
        new_table.add_row(new_row)
        g = vstack([new_table,g])
        if outlocation=='local':
            g.write(f'Clusters/{self.name}/{self.name}_fs4giesler.tsv', format='csv', delimiter='\t', overwrite=True)
        elif outlocation=='remote':
            g.write(f"/home/surodeep/suro_aiu/traceback/cluster_runaway/{self.name}/{self.name}_fs4giesler.tsv", format='csv', delimiter='\t', overwrite=True)
            

        return g
    def fs4giesler3d(self,outlocation=None):
        table = self.fast_stars_in_region()
        mask_vpec3dexists = ~table['v_pec3d'].mask
        table = table[mask_vpec3dexists]
        for star in config['observed_stars'][self.name]:
            table.add_row(self.stars_in_region(star)[0])
        g = Table()
        g['TypeInput'] = np.ones_like(table['e_Plx'].value).astype(int)
        g['RA'] = table['SkyCoord'].ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True)
        g['DE'] = table['SkyCoord'].dec.to_string(unit='degree', sep=' ', precision=3, pad=True)
        g['Plx'] = table['rgeo'].to(u.mas, u.parallax())
        g['e_Plx'] = table['e_Plx']
        # g['RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['RV'] = table['RV']
        # g['e_RV'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['e_RV'] = table['e_RV']
        g['RVdist'] = np.zeros_like(table['e_Plx'].value).astype(int)
        g['pmRA'] = table['pmRA']
        g['e_pmRA'] = table['e_pmRA']
        g['pmDE'] = table['pmDE']
        g['e_pmDE'] = table['e_pmDE']
        g['Source'] = table['Source'].astype(str)

        #this is the first row for the cluster's motion
        new_row = [1,
           self.ra.to_string(unit='hourangle', sep=' ', precision=3, pad=True), #ra
           self.dec.to_string(unit='degree', sep=' ', precision=3, pad=True), #dec
           ((self.all['Dist']*u.pc).to(u.mas, u.parallax())).value, #plx
           self.all['e_Plx'], #e_plx
           self.all['RV'], #RV
        #    0,   #RV
           self.all['e_RV'], #e_RV
        #    0,   #e_RV
           0, #RV_distribution
           self.all['pmRA'], #pmRA
           self.all['e_pmRA'], #e_pmRA
           self.all['pmDE'], #pmDE
           self.all['e_pmDE'], #e_pmDE
           self.name #ID
           ]
        new_table = Table(names=g.colnames, dtype=[col.dtype for col in g.columns.values()])
        new_table.add_row(new_row)
        g = vstack([new_table,g])
        if outlocation=='local':
            g.write(f'Clusters/{self.name}/{self.name}_fs4giesler.tsv', format='csv', delimiter='\t', overwrite=True)
        elif outlocation=='remote':
            g.write(f"/home/surodeep/suro_aiu/traceback/cluster_runaway3d/{self.name}/{self.name}_fs4giesler.tsv", format='csv', delimiter='\t', overwrite=True)
            
        return g        

    def get_stars_in_region(self) -> Table:
        c = ClusterDias(self.name).skycoord
        t1 = searchDR3(c,self.search_arcmin)
        t2 = searchDR3_dist(c,self.search_arcmin)
        t3 = merge_gaia_tables(t1,t2)
        t3.sort('Gmag')
        return t3
    
    def changeParam(self,change: tuple):
        param,new_value = change
        cluster_list = Table.read('dias2021.tsv', format='ascii.ecsv')
        cluster_list_mod = Table.read('suro2024.tsv', format='ascii.ecsv')
        _old_value  = cluster_list[cluster_list['Cluster'] == self.name][0][param]
        cluster_list_mod[param][cluster_list_mod['Cluster'] == self.name] = new_value
        cluster_list_mod.write('suro2024.tsv', format='ascii.ecsv', overwrite=True)
        warnings.simplefilter('ignore', FutureWarning)
        print(f'Changed {param:10} {_old_value:.2f} --> {new_value:.2f}')
    

    def restoreParam(self, param: str):
        # Read the original and modified tables
        cluster_list_original = Table.read('dias2021.tsv', format='ascii.ecsv')
        cluster_list_modified = Table.read('suro2024.tsv', format='ascii.ecsv')
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
        cluster_list_modified.write('suro2024.tsv', format='ascii.ecsv', overwrite=True)
    def prepare_trace(self):
        mainfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway/{self.name}"
        runawayfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway/{self.name}/runaways"
        trajfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway/{self.name}/runaway_trajectories"
        os.mkdir(mainfolder) if not os.path.exists(mainfolder) else None
        os.mkdir(runawayfolder) if not os.path.exists(runawayfolder) else None
        os.mkdir(trajfolder) if not os.path.exists(trajfolder) else None

        def update_config_template(config_file_path, output_file_path, **kwargs):
            # Read the template config file
            with open(config_file_path, 'r') as file:
                config_content = file.read()
            # Replace placeholders with actual values
            for key, value in kwargs.items():
                placeholder = f"val_{key}"
                config_content = config_content.replace(placeholder, str(value))
            # Write the updated content to the output file
            with open(output_file_path, 'w') as file:
                file.write(config_content)
            
        fs4giesler = self.fs4giesler(outlocation='remote')
        lastline = len(fs4giesler)+1
        # Example usage:
        #starno=5
        #template config
        config_file_path = "/home/surodeep/suro_aiu/traceback/cluster_runaway/template_config_trace.conf"
        #output the modified config
        output_file_path = f"/home/surodeep/suro_aiu/traceback/cluster_runaway/{self.name}/{self.name}_trace.conf"
        params = {
            "Threads": 0,
            "Orbits": 1000,
            "Steps": 1000,
            "Width": 100,
            "StepSize": -100,
            "TimeLimit": 100,
            "Criterion": "Hoogerwerf",
            "Limit": round(self.rPhy.value,3),
            #cluster line
            "Star1": f'"/astro/surodeep/traceback/cluster_runaway/{self.name}/{self.name}_fs4giesler.tsv#2"',
            #rest of the stars
            "Star2": f'"/astro/surodeep/traceback/cluster_runaway/{self.name}/{self.name}_fs4giesler.tsv#3..{lastline}"',
            "Assoc": 0,
            "OutFile": f'"/astro/surodeep/traceback/cluster_runaway/{self.name}/runaways/run"'
        }
        update_config_template(config_file_path, output_file_path, **params)
    
        print(f'./traceback ../../cluster_runaway/{self.name}/{self.name}_trace.conf')
        
    def prepare_trace3d(self):
        mainfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway3d/{self.name}"
        runawayfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway3d/{self.name}/runaways"
        trajfolder = f"/home/surodeep/suro_aiu/traceback/cluster_runaway3d/{self.name}/runaway_trajectories"
        os.mkdir(mainfolder) if not os.path.exists(mainfolder) else None
        os.mkdir(runawayfolder) if not os.path.exists(runawayfolder) else None
        os.mkdir(trajfolder) if not os.path.exists(trajfolder) else None

        def update_config_template(config_file_path, output_file_path, **kwargs):
            # Read the template config file
            with open(config_file_path, 'r') as file:
                config_content = file.read()
            # Replace placeholders with actual values
            for key, value in kwargs.items():
                placeholder = f"val_{key}"
                config_content = config_content.replace(placeholder, str(value))
            # Write the updated content to the output file
            with open(output_file_path, 'w') as file:
                file.write(config_content)
            
        fs4giesler = self.fs4giesler3d(outlocation='remote')
        lastline = len(fs4giesler)+1
        # Example usage:
        #starno=5
        #template config
        config_file_path = "/home/surodeep/suro_aiu/traceback/cluster_runaway3d/template_config_trace.conf"
        #output the modified config
        output_file_path = f"/home/surodeep/suro_aiu/traceback/cluster_runaway3d/{self.name}/{self.name}_trace.conf"
        params = {
            "Threads": 0,
            "Orbits": 1000,
            "Steps": 1000,
            "Width": 100,
            "StepSize": -100,
            "TimeLimit": 100,
            "Criterion": "Hoogerwerf",
            "Limit": round(self.rPhy.value,3),
            #cluster line
            "Star1": f'"/astro/surodeep/traceback/cluster_runaway3d/{self.name}/{self.name}_fs4giesler.tsv#2"',
            #rest of the stars
            "Star2": f'"/astro/surodeep/traceback/cluster_runaway3d/{self.name}/{self.name}_fs4giesler.tsv#3..{lastline}"',
            "Assoc": 0,
            "OutFile": f'"/astro/surodeep/traceback/cluster_runaway3d/{self.name}/runaways/run"'
        }
        update_config_template(config_file_path, output_file_path, **params)
    
        print(f'./traceback ../../cluster_runaway3d/{self.name}/{self.name}_trace.conf')
    def runaways_all(self):
        #runaways from giesler traceback
        outputs = os.listdir(f"{config['runaways_path']}{self.name}/runaways/")
        print("getting from",config['runaways_path'] )
        linenos = []
        for output in outputs:
            #print(output)
            if 'run' in output:
                linenos.append(int(output.split("+")[1].replace(".out","")))
        linenos.sort()
        print(linenos)
        fs4giesler = Table.read('/home/surodeep/suro_aiu/traceback/cluster_runaway3d/Basel_8/Basel_8_fs4giesler.tsv', format='ascii.tab')
        runaways_all = fs4giesler[np.array(linenos)-2]
        runaways_all['Source'] = runaways_all['Source'].astype(np.int64)
        runaways_all['RV'] = runaways_all['RV'].astype(np.float64)
        runaways_all['e_RV'] = runaways_all['e_RV'].astype(np.float64)
        
        sir = self.stars_in_region()
        
        return runaways_all
    
    def runaways(self,params=None,temp_threshold=10000):
        params = params or {}
        runaways = estimate_temperature(self.runaways_all(), self.theoretical_isochrone(params=params))
        runaways = runaways[
                            "RA_ICRS_1", "DE_ICRS_1", "rgeo", "Teff", "Temp. Est","v_pec","v_pec3d", "HIP", "TYC2", "Source", "Plx", "e_Plx", "pmRA", "pmDE", "e_pmRA", "e_pmDE", "RUWE", 
                            "Gmag", "BP-RP", "BPmag", "RPmag", "b_rgeo", "B_rgeo", "e_Gmag", "e_BPmag", "e_RPmag", "e_BP-RP", "SkyCoord", 
                            "rmRA","e_rmRA", "rmDE", "e_rmDE", "logg", "RV", "e_RV","rRV", "e_rRV", "FG", "e_FG", "FBP", "e_FBP", "FRP", "e_FRP", "RAVE5", "RAVE6"
                        ]
        
        # runaways['v_pec3d'] = np.sqrt(runaways['v_pec']**2+runaways['RV']**2)
        mask_fast2d = runaways['v_pec'] > 17.6*u.km/u.s
        mask_vpec3ddosentexist = runaways['v_pec3d'].mask
        mask_fast3d = runaways['v_pec3d'] > 25*u.km/u.s
        mask_temp = runaways['Temp. Est'] >= temp_threshold*u.K
        
        runaways = runaways[(mask_fast2d & mask_vpec3ddosentexist) | mask_fast3d & mask_temp]
        runaways.sort('Temp. Est', reverse=True)
        return runaways
    
    def theoretical_isochrone(self, params=None, returnparams=False, parsec_version=2):
        # self.members()
        params = params or {}
        Av = float(params.get('Av', None)) if params.get('Av') is not None else round(float(self.Av.value), 1)
        logage = float(params.get('logage', None)) if params.get('logage') is not None else round(float(self.logage), 1)
        FeH = float(params.get('FeH', None)) if params.get('FeH') is not None else round(float(self.FeH), 1)
        
        theo_iso_path = f"./Clusters/{self.name}/{self.name}_compare_data_out_Av{str(Av)}_logage{str(logage)}_FeH{str(FeH)}.isochrone{parsec_version}"
        # print(Av, logage, FeH)
        if os.path.exists(theo_iso_path):
            theo_iso = Table.read(theo_iso_path, format="ascii")
        else:
            theo_iso = get_theoretical_isochrone(Av=Av, logage=logage, FeH=FeH, parsec_version=parsec_version)
            theo_iso['Gmag'] = theo_iso['Gmag'] + 5 * np.log10(self.distance.value) - 5
            theo_iso['G_BP'] = theo_iso["G_BP"]+ 5 * np.log10(self.distance.value) - 5
            theo_iso['G_RP'] = theo_iso["G_RP"]+ 5 * np.log10(self.distance.value) - 5
            if parsec_version==1.2:
                theo_iso['Teff0'] = 10**theo_iso['logTe']
            theo_iso = theo_iso["Mass", "Teff0", "BP-RP", "Gmag", "G_BP", "G_RP", "logg", "logAge", "logL", "logTe", "Mini"]
            theo_iso.write(theo_iso_path, format="ascii", overwrite=True)
        if returnparams:
            return theo_iso, (Av, logage, FeH)
        else:
            return theo_iso
    
    def plot_traceback_clean(self, star_tables=[]):
        warnings.simplefilter('ignore', ErfaWarning)
        if len(star_tables) == 0:
            star_tables.append(self.runaways())
        # Open the FITS file and extract the image and WCS
        cluster_10pc_fits_path = f'./Clusters/{self.name}/{self.name}_extra10pc.fits'
        if not os.path.exists(cluster_10pc_fits_path):
            get_search_region(self)
        with fits.open(cluster_10pc_fits_path) as fits_file:
            image = fits_file[0]
            wcs = WCS(image.header)
            fig, ax = plt.subplots(subplot_kw={'projection': wcs}, figsize=(10, 10))

            ax.imshow(image.data, cmap='gray', alpha=0.5)
            ax.set_xlabel('Right Ascension (hms)', color="white")
            ax.set_ylabel('Declination (degrees)', color="white")

            # Set the background color to black
            fig.patch.set_facecolor('black')
            ax.set_facecolor('black')
            # Set text colors to white
            ax.title.set_color('white')
            ax.tick_params(axis='both', colors='white')
            # Plot the cluster region
            c = self.skycoord
            radius = self.r50
            region = CircleSkyRegion(c, radius)
            region_pix = region.to_pixel(wcs)
            region_pix.plot(ax=ax, color='red', lw=2)
            def plot_traces(ax, allrun, alpha=0.5):
                allrun_coord_now = SkyCoord(ra=allrun['RA_ICRS_1'], 
                        dec=allrun['DE_ICRS_1'],
                        distance=allrun['rgeo'], 
                        pm_ra_cosdec=allrun['rmRA'],
                        pm_dec=allrun['rmDE'],
                        obstime=Time('J2000')+500*u.kyr)

                allrun_coord_earlier = allrun_coord_now.apply_space_motion(dt=-100*u.kyr)
                # Calculate the pixel coordinates of the runaway stars
                allrun_pixels_now = wcs.world_to_pixel(allrun_coord_now)
                allrun_pixels_earlier = wcs.world_to_pixel(allrun_coord_earlier)

                # Plot the current positions as scatter points
                scatter_main = ax.scatter(allrun_pixels_now[0], allrun_pixels_now[1], 
                                        c=allrun['Temp. Est'], cmap='rainbow_r', edgecolor='white',linewidth=0.5,
                                        zorder=5,alpha=alpha,
                                        s=30,norm=plt.Normalize(3000, 18000))
                # Plot the lines showing motion errors
                if len(allrun)>0:
                    runaway_00, runaway_apdp,runaway_apdm,runaway_amdp,runaway_amdm = [SkyCoord(ra=allrun['RA_ICRS_1'],dec=allrun['DE_ICRS_1'], pm_ra_cosdec=allrun['rmRA'],pm_dec=allrun['rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                                        SkyCoord(ra=allrun['RA_ICRS_1'],dec=allrun['DE_ICRS_1'], pm_ra_cosdec=(allrun['rmRA']+allrun['e_rmRA']),pm_dec=allrun['rmDE']+allrun['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                                        SkyCoord(ra=allrun['RA_ICRS_1'],dec=allrun['DE_ICRS_1'], pm_ra_cosdec=(allrun['rmRA']+allrun['e_rmRA']),pm_dec=allrun['rmDE']-allrun['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                                        SkyCoord(ra=allrun['RA_ICRS_1'],dec=allrun['DE_ICRS_1'], pm_ra_cosdec=(allrun['rmRA']-allrun['e_rmRA']),pm_dec=allrun['rmDE']+allrun['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr)),
                                                                        SkyCoord(ra=allrun['RA_ICRS_1'],dec=allrun['DE_ICRS_1'], pm_ra_cosdec=(allrun['rmRA']-allrun['e_rmRA']),pm_dec=allrun['rmDE']-allrun['e_rmDE'], frame='icrs',obstime=(Time('J2000')+1*u.Myr))]

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
                        ax.add_patch(ellipse)

                    for c in [earlier_runaway_apdp_px,earlier_runaway_apdm_px,earlier_runaway_amdp_px,earlier_runaway_amdm_px]:
                        delta_x = c[0] - runaway_00_px[0]
                        delta_y = c[1] - runaway_00_px[1]
                        # Draw the vectors
                        ax.quiver(runaway_00_px[0], runaway_00_px[1], delta_x, delta_y, angles='xy', scale_units='xy', scale=1, color='limegreen', width=0.001)
                ############################
            
            
            # for start, end in zip(np.transpose(allrun_pixels_now), np.transpose(allrun_pixels_earlier)):
            #     ax.plot([start[0], end[0]], [start[1], end[1]], color='blue')
                return
        
        for star_table in star_tables:
            star_table['rmRA'] = star_table['pmRA']-self.pm_ra_cosdec
            star_table['rmDE'] = star_table['pmDE']-self.pm_dec
            star_table['e_rmRA'] = star_table['e_pmRA']+self.e_pm_ra_cosdec
            star_table['e_rmDE'] = star_table['e_pmDE']+self.e_pm_dec
            star_table['rRV'] = star_table['RV']-self.RV
            star_table['e_rRV'] = star_table['e_RV']+self.e_RV
            star_table['Temp. Est'] = 0
            plot_traces(ax, star_table)

        if len(self.runaways())>0:
            plot_traces(ax, self.runaways(),alpha=1)

        
        plt.tight_layout()
        plt.show()

class Isochrone:
    def __init__(self, cluster, Av=None, logage=None, FeH=None):
        self.cluster = Cluster(cluster.name)
        self.clusterdias = ClusterDias(cluster.name)
        self.Av = round(Av, 1) if Av is not None else cluster.Av.value
        self.logage = round(logage, 1) if logage is not None else cluster.logage
        self.FeH = round(FeH, 1) if FeH is not None else cluster.FeH
        self.theoretical_isochrone, self.params = self.cluster.theoretical_isochrone(
            {'Av': self.Av, 'logage': self.logage, 'FeH': self.FeH}, returnparams=True
        )

    def plot(self, ax):
        if self.Av == self.cluster.Av.value and self.logage == self.cluster.logage and self.FeH == self.cluster.FeH:
            label = f'Av={self.Av:.1f}, logage={self.logage:.1f}, FeH={self.FeH:.1f} (Teff)'
        
        elif self.Av == self.clusterdias.Av.value and self.logage == self.clusterdias.logage and self.FeH == self.clusterdias.FeH:
            label = f'Av={self.Av:.1f}, logage={self.logage:.1f}, FeH={self.FeH:.1f} (Dias)'

        else:
            label = f'Av={self.Av:.1f}, logage={self.logage:.1f}, FeH={self.FeH:.1f}'

        ax.plot(
            self.theoretical_isochrone['BP-RP'],
            self.theoretical_isochrone['Gmag'],
            label=label
        )

def plot_cmd(cluster, isochrones=[], **kwargs):
    """
    ### Example usage
    >>> cl = Cluster("Berkeley_97")
    >>> cld = ClusterDias("Berkeley_97")
    >>> isochrone1 = Isochrone(cld)
    >>> isochrone2 = Isochrone(cl, Av=3, logage=7)
    >>> plot_cmd(cl, isochrones=[isochrone1,isochrone2])
    """
    plt.clf()
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)")
    ax.set_ylabel(r"$G$ (mag)")
    ax.set_title(f"CMD for {cluster.name}")

    #main isochrone for temp
    isochrones.reverse()
    isochrones.append(Isochrone(cluster))
    isochrones.reverse()
        
    for isochrone in isochrones:
        isochrone.plot(ax)

    # Plot cluster members
    mymembers = cluster.mymembers
    ax.errorbar(
        mymembers['BP-RP'], mymembers['Gmag'],
        color='black', zorder=2, fmt='o',
        xerr=mymembers['e_BP-RP'] + 0.02, yerr=mymembers['e_Gmag'],
        label=rf'{len(mymembers)} cluster members'
    )

    # Plot stars in region
    cluster = Cluster(cluster.name)
    stars_in_region = cluster.stars_in_region()
    ax.scatter(
        stars_in_region['BP-RP'], stars_in_region['Gmag'],
        s=2, color='grey', zorder=1, label=f"{len(stars_in_region)} stars in the region"
    )
    
    #scatter find_cluster
    cir = find_cluster(cluster.stars_in_region())
    ax.scatter(
        cir['BP-RP'], cir['Gmag'],
        s=10, color='blue', zorder=3, label=f"{len(cir)} kinematic members found"
    )

    # Plot runaways
    theoretical_isochrone_temp = cluster.theoretical_isochrone()
    runaways = cluster.runaways()
    runaways = estimate_temperature(runaways, theoretical_isochrone_temp)
    scatter_runaways = ax.scatter(
        runaways['BP-RP'], runaways['Gmag'],
        s=50, zorder=4,
        c=runaways['Temp. Est'],
        cmap='spring_r', norm=plt.Normalize(4000, 23000),
        label=f'{len(runaways)} runaway(s)'
    )
    


    # Add colorbar
    colorbar = fig.colorbar(scatter_runaways, ax=ax)
    colorbar.set_label('Temperature (K)')

    # Add cluster parameters table
    cluster_table = [
        ['N', len(cluster.mymembers)],
        [r'$[Fe/H]$', cluster.FeH],
        ['log(Age)', cluster.logage],
        ['Av (mag)', round(cluster.Av.value, 2)],
        ['Dist. (pc)', str(round(cluster.distance.value)) + "$\pm$" + f'{cluster.all["e_Dist"]}']
    ]

    if 'FeH' in kwargs and kwargs['FeH'] != cluster.FeH:
        cluster_table[1][1] = f'{cluster.FeH:.2f} --> {kwargs["FeH"]}'
    if 'logage' in kwargs and kwargs['logage'] != cluster.logage:
        cluster_table[2][1] = f'{cluster.logage:.2f} --> {kwargs["logage"]}'
    if 'Av' in kwargs and kwargs['Av'] != cluster.Av:
        cluster_table[3][1] = f'{cluster.Av.value:.2f} --> {kwargs["Av"]}'

    table_bbox = [0.0, 0.84, 0.44, 0.16]  # [left, bottom, width, height]
    table = ax.table(cellText=cluster_table, cellLoc='right', loc='upper left', bbox=table_bbox)

    for key, cell in table._cells.items():
        cell.set_linewidth(0.5)
        cell.set_edgecolor('lightgray')

    # Set plot limits and invert y-axis
    ax.set_ylim(bottom=min(cluster.theoretical_isochrone()['Gmag']) - 2.5, top=18)
    ax.set_xlim(left=min(cluster.theoretical_isochrone()['BP-RP']) - 0.5, right=4)
    ax.invert_yaxis()
    ax.legend()
    return 

def get_main_name(source):
    warnings.filterwarnings("ignore", category=UserWarning)
    t = Simbad.query_object(f"Gaia DR3 {source}")
    if t is None:
        return f"Gaia DR3 {source}"
    else:
        return (t['MAIN_ID'].value[0])
    
def estimate_temperature(stars, theoretical_isochrone):
    stars['Temp. Est'] = stars['Teff']
    
    Ttheo = theoretical_isochrone['Teff0']
    bprptheo = theoretical_isochrone['BP-RP']
    gmagtheo = theoretical_isochrone['Gmag']
    for star in stars:
        differences_bprp = abs(bprptheo - star['BP-RP'])
        differences_gmag = abs(gmagtheo - star['Gmag'])
        # differences = differences_bprp**2+differences_gmag**2 #method 1
        differences = differences_bprp #method 2
        closest_star_index = np.argmin(differences)
        new_closest_star_temperature = Ttheo[closest_star_index]
        star['Temp. Est']=new_closest_star_temperature
    return stars
def find_cluster(stars_in_region: astropy.table.Table, sigma: float = config['find_cluster']['sigma_clip']) -> astropy.table.Table: 
    """
    Takes input as an astropy table with one column called SkyCoord conatining the ra,dec,distance,pm_ra_cosdec,pm_dec.
    Sigma clips the outliers unit the desired sigma value, saving only the cluster members.
    make sure to put only the ones with >10 plx_quality to ensure effictiveness
    ## input:
    astropy table with one column called SkyCoord conatining the ra,dec,distance,pm_ra_cosdec,pm_dec
    ## returns:
    astropy table with cluster members, sorted according to decreadsing brightness
    ## Example Usage:
    >>> cluster = Cluster('Berkeley_97') 
    >>> sir = cluster.stars_in_region() #get a stars in region table
    >>> my_cluster = find_cluster(sir)
    >>> display(my_cluster)
    """
    mask_clip_RA = ~sigma_clip(stars_in_region['SkyCoord'].ra.value,sigma=1.1*sigma).mask
    mask_clip_DE = ~sigma_clip(stars_in_region['SkyCoord'].dec.value,sigma=1.1*sigma).mask
    mask_clip_pmRA = ~sigma_clip(stars_in_region['SkyCoord'].pm_ra_cosdec.value,sigma=sigma).mask
    mask_clip_pmDE = ~sigma_clip(stars_in_region['SkyCoord'].pm_dec.value,sigma=sigma).mask
    mask_clip_rgeo = ~sigma_clip(stars_in_region['SkyCoord'].distance.value,sigma=sigma).mask
    # mask_plx_quality = ~sigma_clip(stars_in_region['Plx']/stars_in_region['e_Plx'].value,sigma=sigma).mask
    my_stars = stars_in_region[
                            #    mask_plx_quality&
                               mask_clip_RA&
                               mask_clip_DE&
                               mask_clip_pmRA&
                               mask_clip_pmDE&
                               mask_clip_rgeo
                               ]
    # print(f"{len(my_stars)} kinematic members")
    my_stars.sort("Gmag")
    return my_stars
def search_arcmin(distance, radius:Angle, extra=config['Cluster']['search_extent']):
    theta = radius
    D = distance
    r = np.tan(theta) * D #physical radius
    search_arcmin = np.arctan((r + extra * u.pc) / D)
    search_arcmin = search_arcmin.to(u.arcminute)
    return search_arcmin.round(3)
def searchDR3(skycoordinate: SkyCoord, radius: Angle) -> Table:
    start_time = time.time()
    print(f"Searching in Gaia DR3 I/355/gaiadr3 {radius} around {skycoordinate.ra, skycoordinate.dec, skycoordinate.distance}")

    filters = {'Gmag': '<17', 'Plx': '>0.3'}
    stars_fromDR3 = Vizier(columns=["*", "+_r"], row_limit=-1).query_region(
        skycoordinate, 
        radius=radius, 
        catalog="I/355/gaiadr3",
        column_filters=filters
    )[0]
    
    stars_fromDR3['SkyCoord1'] = SkyCoord(
        ra=stars_fromDR3['RA_ICRS'],
        dec=stars_fromDR3['DE_ICRS'],
        distance=(stars_fromDR3['Plx']).to(u.pc, equivalencies=u.parallax()),
        pm_ra_cosdec=stars_fromDR3['pmRA'],
        pm_dec=stars_fromDR3['pmDE'],
        obstime=(Time('J2000') + 1 * u.Myr)
    )

    end_time = time.time()
    print(f"found {len(stars_fromDR3):,} sources in {end_time - start_time:.2f} seconds")
    return stars_fromDR3
def searchDR3_dist(skycoordinate: SkyCoord, radius: Angle) -> Table:
    start_time = time.time()
    print(f"Searching in Gaia DR3 distances I/352/gedr3dis {radius} around {skycoordinate.ra, skycoordinate.dec, skycoordinate.distance}")


    stars_fromDR3_dist = Vizier(columns=["*", "+_r"], row_limit=-1).query_region(
        skycoordinate, 
        radius=radius, 
        catalog="I/352/gedr3dis"
    )[0]
    
    stars_fromDR3_dist['SkyCoord2'] = SkyCoord(
        ra=stars_fromDR3_dist['RA_ICRS'],
        dec=stars_fromDR3_dist['DE_ICRS'],
        distance=(stars_fromDR3_dist['rgeo']),
        obstime=(Time('J2000') + 1 * u.Myr)
    )

    end_time = time.time()
    print(f"found {len(stars_fromDR3_dist):,} sources in {end_time - start_time:.2f} seconds")
    filtered_catalog = stars_fromDR3_dist[skycoordinate.separation_3d(stars_fromDR3_dist['SkyCoord2']) < skycoordinate.distance*config['Cluster']['distance_tolerance']]
    return filtered_catalog
def merge_gaia_tables(stars_fromDR3: Table, stars_fromDR3_dist: Table) -> Table:
    start_time = time.time()
    print("Starting merge of DR3 and distance catalog data")

    warnings.simplefilter('ignore', MergeConflictWarning)
    merged = join(stars_fromDR3, stars_fromDR3_dist, keys='Source', join_type='inner')

    # Order the columns
    merged = merged['RA_ICRS_1', 'DE_ICRS_1', 'e_RA_ICRS', 'e_DE_ICRS', '_r_1',
                                       'HIP', 'TYC2', 'Source', 'rgeo', 'Plx', 'e_Plx',
                                       'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE',
                                       'RUWE', 'Teff', 'logg', 'Gmag', 'BP-RP', 'BPmag', 'RPmag', 'RV', 'e_RV',
                                       'b_rgeo', 'B_rgeo', 'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP', 'RAVE5', 'RAVE6']

    # Adding row for uncertainties in Gmag and BPmag and RPmag
    # values for Gaia G, G_BP, G_RP zero point uncertainties
    #refer https://astronomy.stackexchange.com/questions/38371/how-can-i-calculate-the-uncertainties-in-magnitude-like-the-cds-does
    #refer https://www.cosmos.esa.int/web/gaia/edr3-passbands
    sigmaG_0 = 0.0027553202
    sigmaGBP_0 = 0.0027901700
    sigmaGRP_0 = 0.0037793818

    merged['e_Gmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FG'] / merged['FG'])**2 + sigmaG_0**2)
    merged['e_BPmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FBP'] / merged['FBP'])**2 + sigmaGBP_0**2)
    merged['e_RPmag'] = np.sqrt((-2.5 / np.log(10) * merged['e_FRP'] / merged['FRP'])**2 + sigmaGRP_0**2)
    merged['e_BP-RP'] = merged['e_BPmag'] + merged['e_RPmag']

    merged['SkyCoord'] = SkyCoord(
        ra=merged['RA_ICRS_1'],
        dec=merged['DE_ICRS_1'],
        distance=(merged['rgeo']),
        pm_ra_cosdec=merged['pmRA'],
        pm_dec=merged['pmDE'],
        obstime=(Time('J2000') + 1 * u.Myr)
    )

    end_time = time.time()
    print(f"{len(merged):,} sources found by merging in {end_time - start_time:.2f} seconds")
    return merged
def get_theoretical_isochrone(Av=None,logage=None,FeH=None,parsec_version=2):
    print(f"getting isochrone from cmd3.7 (PARSEC{parsec_version}) with:\nAv:{Av:.2f}\nlogage:{logage:.2f}\nmetallicity:{FeH:.2f}")
    start_time = time.time()
    s = Service()
    # options = webdriver.ChromeOptions()
    browser = webdriver.Chrome(service=s)#, options=options)
    browser.get('http://stev.oapd.inaf.it/cgi-bin/cmd')

    #Evolutionary Tracks #from config
    if parsec_version==2:
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[3]/td[1]/input[1]").click() #PARSEC version 2.0
    elif parsec_version==1.2:
        browser.find_element(By.XPATH,"/html/body/form/div/fieldset[1]/table/tbody/tr[4]/td/input").click() #PARSEC version 1.2S
    
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
    initialMetallicity_field.send_keys(str(FeH))
    finalMetallicity_field = browser.find_element(By.XPATH,"/html/body/form/div/fieldset[7]/table/tbody/tr[7]/td[3]/input")
    finalMetallicity_field.clear()
    finalMetallicity_field.send_keys(str(FeH))
    #Submit #dosen't change
    browser.find_element(By.XPATH,"/html/body/form/div/input[4]").click()
    browser.find_element(By.XPATH,"/html/body/form/fieldset[1]/p[1]/a").click()
    data_out = browser.find_element(By.XPATH,'/html/body/pre').text
    with open('temp_isochrone', 'w') as f:
        f.write(data_out)
    browser.close()
    theoretical_data = Table.read("temp_isochrone", format='ascii')
    # os.remove("temp_isochrone")
    # Define the new column names
    if parsec_version==2:
        new_column_names = ['Zini', 'MH', 'logAge', 'Mini', 'int_IMF', 'Mass', 'logL', 'logTe', 'logg', 'label', 'McoreTP', 'C_O', 'period0', 'period1', 'period2', 'period3', 'period4', 'pmode', 'Mloss', 'tau1m', 'X', 'Y', 'Xc', 'Xn', 'Xo', 'Cexcess', 'Z', 'Teff0', 'omega', 'angvel', 'vtaneq', 'angmom', 'Rpol', 'Req', 'mbolmag', 'G_fSBmag', 'G_BP_fSBmag', 'G_RP_fSBmag', 'G_fSB', 'G_f0', 'G_fk', 'G_i00', 'G_i05', 'G_i10', 'G_i15', 'G_i20', 'G_i25', 'G_i30', 'G_i35', 'G_i40', 'G_i45', 'G_i50', 'G_i55', 'G_i60', 'G_i65', 'G_i70', 'G_i75', 'G_i80', 'G_i85', 'G_i90', 'G_BP_fSB', 'G_BP_f0', 'G_BP_fk', 'G_BP_i00', 'G_BP_i05', 'G_BP_i10', 'G_BP_i15', 'G_BP_i20', 'G_BP_i25', 'G_BP_i30', 'G_BP_i35', 'G_BP_i40', 'G_BP_i45', 'G_BP_i50', 'G_BP_i55', 'G_BP_i60', 'G_BP_i65', 'G_BP_i70', 'G_BP_i75', 'G_BP_i80', 'G_BP_i85', 'G_BP_i90', 'G_RP_fSB', 'G_RP_f0', 'G_RP_fk', 'G_RP_i00', 'G_RP_i05', 'G_RP_i10', 'G_RP_i15', 'G_RP_i20', 'G_RP_i25', 'G_RP_i30', 'G_RP_i35', 'G_RP_i40', 'G_RP_i45', 'G_RP_i50', 'G_RP_i55', 'G_RP_i60', 'G_RP_i65', 'G_RP_i70', 'G_RP_i75', 'G_RP_i80', 'G_RP_i85', 'G_RP_i90']
    elif parsec_version==1.2:
        new_column_names = ["Zini", "MH", "logAge", "Mini", "int_IMF", "Mass", "logL", "logTe", "logg", "label", "McoreTP", "C_O", "period0", "period1", "period2", "period3", "period4", "pmode", "Mloss", "tau1m", "X", "Y", "Xc", "Xn", "Xo", "Cexcess", "Z", "mbolmag", "Gmag", "G_BPmag", "G_RPmag"]

    # Rename the columns
    for old_name, new_name in zip(theoretical_data.colnames, new_column_names):
        theoretical_data.rename_column(old_name, new_name)
    # Calculate BP-RP and add column at the end
    if parsec_version==2:
        theoretical_data['Gmag'] = theoretical_data['G_fSBmag']
        theoretical_data['G_BP'] = theoretical_data['G_BP_fSBmag']
        theoretical_data['G_RP'] = theoretical_data['G_RP_fSBmag']
        theoretical_data['BP-RP'] = theoretical_data['G_BP_fSBmag']-theoretical_data['G_RP_fSBmag']
    elif parsec_version==1.2:
        # theoretical_data['Gmag'] = theoretical_data['Gmag']
        theoretical_data['G_BP'] = theoretical_data['G_BPmag']
        theoretical_data['G_RP'] = theoretical_data['G_RPmag']
        theoretical_data['BP-RP'] = theoretical_data['G_BPmag']-theoretical_data['G_RPmag']
    end_time = time.time()
    print(f"isochrone downloaded in {end_time-start_time:.1f}s")
    
    return theoretical_data
def get_search_region(cluster, extra=10,display=True,**kwargs):
    """
    Plots and saves the fits file for the region around the given cluster.

    Parameters:
    - extra (float): Additional extent for the search region (default is from config['Cluster']['search_extent']).

    Returns:
    None
    """
    search_arcmin = cluster.search_arcmin
    
    # Define the file path
    fits_file_path = f'./Clusters/{cluster.name}/{cluster.name}_extra{extra}pc.fits'
    
    # Check if the file already exists
    if os.path.exists(fits_file_path):
        print(f'fits image exists in {cluster.name} folder')
        # File exists, no need to download, use the existing file
        images = [fits.open(fits_file_path)]
        # Extract the WCS information
        wcs = WCS(images[0][0].header)
    else:
        # print("downloading fits image 10pc around cluster center...")
        start_time = time.time()
        
        # File doesn't exist, get the image data from SkyView
        images = SkyView.get_images(position=cluster.skycoord,
                                    survey='DSS',
                                    radius=2*search_arcmin,
                                    **kwargs)
        # Extract the WCS information
        wcs = WCS(images[0][0].header)
        hdu = fits.PrimaryHDU(data=images[0][0].data, header=images[0][0].header)
        hdulist = fits.HDUList([hdu])
        end_time = time.time()
        
        # Save the fits file
        hdulist.writeto(fits_file_path, overwrite=True)
        print(f"image downloaded in {(end_time-start_time):1f}s")
        