def get_cluster(cluster_name):
    import os
    from astropy.table import Table, Column
    matched_table = Table.read("Clusters_from_dias_and_a99", format='ascii.ecsv')
    mask = (matched_table['Cluster'] == cluster_name)
    cluster = matched_table[mask]
    if os.path.exists(cluster_name):
        print('Folder present for ',cluster_name)
    else:
        os.mkdir(cluster_name)
    return cluster
def refine_my_members(stars_in_region,cluster_name,v_dispersion=None,parallax_quality_threshold=None,diameter=None,parallax_range=None):
    import os
    import pandas
    from astropy.table import Table, Column
    import numpy as np

    #PARALLAX QUALITY
    if parallax_quality_threshold is not None:
        mask_parallaxQuality = [value >= parallax_quality_threshold for value in stars_in_region['Plx']/stars_in_region['e_Plx'] ]
        refined_cluster_members_plxquality = stars_in_region[mask_parallaxQuality]
        print('parallax quality filtered')
    else:
        refined_cluster_members_plxquality = stars_in_region
    #PARALLAX RANGE
    if parallax_range is not None:
        mask_parallaxRange = [(parallax_range[0] < value < parallax_range[1]) for value in refined_cluster_members_plxquality['Plx']]
        refined_cluster_members_plxquality_plxrange = refined_cluster_members_plxquality[mask_parallaxRange]
        print('parallax range filtered')

    else:
        refined_cluster_members_plxquality_plxrange = refined_cluster_members_plxquality
    #DIAMETER
    if diameter is not None:
        mask_clusterDiameter = [d <= diameter for d in refined_cluster_members_plxquality_plxrange['_r']]
        refined_cluster_members = refined_cluster_members_plxquality_plxrange[mask_clusterDiameter]
        print(f'distance from center (angular) filtered within {diameter}')

    else:
        refined_cluster_members = refined_cluster_members_plxquality_plxrange

    refined_cluster_members.to_pandas().to_excel(os.path.join(cluster_name,f'{cluster_name}_my_refined_members.xlsx'), index=False)
    cluster_members = refined_cluster_members
    df = cluster_members.to_pandas()

    df['v_RA'] = (df['pmRA']-df['pmRA'].mean())*4.74*df['rgeo']/1000
    position = df.columns.get_loc('pmRA') 
    df.insert(position, 'v_RA', df.pop('v_RA'))

    df['v_DE'] = (df['pmDE']-df['pmDE'].mean())*4.74*df['rgeo']/1000
    position = df.columns.get_loc('v_RA')
    df.insert(position, 'v_DE', df.pop('v_DE'))

    df['v'] = (df['v_RA']**2+df['v_DE']**2)**0.5
    position = df.columns.get_loc('v_DE') 
    df.insert(position, 'v', df.pop('v'))

    mean_v_RA = df['v_RA'].mean()
    std_v_RA = df['v_RA'].std()
    mean_v_DE = df['v_DE'].mean()
    std_v_DE = df['v_DE'].std()

    mean_v = (mean_v_RA**2+mean_v_DE**2)**0.5

    if v_dispersion is not None:
        mask_pmRA = (df['v_RA'] >= mean_v_RA - v_dispersion) & (df['v_RA'] <= mean_v_RA + v_dispersion)
        mask_pmDE = (df['v_DE'] >= mean_v_DE - v_dispersion) & (df['v_DE'] <= mean_v_DE + v_dispersion)
        mask_v    = (df['v'] >= mean_v - v_dispersion) & (df['v'] <= mean_v + v_dispersion)

        # mask = mask_pmRA & mask_pmDE
        mask = mask_v
        print("Before removing outliers")
        print(df['v_RA'].mean(),df['v_DE'].mean())
        print(df['v_RA'].std(),df['v_DE'].std())
        print('proper motions')
        print(df['pmRA'].mean(),df['pmDE'].mean())
        print(df['pmRA'].std(),df['pmDE'].std())
        filtered_df = df[mask].reset_index(drop=True)
        print('kinematic outliers filtered')

    else:
        filtered_df = df
        print('kinematic outliers not removed')
    print("After removing outliers")
    print(filtered_df['v_RA'].mean(),filtered_df['v_DE'].mean())
    print(filtered_df['v_RA'].std(),filtered_df['v_DE'].std())
    print('proper motions')
    print(filtered_df['pmRA'].mean(),filtered_df['pmDE'].mean())
    print(filtered_df['pmRA'].std(),filtered_df['pmDE'].std())
    refined_cluster_members = Table.from_pandas(filtered_df)
    refined_cluster_members = refined_cluster_members[np.argsort(refined_cluster_members['Gmag'])]
    cluster_members = refined_cluster_members
    return cluster_members
# Visibility check
def observe_from_gsh(clustername,obstime=None):
    import matplotlib.pyplot as plt
    import numpy as np
    import astropy
    from astropy.table import Table
    clusters = Table.read('Clusters_from_dias_and_a99', format='ascii.ecsv')
    from astropy.visualization import astropy_mpl_style, quantity_support
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    from astropy.coordinates import get_sun, get_body
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





#deleteme
def refine_my_members_new(stars_in_region,cluster_name,v_dispersion=None,parallax_quality_threshold=None,diameter=None,parallax_range=None):
    import os
    import pandas
    from astropy.table import Table, Column
    import numpy as np

    #PARALLAX QUALITY
    if parallax_quality_threshold is not None:
        mask_parallaxQuality = [value >= parallax_quality_threshold for value in stars_in_region['Plx']/stars_in_region['e_Plx'] ]
        refined_cluster_members_plxquality = stars_in_region[mask_parallaxQuality]
        print('parallax quality filtered')
    else:
        refined_cluster_members_plxquality = stars_in_region
    #PARALLAX RANGE
    if parallax_range is not None:
        mask_parallaxRange = [(parallax_range[0] < value < parallax_range[1]) for value in refined_cluster_members_plxquality['Plx']]
        refined_cluster_members_plxquality_plxrange = refined_cluster_members_plxquality[mask_parallaxRange]
        print('parallax range filtered')

    else:
        refined_cluster_members_plxquality_plxrange = refined_cluster_members_plxquality
    #DIAMETER
    if diameter is not None:
        mask_clusterDiameter = [d <= diameter for d in refined_cluster_members_plxquality_plxrange['_r']]
        refined_cluster_members = refined_cluster_members_plxquality_plxrange[mask_clusterDiameter]
        print(f'distance from center (angular) filtered within {diameter}')

    else:
        refined_cluster_members = refined_cluster_members_plxquality_plxrange
    # deleteme
    # refined_cluster_members.to_pandas().to_excel(os.path.join(cluster_name,f'{cluster_name}_my_refined_members.xlsx'), index=False)
    cluster_members = refined_cluster_members
    df = cluster_members.to_pandas()

    df['v_RA'] = (df['pmRA']-df['pmRA'].mean())*4.74*df['rgeo']/1000
    position = df.columns.get_loc('pmRA') 
    df.insert(position, 'v_RA', df.pop('v_RA'))

    df['v_DE'] = (df['pmDE']-df['pmDE'].mean())*4.74*df['rgeo']/1000
    position = df.columns.get_loc('v_RA')
    df.insert(position, 'v_DE', df.pop('v_DE'))

    df['v'] = (df['v_RA']**2+df['v_DE']**2)**0.5
    position = df.columns.get_loc('v_DE') 
    df.insert(position, 'v', df.pop('v'))

    mean_v_RA = df['v_RA'].mean()
    std_v_RA = df['v_RA'].std()
    mean_v_DE = df['v_DE'].mean()
    std_v_DE = df['v_DE'].std()

    mean_v = (mean_v_RA**2+mean_v_DE**2)**0.5

    if v_dispersion is not None:
        mask_pmRA = (df['v_RA'] >= mean_v_RA - v_dispersion) & (df['v_RA'] <= mean_v_RA + v_dispersion)
        mask_pmDE = (df['v_DE'] >= mean_v_DE - v_dispersion) & (df['v_DE'] <= mean_v_DE + v_dispersion)
        mask_v    = (df['v'] >= mean_v - v_dispersion) & (df['v'] <= mean_v + v_dispersion)

        # mask = mask_pmRA & mask_pmDE
        mask = mask_v
        print("Before removing outliers")
        print(df['v_RA'].mean(),df['v_DE'].mean())
        print(df['v_RA'].std(),df['v_DE'].std())
        print('proper motions')
        print(df['pmRA'].mean(),df['pmDE'].mean())
        print(df['pmRA'].std(),df['pmDE'].std())
        filtered_df = df[mask].reset_index(drop=True)
        print('kinematic outliers filtered')

    else:
        filtered_df = df
        print('kinematic outliers not removed')
    print("After removing outliers")
    print(filtered_df['v_RA'].mean(),filtered_df['v_DE'].mean())
    print(filtered_df['v_RA'].std(),filtered_df['v_DE'].std())
    print('proper motions')
    print(filtered_df['pmRA'].mean(),filtered_df['pmDE'].mean())
    print(filtered_df['pmRA'].std(),filtered_df['pmDE'].std())
    refined_cluster_members = Table.from_pandas(filtered_df)
    refined_cluster_members = refined_cluster_members[np.argsort(refined_cluster_members['Gmag'])]
    cluster_members = refined_cluster_members
    return cluster_members

def calculate_search_arcmin(cluster, extra=10):
    from IPython.display import Math
    import astropy.units as u
    import numpy as np
    

    theta = cluster['Diameter'][0] * u.arcmin / 2  # radius of the cluster in arcmin
    D = cluster['Dist'][0] * u.pc
    r = np.tan(theta) * D
    search_arcmin = np.arctan((r + extra * u.pc) / D)
    search_arcmin = search_arcmin.to(u.arcminute)
    r, search_arcmin = r.round(3), search_arcmin.round(3)

    display(Math(r'\mathrm{Dist.}' + f'= {D}' + r'\mathrm{\ Radius}' + f'= {r}' +
                 r'\mathrm{\ search \ arcmin}' + f'= {search_arcmin}' +
                 r'\mathrm{\ Cluster \ Radius}' + f'= {theta}'))

    return search_arcmin