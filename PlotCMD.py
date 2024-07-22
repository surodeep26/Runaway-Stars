# Import necessary libraries
# %run runaway_functionsv3

import matplotlib.pyplot as plt
from runaway_functionsv3 import Cluster, ClusterDias, estimate_temperature

# Function to plot isochrones
def plot_isochrones(ax, cluster, **kwargs):
    cluster_dias = ClusterDias(cluster.name)
    theoretical_isochrone_dias, (Av0, logage0, FeH0) = cluster_dias.theoretical_isochrone(returnparams=True)
    isochrone_theo_dias = ax.plot(
        theoretical_isochrone_dias['BP-RP'], 
        theoretical_isochrone_dias['Gmag'], 
        label=f'Av={Av0}, logage={logage0}, FeH={FeH0} (Dias, for Teff)'
    )[0]
    
    theoretical_isochrone_temp, (Av, logage, FeH) = cluster.theoretical_isochrone(kwargs, returnparams=True)
    if not (Av0 == Av and logage0 == logage and FeH0 == FeH):
        isochrone_theo_dias.set_label(f'Av={Av0}, logage={logage0}, FeH={FeH0} (Dias)')
        isochrone_theo = ax.plot(
            theoretical_isochrone_temp['BP-RP'], 
            theoretical_isochrone_temp['Gmag'], 
            label=f'Av={Av}, logage={logage}, FeH={FeH} (for Teff)'
        )[0]
    return isochrone_theo_dias

# Function to plot cluster members
def plot_cluster_members(ax, cluster):
    mymembers = cluster.mymembers
    scatter_members = ax.errorbar(
        mymembers['BP-RP'], mymembers['Gmag'], 
        color='black', zorder=2, fmt='o',
        xerr=mymembers['e_BP-RP']+0.02, yerr=mymembers['e_Gmag'],
        label=rf'{len(mymembers)} cluster members'
    )
    return scatter_members

# Function to plot stars in the region
def plot_stars_in_region(ax, cluster):
    stars_in_region = cluster.stars_in_region()
    scatter_sir = ax.scatter(
        stars_in_region['BP-RP'], stars_in_region['Gmag'],
        s=2, color='grey', zorder=1, label=f"{len(stars_in_region)} stars in the region"
    )
    return scatter_sir

# Function to plot runaways
def plot_runaways(ax, cluster, theoretical_isochrone_temp):
    runaways = cluster.runaways()
    runaways = estimate_temperature(runaways, theoretical_isochrone_temp)
    scatter_runaways = ax.scatter(
        runaways['BP-RP'], runaways['Gmag'],
        s=30, zorder=4,
        c=runaways['Temp. Est'],
        cmap='spring_r', norm=plt.Normalize(4000, 23000),
        label=f'{len(runaways)} runaway(s)'
    )
    return scatter_runaways

# Function to add colorbar
def add_colorbar(fig, scatter_runaways, ax):
    colorbar = fig.colorbar(scatter_runaways, ax=ax)
    colorbar.set_label('Temperature (K)')
    return colorbar

# Function to add cluster parameters table
def add_cluster_parameters_table(ax, cluster):
    cluster_dias = ClusterDias(cluster.name)
    print()
    cluster_table = [
        ['N', len(cluster.mymembers)],
        [r'$[Fe/H]$', cluster.FeH],
        ['log(Age)', cluster.logage],
        ['Av (mag)', round(cluster.Av.value, 2)],
        ['Dist. (pc)', str(round(cluster.distance.value))+"$\pm$"+f'{cluster.all["e_Dist"]}']
    ]

    if cluster.FeH != cluster_dias.FeH:
        cluster_table[1][1] = f'{cluster_dias.FeH:.1f} --> {round(float(cluster.FeH),1)}'
    if cluster.logage != cluster_dias.logage:
        cluster_table[2][1] = f'{cluster_dias.logage:.1f} --> {round(float(cluster.logage),1)}'
    if cluster.Av != cluster_dias.Av:
        cluster_table[3][1] = f'{cluster_dias.Av.value:.1f} --> {round(float(cluster.Av.value),1)}'
    if cluster.distance != cluster_dias.distance:
        cluster_table[4][1] = f'{cluster_dias.distance.value:.1f} --> {round(float(cluster.distance.value),1)}'

    table_bbox = [0.0, 0.84, 0.44, 0.16]  # [left, bottom, width, height]
    table = ax.table(cellText=cluster_table, cellLoc='right', loc='upper left', bbox=table_bbox)

    for key, cell in table._cells.items():
        cell.set_linewidth(0.5)
        cell.set_edgecolor('lightgray')
    return table

# Function to add additional isochrones
def add_isochrone(ax, cluster, Av, logage, FeH, parsec_version=2):
    theoretical_isochrone, _ = cluster.theoretical_isochrone({'Av': Av, 'logage': logage, 'FeH': FeH}, returnparams=True, parsec_version=parsec_version)
    isochrone_plot = ax.plot(
        theoretical_isochrone['BP-RP'], 
        theoretical_isochrone['Gmag'], 
        label=f'Av={Av}, logage={logage}, FeH={FeH}'
    )[0]
    return isochrone_plot



# Main function to plot CMD
def plot_cmd(cluster, multiple=False, **kwargs):
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.set_xlabel(r"$G_{BP}-G_{RP}$ (mag)")
    ax.set_ylabel(r"$G$ (mag)")
    ax.set_title(f"CMD for {cluster.name}")

    isochrone_theo_dias = plot_isochrones(ax, cluster, **kwargs)
    plot_cluster_members(ax, cluster)
    plot_stars_in_region(ax, cluster)

    theoretical_isochrone_temp = cluster.theoretical_isochrone(kwargs)
    scatter_runaways = plot_runaways(ax, cluster, theoretical_isochrone_temp)
    add_colorbar(fig, scatter_runaways, ax)
    add_cluster_parameters_table(ax, cluster)

    # add_isochrone(ax, cluster, Av=3.1, logage=7.1, FeH=-0.1, parsec_version=1.2)
    # add_isochrone(ax, cluster, Av=3.1, logage=7.0, FeH=-0.1, parsec_version=2)
    
    
    ax.set_ylim(bottom=min(theoretical_isochrone_temp['Gmag'])-4, top=18)
    ax.set_xlim(left=min(theoretical_isochrone_temp['BP-RP'])-0.5, right=3)
    ax.invert_yaxis()
    ax.legend()
    
    # Optional: make points clickable (code to be implemented)
    # make_points_clickable(ax, scatter_runaways)

    return fig, ax