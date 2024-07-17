from astropy import units as u
from astropy.coordinates import SkyCoord
def calculate_circle(coordinates,extra=10, output=False):
    """

    """
    ra = coordinates.ra
    dec = coordinates.dec
    dist = coordinates.distance

    search_arcmin = np.arctan((extra * u.pc) / dist)
    
    search_arcmin = (search_arcmin.to(u.arcminute)).round(3)

    return search_arcmin