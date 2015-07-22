import numpy as np

from point_charge_ewald.site import Site

class Unit_Cell:

    def __init__( self, h, unit_len, sites ):
        self.h = h
        self.unit_len = unit_len
        self.site_coordinates = { label: self.load_site_coordinates( filename ) for label, filename in sites.items() }

    def load_site_coordinates( self, filename ):
        return np.loadtxt( filename )

