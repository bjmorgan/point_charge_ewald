import numpy as np

from point_charge_ewald.atom import Atom

from random import sample

def flatten( l ):
    return [ item for sublist in l for item in sublist ]

class Species():
    
    def __init__( self, label, number, q, fixed, allowed_sites, sites ):
        self.label = label
        self.number = number
        self.q = q
        self.fixed = fixed
        self.allowed_sites = allowed_sites
        self.assign_atoms_to_sites( sites )
 
    def assign_atoms_to_sites( self, sites ):
        sites_to_occupy = [ s for s in sites if not s.is_occupied 
                                             if s.label in self.allowed_sites ] 
        print( "{} available sites for species {}".format( len( sites_to_occupy ), self.label ) )
        assert( len( sites_to_occupy ) >= self.number )
        self.atoms = []
        for site in sample( sites_to_occupy, self.number ):
            this_atom = Atom( self.label, self.q, self.fixed, self.allowed_sites )
            site.occupied_by = this_atom
            self.atoms.append( this_atom )
    