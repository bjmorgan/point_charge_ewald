import numpy as np

from point_charge_ewald.ewald import boxreset
from point_charge_ewald.site import Site
from point_charge_ewald.general import ener_fortran

from random import choice
from collections import Counter
import copy

class Cell():

    def __init__( self, unit_cell, nunitcell ):

        self.h = unit_cell.h
        self.nunitcell = nunitcell
        self.unit_cell = unit_cell
        self.unit_len = self.unit_cell.unit_len
        self.boxlen = self.unit_len * self.nunitcell
        self.halfboxrec = 2.0 / self.boxlen
        self.fullh, self.fullhi, self.fourpicell = boxreset( self.h, self.boxlen )
        self.construct_sites( self.unit_cell, self.nunitcell )

    def move_atom_to_new_site( self ):
        new_cell = copy.deepcopy( self )
        site_i = choice( [ s for s in new_cell.sites if s.is_occupied_by_a_mobile_atom ] )
        atom_i = site_i.occupied_by
        possible_sites_to_move_to = [ s for s in new_cell.sites if s.can_be_occupied_by_atom( atom_i )
                                                                    if not s.is_blocked ]
        site_j = choice( possible_sites_to_move_to )
        site_i.occupied_by, site_j.occupied_by = site_j.occupied_by, site_i.occupied_by                                                      
        return new_cell

    def construct_sites( self, unit_cell, nunitcell ):
        self.sites_by_label = {}
        self.sites = []
        for label, coords in unit_cell.site_coordinates.items():
            these_sites = []
            for coord in coords:
                for i in range( self.nunitcell[0] ):
                    for j in range( self.nunitcell[1] ):
                        for k in range( self.nunitcell[2] ):
                            these_sites.append( Site( label, self.unit_len * np.array( [ i, j, k ] ) + coord * self.unit_len ) )
            self.sites_by_label[ label ] = these_sites
            self.sites.extend( these_sites )

    @property
    def hlab2( self ):
        return self.h

    @property
    def occupied_sites( self ):
        return [ s for s in self.sites if s.is_occupied ]

    @property
    def charge( self ):
        return( sum( [ s.charge for s in self.occupied_sites ] ) )

    def dr2( self, r1, r2 ): # apply minimum image convention and convert to lab coordiates.
        dr = r2 - r1

        dr[0] = dr[0] - self.boxlen[0] * int( dr[0] * self.halfboxrec[0] )
        dr[1] = dr[1] - self.boxlen[1] * int( dr[1] * self.halfboxrec[1] )
        dr[2] = dr[2] - self.boxlen[2] * int( dr[2] * self.halfboxrec[2] )

        dr_lab = [ self.h[0,0]*dr[0] + self.h[0,1]*dr[1] + self.h[0,2]*dr[2], 
                   self.h[1,0]*dr[0] + self.h[1,1]*dr[1] + self.h[1,2]*dr[2],
                   self.h[2,0]*dr[0] + self.h[2,1]*dr[1] + self.h[2,2]*dr[2] ]

        return np.dot( dr_lab, dr_lab )

    @property
    def sorted_occupied_sites( self ):
        return sorted( self.occupied_sites, key=lambda x: x.occupied_by.label )

    @property
    def site_labels( self ):
        return [ label for label in self.sites_by_label ]

    @property
    def species_labels( self ):
        return sorted( set( [ site.occupied_by.label for site in self.occupied_sites ] )  )

    @property
    def species_numbers( self ):
        return Counter( [ site.occupied_by.label for site in self.occupied_sites ] )

    def write_xyz( self, filename, title = 'Title' ):
        with open( filename, 'w' ) as f:
            f.write( str( sum( self.species_numbers.values() ) ) + "\n" )
            f.write( "{}\n".format( title ) )
            for site in self.sorted_occupied_sites:
                f.write( "{} {} {} {}\n".format( site.occupied_by.label, *( str(f) for f in site.r ) ) )

    def write_poscar( self, filename, title = 'Title' ):
        with open( filename, 'w' ) as f:
            f.write( "{}\n".format( title ) )
            f.write( "0.52918\n" )
            for row in self.fullh.transpose(): # WARNING ! If this is *not* orthorhombic, check the row / column convention
                f.write( ' '.join( [ str( f ) for f in row ] ) + "\n" )
            f.write( ' '.join( self.species_labels ) + "\n" )
            f.write( ' '.join( [ str(n) for n in self.species_numbers.values() ] ) + "\n" )
            f.write( "Cartesian\n" )
            for site in self.sorted_occupied_sites:
                f.write( ' '.join( [ str(f) for f in site.r ] ) + "\n" )
        
    def coulomb_energy( self, ewald_params ):
        return ener_fortran( self, ewald_params )
 
