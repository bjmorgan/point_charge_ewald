#! /usr/bin/env python3

import numpy as np

from point_charge_ewald.ewald import accept_move
from point_charge_ewald.general import ener_fortran
from point_charge_ewald.cell import Cell
from point_charge_ewald.unit_cell import Unit_Cell
from point_charge_ewald.ewald_params import Ewald_Params
from point_charge_ewald.species import Species

def setup():

    h = np.array( [ [ 1.0, 0.0, 0.0 ],  # unit cell matrix
                    [ 0.0, 1.0, 0.0 ],
                    [ 0.0, 0.0, 1.0 ] ] )

    unitlen = np.array( [ 7.8, 7.8, 7.8 ] ) # unit cell lengths

    nunitcell = np.array( [ 4, 4, 4 ] ) # number of unit cells in the full simulation cell

    sites = { 'F' :  'F_B1.mat',    # labels and input files for sites in the unit cell
              'Li' : 'Li_B1.mat' }

    unit_cell = Unit_Cell( h, unitlen, sites )
    cell = Cell( unit_cell, nunitcell )

    ewald_params = Ewald_Params.default( cell )

    # BEN TODO: This populates the cell sites, so should be a Cell method
    # populate the cell sites with some atoms
    species = [ Species( label='F',  number=256, q=-1,  fixed=False, allowed_sites=( 'F' ),  sites=cell.sites ),
                Species( label='Li', number=256, q=-+1, fixed=False, allowed_sites=( 'Li' ), sites=cell.sites ) ]

    print( cell.charge )
    assert( cell.charge == 0 ) # check cell is charge neutral 

    return cell, ewald_params

if __name__ == "__main__":
    iterations = 100000
    temperature = 0.0
    kT = 3.2e-6 / 300.0 * temperature # atomic units

    cell, ewald_params = setup()
    coulomb_energy = ener_fortran( cell, ewald_params )
    print( -1, coulomb_energy )
    for it in range( iterations ):
        new_cell = cell.move_atom_to_new_site()
        new_coulomb_energy = ener_fortran( new_cell, ewald_params )
        if accept_move( new_coulomb_energy, coulomb_energy, kT ):
            coulomb_energy = new_coulomb_energy
            cell = new_cell
            print( it, coulomb_energy )
        with open( 'lattice_struc.xyz', 'w' ) as f:
            sorted_occupied_sites = sorted( cell.occupied_sites, key=lambda x: x.occupied_by.label )
            for site in sorted_occupied_sites:
                f.write( "{} {} {} {}\n".format( site.occupied_by.label, *( str(f) for f in site.r ) ) )
