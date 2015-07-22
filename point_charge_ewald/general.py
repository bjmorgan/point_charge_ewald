from point_charge_ewald.ewald import ener

def ener_fortran( cell, ewald_params ):
    charged_sites = [ s for s in cell.sites if s.charge != 0 ]
    q = [ s.charge for s in charged_sites ]
    r = [ s.r for s in charged_sites ]
    return( ener( q = q,
                  r = r,
                  kmax = ewald_params.kmax,
                  eta = ewald_params.eta,
                  rsqmax = ewald_params.rsqmax, 
                  rksqmax = ewald_params.rksqmax,
                  boxlen = cell.boxlen,
                  fullhi = cell.fullhi,
                  fourpicell = cell.fourpicell,
                  hlab2 = cell.hlab2 ) )

