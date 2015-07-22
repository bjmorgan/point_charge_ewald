from point_charge_ewald.ewald import kset

class Ewald_Params():

    def __init__( self, etainpt, rcut, conv1, convfac, cell ):
        self.etainpt = etainpt
        self.rcut = rcut
        self.conv1 = conv1
        self.convfac = convfac
        self.rsqmax = rcut**2.0
        self.eta = etainpt / (2.0 * rcut )
        self.kmax, self.rksqmax = kset( conv1, convfac, self.eta, cell.fullhi ) 

    @classmethod
    def default( cls, cell ):
        return cls( etainpt = 5.6, rcut = 15.0, conv1 = 1.0e-7, convfac = 0.1 , cell = cell )