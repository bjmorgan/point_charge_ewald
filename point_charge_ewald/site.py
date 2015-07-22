import numpy as np

class Site:

    def __init__( self, label, r ):
        self.label = label
        self.r = r
        self.occupied_by = None

    @property
    def is_occupied( self ):
        if self.occupied_by:
            return True
        else:
            return False

    @property
    def is_occupied_by_a_mobile_atom( self ):
        if self.is_occupied:
            if self.occupied_by.fixed is True:
                return False
            else:
                return True
        else:
            return False

    def can_be_occupied_by_atom( self, atom ):
        return self.label in atom.allowed_sites

    @property
    def is_blocked( self ):
        if self.is_occupied:
            if self.occupied_by.fixed is True:
                return True
        return False

    @property
    def charge( self ):
        if self.is_occupied:
            return self.occupied_by.q
        else:
            return 0.0

