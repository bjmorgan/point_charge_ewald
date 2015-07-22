class Atom():

    num = 0
  
    def __init__( self, label, q, fixed, allowed_sites ):
        self.label = label
        self.q = q
        self.fixed = fixed
        self.allowed_sites = allowed_sites