class CrossSectToken(object):
    # CROSSSECTTOKEN Data generated upon drawing a cross-section
    #    Objects of this class contain data that was generated when a 
    #    cross-section is drawn. 
    token = [] # Untyped token created by the drawing tool
    innerCoord # x,y coordinate that is located somewhere inside the cross-section (does not have to be the center        
    
    def __init__(self, innerCoord, token):
        # CROSSSECTTOKEN Constructor of a CrossSectToken with a token
        #    This cross-section token object consists of a coordinate
        #    and a token
        self.innerCoord = innerCoord # TO DO: How to deal with units??
        self.token      = token

