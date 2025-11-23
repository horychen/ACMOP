from pylab import np
class Location2D(object):
    # LOCATION2D Indicates a cross section's location
    def __init__(self, anchor_xy, deg_theta):

        self.anchor_xy = anchor_xy  #Distance from global origin xy coordinate to component's origin xy coordinate
        
        self.theta = deg_theta * np.pi/180  # Angles about global xy axes to 
                                            # rotate component's xy axes in radians

    def transformCoords(self, points, deg_addTheta=0):
        # This function takes in an nx2 array of coordinates of the form
        # [x,y] and returns rotated and translated coordinates. The
        # translation and rotation are described by obj.anchor_xy and
        # obj.theta. The optional "addTheta" argument adds an
        # additional angle of "addTheta" to the obj.theta attribute.

        transCoords = []
        for point in points:

            # 旋转方向和eMach是反一下的，我是Park变换。
            cosT = np.cos( self.theta + (deg_addTheta*np.pi/180) )
            sinT = np.sin( self.theta + (deg_addTheta*np.pi/180) )
            
            transCoords.append( [ points[0]*cosT + points[1]*sinT, 
                                  points[0]*-sinT + points[1]*cosT ] )
        return np.array(transCoords)
