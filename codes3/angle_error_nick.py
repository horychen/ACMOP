import numpy as np
def angle_error(alpha_star, alpha_actual):
    """
    This function computes the angle error between a desired and actual angle.
    
    Input Variables:
        
        alpha_star = an array like object containing all of the desired vector
                     angles in degrees
        
        alpha_actual = an array like object containing all of the actual vector
                       angles in degrees
                       
    Note: alpha_star and alpha_actual must be the same dimension                  
                       
    Output Variables:
        
        angle_err: the angle error in degrees
        
        
    Example:
        
        alpha_star = [5, 60, 359]
        alpha_actual = [358, 65, 6]
        
        err = angle_error(alpha_star, alpha_actual)
        
        print(err) # --> [-7, 5, 7]
    
    """
    
    #ensure that input variables are numpy arrays
    alpha_star = np.atleast_1d(alpha_star)
    alpha_actual = np.atleast_1d(alpha_actual)
    
    N = alpha_star.shape[0] #determine number of input angles
    
    #first preallocate arrays for unit vectors directed at every angle
    vectors_star = np.zeros((N,2))
    vectors = np.zeros((N,2))
    
    #unit vectors for desired angle
    vectors_star[:,0] = np.cos( np.deg2rad(alpha_star) )
    vectors_star[:,1] = np.sin( np.deg2rad(alpha_star) )
    
    #unit vectors for actual angle
    vectors[:,0] = np.cos( np.deg2rad(alpha_actual) )
    vectors[:,1] = np.sin( np.deg2rad(alpha_actual) )
    
    #determine angle between vectors in degrees (note that this is only the angle magnitude):
    #This is just doing the dot product between all corresponding unit vectors
    error_angle_mag = np.rad2deg( np.arccos( (vectors_star*vectors).sum(axis = 1) ) )
    
    #to determine the error angle direction we can use the cross product between the vectors
    #positive cross product means the desired vector lags the actual vector and the error angle
    #is positive
    sign = np.sign( np.cross(vectors_star, vectors) )

    angle_err = sign*error_angle_mag

    return angle_err

def compute_angle_error(alpha_star, alpha_actual):
    
    #both angle arrays are assumed to be in degrees
    N = alpha_star.shape[0] #determine number of input angles
    
    #first preallocate arrays for unit vectors directed at every angle
    vectors_star = np.zeros((N,2))
    vectors =      np.zeros((N,2))
    
    #unit vectors for desired angle
    vectors_star[:,0] = np.cos( np.deg2rad(alpha_star) )
    vectors_star[:,1] = np.sin( np.deg2rad(alpha_star) )
    #unit vectors for actual angle
    vectors[:,0] = np.cos( np.deg2rad(alpha_actual) )
    vectors[:,1] = np.sin( np.deg2rad(alpha_actual) )
    
    #determine angle between vectors in degrees (note that this is only the angle magnitude):
    #This is just doing the cross product between all corresponding unit vectors
    error_angle_mag = np.rad2deg( np.arccos((vectors_star*vectors).sum(axis = 1)) )
    
    #to determine the error angle direction we can use the cross product between the vectors
    #positive cross product means the desired vector lags the actual vector and the error angle
    #is positive
    sign = np.sign( np.cross(vectors_star, vectors) )
    
    return sign*error_angle_mag

def main():
    
    alpha_star = [5, 60, 359, 150]
    alpha_actual = [358, 65, 6]
    
    err = angle_error(alpha_star, alpha_actual)
    print(err) # -> [-7, 5, 7]

    err = compute_angle_error(np.array(alpha_star), np.array(alpha_actual))
    print(err) # -> [-7, 5, 7]

if __name__ == '__main__':
    
    main()
