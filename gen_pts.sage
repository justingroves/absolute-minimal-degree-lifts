def gen_pts(Ellpt_Crv, num_pts = 0, rnd = False):
    
    """
    Generates a list of affine points from an elliptic or hyperelliptic curve.
    
    Arguments:
        
        Ellpt_Crv - an Elliptic Curve or Hyperellitpic Curve object
            
        num_pts - an Integer, the number of desired points (if 0, return all points)
            
        rnd - a Boolean, if the points are randomly selected or not
    
    Returns:
        
        selected_points - a list of (affine) points on Ellpt_Crv
            
    """
    #Throw an Exception if not enough points on curve
    if num_pts >= len(Ellpt_Crv.points()):
        raise Exception("Can't generate more distinct affine points than are on the curve")    
    
    
    #Remove the point at infinity to only select affine points
    list_of_points = Ellpt_Crv.points().copy()
    for points in list_of_points:
        if points[2] == 0:
            list_of_points.remove(points)
            break
        
    #Initialize return list
    selected_points = []
    
    
    if rnd == True: #Select the affine points randomly
        
        import random #For random.sample()
        
        for pts_on_curve in random.sample(list_of_points, k = num_pts):
            selected_points.append([pts_on_curve[0],pts_on_curve[1]])
    else: #Select the affine points as ordered in SAGE
        if num_pts == 0:
            num_pts = len(list_of_points)
        for i in range(num_pts):
            selected_points.append([list_of_points[i][0],list_of_points[i][1]])
    
    return selected_points
    
def lift_pts(lftd_crv, list_pts):
    
    """
    Takes (affine) points on an elliptic or hyperelliptic curve over a finite field and lifts
    them to (affine) points on the lifted elliptic or hyperelliptic curve over the ring of 
    Witt vectors over the finite field.
    
    Arguments:
        
        lftd_crv - a list specifying the elliptic curve lift; namely the output of the
            lift function (Finotti) or hyper_lift function
            
        list_pts - a list of (affine) points on an elliptic curve (the same curve as
            provided in lftd_crv) in [ (x0, y0), (x1, y1), ... , (xn, yn)] format
            
    Returns:
    
        lifted_points - a list of the lifted (affine) points on the lifted elliptic
            curve over the ring of Witt vectors in [ (f(x0), f(y0)), (f(x1), f(y1)), 
            ... , (f(xn), f(yn))] where f is the elliptic curve lift
    """
    
    #Initialize the lifted coordinates
    x_cor = [None]*len(lftd_crv[1])
    y_cor = [None]*len(lftd_crv[1])
    
    lifted_points = [] #Initialize the list of lifted points
    
    if len(lftd_crv) == 4:
        typ_of_crv = 1 #Ordinary elliptic curve
    else:
        typ_of_crv = 0 #Hyperelliptic curve
    
    for points in list_pts: #Compute the lift of each point provided
        
        for i in range(len(x_cor)): #Compute the coordinates of the Witt vectors
            
            #See notes on lift.sage for where these formula come from
            x_cor[i] = lftd_crv[1+typ_of_crv][i](points[0])
            y_cor[i] = lftd_crv[2+typ_of_crv][i](points[0]) * points[1]
            
            
        lifted_points.append((x_cor.copy(), y_cor.copy()))
        
    return lifted_points