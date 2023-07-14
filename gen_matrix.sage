def gen_matrix(lifted_pts, deg, dvsr_deg):
    
    """
    This function takes a list of lifted points on a lifted hyperelliptic curve along with the
    degree of the base hyperelliptic curve and the degree of the divisor and forms the resulting
    generator matrix. This matrix is the generating matrix for an algebraic code.
    
    Arguments:
        
        lifted_pts - a list of Witt vector points that lie on a lifted elliptic curve 
            with unique reduction modulo the charaterstic of the base field.
            
        deg - an integer, the degree of f(x) where y^2 = f(x) is the base hyperelliptic curve
        
        dvsr_deg - an integer, the degree of the divisor D (namely the point at infinity)
    
    Returns:
    
        gen_mat - the generator matrix for the algebraic code over the ring of Witt vectors
    
    """
    
    F = lifted_pts[0][0][0].parent()
    p = F.characteristic()
    r = F.degree()
    ell = len(lifted_pts[0][0]) #Length of the Witt vectors
    
    if deg == 3: #Elliptic Curve Case
        
        gen_mat = [ [ [F(1)] + [F(0) for _ in range(ell-1)]  for _ in range(len(lifted_pts)) ] ] #Initalize the generator matrix
        
        pwr_x = 1
        for i in range(2,dvsr_deg+1):
            new_row = []
            if is_even(i):
                for pts in lifted_pts:
                    new_row.append(WittPower(pts[0], pwr_x))
                pwr_x += -1
            else:
                for pts in lifted_pts:
                    if i > 1:
                        new_row.append(WittProd( WittPower(pts[0], pwr_x), pts[1]))
                    else:
                        new_row.append(pts[1])
                pwr_x += 2
            
            for j in range(r):
                gen_mat.append( [ [ F.gen()^(j)*coor for coor in output] for output in new_row] )
    
    else: #Hyperelliptic Curve Case
    
        #Form the basis of the vector space L(n_P(P_infinity)) where n_P = dvsr_deg based on
        #the powers of x and y (e.g. (1,0) corresponds to (x^1)*(y^0) in the basis)
        if is_even(deg):
            basis = [(1,0)]
            while 2*basis[-1][0] + deg*basis[-1][1] <= dvsr_deg:
                if basis[-1][0] < (deg-2)/2:
                    basis.append((basis[-1][0]+1, basis[-1][1]))
                else:
                    basis.append((0, basis[-1][1]+1))
        else:
            basis = [(1,0)]
            while 2*basis[-1][0] + deg*basis[-1][1] <= dvsr_deg:
                if 2*basis[-1][0] + deg*basis[-1][1] < deg-1:
                    basis.append( (basis[-1][0]+1, 0))
                else:
                    if basis[-1][1] == 0:
                        if basis[-1][0] == (deg-1)/2:
                            basis.append( (0,1) )
                        else:
                            basis.append( (basis[-2][0]+1,1) )
                    else:
                        basis.append( (basis[-2][0]+1,0) )
        
        del basis[-1] #Delete the last entry that has degree just past dvsr_deg
        
        gen_mat = [ [ [F(1)] + [F(0) for _ in range(ell-1)]  for _ in range(len(lifted_pts)) ] ] #Initalize the generator matrix
        
        for pwrs in basis:
            new_row = []
            if pwrs[0] == 0 or pwrs[1] == 0:
                if pwrs[1] == 0:
                    for pts in lifted_pts:
                        new_row.append(WittPower(pts[0],pwrs[0]))
                else:
                    for pts in lifted_pts:
                        new_row.append(WittPower(pts[1],pwrs[1]))
            else:
                for pts in lifted_pts:
                    new_row.append(WittProd( WittPower(pts[0],pwrs[0]), WittPower(pts[1], pwrs[1]) ) )
            
            #If F is an extension of a finite field with p elements
            for j in range(r):
                gen_mat.append( [ [ F.gen()^(j)*coor for coor in output] for output in new_row] )
                
    return gen_mat