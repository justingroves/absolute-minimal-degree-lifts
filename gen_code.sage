load('hyper_lift.sage')
load('ReducedForm.sage')

def gen_code(f, n, dvsr_deg, num_pts = 0, BooleanGrayMap = True, rnd = False, check = True, verbose = False, return_num_pts = False):
    
    """
    Generates a generating matrix of a linear code over Z/p^(n+1) Z where p is the
    characteristic of the base field of the polynomial ring over F and (n+1) is the
    length of the Witt vectors from the lift of f. This is done
    by lifting the hyperelliptic curve y^2 = f to a hyperelliptic curve with minimal
    degree over W_(n+1)(F_q), the ring of Witt vectors of length n+1 over the base field
    of the polynomial ring of F (which has size p^r). It then takes a collection of points
    on the base curve and lifts them to the lifted elliptic curve, forms the basis for the
    vector space L(dvsr_deg(P_infinity)) and generates the algebraic code. Then it applies
    the trace map and isomorphism map to make the generating matrix over Z/p^(n+1) Z.
    
    Arguments:
    
        f - a hyperelliptic curve in one variable over a finite field
        
        n - an integer of what length of Witt vector to lift to (will lift to n+1)
        
        dvsr_deg - an integer, the degree of the divisor of the point at infinity of the curve
        
        num_pts - an integer, the number of distinct (affine) points to select from the base curve. If
            set to 0, will select all distinct points
           
        BooleanGrayMap - a Boolean, indicating which Gray map to use
        
        rnd - a Boolean, to indicate to select points on the base curve randomly or not
        
        check - a Boolean, to check if f meets the requirements to be lifted
        
    Returns:
    
        G - a generating matrix for a linear code over Z/p^(n+1) Z of length at most (num_pts)
    
    """
    
    PR = f.parent()
    F = PR.base()
    p = F.characteristic()
    deg = f.degree()
    
    if deg == 3 and F.characteristic() != 3:
        x = PR.gen(0)
        PR.<x,y> = PolynomialRing(F)
        a, b = WeierstrassForm(f - y^2)
        base_curve = EllipticCurve(F,[a,b])
        lifted_curve = lift(a,b,n)
    else:
        x = PR.gen(0)
        PR.<x> = PolynomialRing(F)
        f = PR(f)
        base_curve = HyperellipticCurve(f)
        lifted_curve = hyper_lift(f,n,check = check)
        
    base_pts = gen_pts(base_curve, num_pts = num_pts, rnd = rnd)
    
    if verbose:
        if is_odd(deg):
            alpha = 1
        else:
            alpha = 2

        genus = (deg - alpha) // 2
    
        if (dvsr_deg >= len(base_pts)) or (dvsr_deg <= (2*genus - 2)):
            print("WARNING: generated code may not be a linear code")
        
    lifted_pts = lift_pts(lifted_curve, base_pts)
    
    gen_mat = gen_matrix(lifted_pts, deg, dvsr_deg)
    
    for i in range(len(gen_mat)):
        for j in range(len(gen_mat[i])):
            gen_mat[i][j] = WittToResidue(gen_mat[i][j])
    
    if False:
        if BooleanGrayMap:
            for i in range(len(gen_mat)):
                gen_mat[i] = GrayMapHengYue(gen_mat[i])
        else:
            for i in range(len(gen_mat)):
                gen_mat[i] = GrayMapYildezOzger(gen_mat[i])

        Gray_mat = []
        for i in range(len(gen_mat)):
            new_row = []
            for j in range(len(gen_mat[i])):
                new_row.extend(gen_mat[i][j])
            Gray_mat.append(new_row)

        del gen_mat

        Gray_mat = matrix(Gray_mat)
    
        G = LinearCode(Gray_mat)
    
    row_dim = len(gen_mat)
    col_dim = len(gen_mat[0])
    
    MS = MatrixSpace(Integers(p^(n+1)), row_dim, col_dim)
    gen_mat = MS(gen_mat)
    gen_mat = ReducedForm(gen_mat)
    
    if return_num_pts == False:
        return gen_mat
    else:
        return gen_mat, len(lifted_pts)