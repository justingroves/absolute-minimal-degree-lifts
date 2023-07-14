load('lift.sage')
load('gen_pts.sage')
load('isom.sage')
load('graymap.sage')
load('witt.sage')
load('gen_matrix.sage')

def hyper_lift(f, n, pols = [], check = False):
    
    #This functions computes the minimal hyperelliptic curve lift of an elliptic curve
    #It takes in the function f(x_0) where y^2 = f(x_0) is a hyperelliptic curve, so
    #f(x_0) must have distinct roots and is assumed to be monic
    
    if not check:
        if f.discriminant() == 0:
            raise ValueError("The curve y^2 = f(x) is not a hyperelliptic curve.")
    
    F = (f.parent().base().fraction_field())
    p = F.characteristic()
    d = f.degree()
    
    Pres = f.parent() #Get x0 for y0^2 = f(x0)
    x0 = Pres.gen()
    
    if len(pols) == 0:
        pols = etapols(p,n)
        
    res_coef = [ [coef] for coef in list(f)[0:len(list(f))-1] ] #Initiate list of Witt vectors with coefficients of f(x0)
    resF = [x0] #Initiate Witt vector x
    resH = [Pres.one()] #Initiate Witt vector y (all entries will be multiplied by y0)
    
    #Determine the Hasse Invariant
    ff = f^((p-1) // 2)
    HI = F(ff.coefficients(sparse = False)[p-1])
    
    if not check:
        if HI == F(0):
            raise ValueError("The coefficient of x^(p-1) in f(x)^((p-1)/2) cannot be 0.")
    
    F_is = [x0]
    
    #The main loop, determine the Witt vectors up to the n-th coordinate
    for i in range(0,n):
        
        #Create the polynomial ring with all unknowns
        M = (d * p^(i+1) - (d - 2))/2 #max degree of F_i
        k = floor(M / p)
        Pi = PolynomialRing(F, ['x_0', 'y_0'] + [ 'a%s_n'%ii for ii in range(d)] + 
                            ['c_%s'%(p*ii) for ii in range(k+1) ] + [ 'Hn' ])

        if i == 0: #If computing the 1st coordinate
            tmppf = ff
        else: #If computing past the 1st coordinate
            tmppf = (tmppf^p)*ff
        
        #Compute dF_n/dx0
        Fi = HI^(-(p^(i+1)-1) / (p-1))*tmppf - sum( [ ((F_is[j])^(p^(i+1-j)-1))*(F_is[j].derivative())
                                                      for j in range(i+1)])
        
        for j in [ ii for ii in range(Fi.degree()+1) if (ii + 1) % p == 0]:
            if Fi.list()[j] != 0:
                raise ValueError("The derivative dF_i/dx0 is not a true derivative.")
        
        Fi = formal_integral(Fi) #integrate dF_n/dx0
        
        # F_is.append(Fi) #Attach F_i for future coordinate lifts
        
        Fi = Fi(Pi.gen(0))
        Fi += sum( Pi.gen((d+2)+j)*(Pi.0)^(p*j) for j in range(k+1)) # if j != i) #REMOVE IF J != I TO CHANGE BACK
        
        #Note: by assumption, f(x) is monic so a_d is 1 in W_l(F_q)
        va = [ [Pi(x) for x in res_coef[j] ] + [Pi.gen(2+j)] for j in range(d)] #Witt vectors a_i with unknown a_i_n
        vF = [ x(Pi.0) for x in resF ] + [Fi] #Witt vector x with unknown x_n
        vG = [ Pi.1*(x(Pi.0)) for x in resH ] + [Pi.1*Pi.gen(d+2+(k+1))] #Witt vector y with unknown H_n
        vone = [Pi.one()] + [Pi.zero() for j in range(i+1)] #Additive identity in W_n(k)
        
        #Use algorithms to compute Witt vector product quickly
        vvars = vF + vG

        GTx = GT( [ [vone, d, 0] ] + [ [ va[d-(j+1)] , d-(j+1), 0] for j in range(d)], 
                 pols=pols, vvars=vvars) #Greenburg Transform of f(x_0)
        GTy = GT( [ [vone, 0, 2] ], pols=pols, vvars=vvars) #Greenburg Transform of y

        #Grab unknown coordinate of Witt vector
        RHS = GTx[i+1]
        LHS = GTy[i+1]
        
        LHS = LHS.coefficient({Pi.gen(d+2+(k+1)):0}) #Remove coefficient of H_n (will divide by them later)
        
        #Converty y_0^2 into f(x_0) in LHS
        deg = LHS.degree(Pi.1)
        tmppf2 = 1
        tmpLHS = LHS.coefficient({Pi.1:0})
        for evn_num in [ j for j in range(2,deg+1) if j % 2 == 0]:
            tmppf2 *= f(Pi.0)
            tmpLHS += LHS.coefficient({Pi.1:evn_num})*tmppf2
        
        RHS -= tmpLHS #Subtract part of LHS now in terms of x_0 and not y_0 to RHS
        
        tmppf2 = (tmppf*f)(Pi.0) #Portion of LHS to divide RHS by
        
        RHS *= Pi(1/2) #Divide RHS by 2 (from coefficient of H_n)
        deg1 = RHS.degree(Pi.0) #Highest degree of RHS in x_0
        deg2 = (d*(p^(i+1)+1) / 2) #Degree of (y_0^2)^((p-1)/2)*(y_0^2), i.e. the divisor
        quo = 0
        
        #Perform long division while degree of the remainder is greater than the degree of the divisor
        while deg1 >= deg2:
            lterm = RHS.coefficient({Pi.0:deg1})*((Pi.0)^Integer(deg1-deg2))
            RHS -= lterm*tmppf2
            quo += lterm
            deg1 = RHS.degree(Pi.0) #Degree of the remainder
            
        #Get coefficients of the remainder (need to force them to be 0)
        vrem = []
        for jj in range(RHS.degree(Pi.0)+1):
            vrem += [ RHS.coefficient({Pi.0:jj})]
            
        neqts = len(vrem)
        
        #Rows: correspond to the coefficients a_0,n, a_1,n; etc.
        #Columns: correspond to what power of x0 the coefficient occurs
        #d+1 : a_0, a_1, etc.; k+1 : c_0, c_p, c_(2*p), etc.
        Mat = Matrix(F, [ [F(vrem[j].coefficient({Pi.gen(kk):1})) for j in range(neqts) ]
                        for kk in [ ii for ii in range(2,(d+1)+(k+1)+1)]])
        
        #Vector of coefficients all evaluated at 0 for all unknowns
        vec = vector(F, [ -vrem[j]([Pi(0) for _ in range(len(Pi.gens()))]) for j in range(neqts)])
        
        vsol = Mat.solve_left(vec) #Solve the corresponding system
        
        #Vector of (now) known values to substitute to get a_i,n, x_n, y_n
        eval_vec = [x0, 0] + [vsol[i] for i in range(len(vsol))] + [0]
        
        for j in range(d):
            res_coef[j].append(eval_vec[j+2]) #Append new a_i,n
        resF.append(Fi(eval_vec))
        resH.append(quo(eval_vec))
        
        F_is.append(Fi(eval_vec))
     
    res_coef.append([F(1)] + [F(0) for j in range(n)])
    
    return res_coef, resF, resH