def WittToResidue(vector):
    
    """
    Compute the representation of elements from W_l(F_q) to an element in Z/p^l Z, 
    where q = p^r for some prime p.  If r = 1, this is an isomorphism.  If r > 1, then
    this first computes the trace of each vector component coordinate wise, then 
    applies the isomorphism.
        
        Arguments:
        
            vector - A finite length Witt vector
            
        Returns:
        
            sum - An element of Z/p^l Z that is the image of vector under the isomorphism
    
    """
    
    vector_copy = vector.copy();

    F_q = vector_copy[0].parent(); #Get parent field F_q
    p = F_q.characteristic(); 
    ell = len(vector_copy); #Get length of Witt vector
    
    #If F_q is a finite field extension, compute the trace
    if F_q.degree() > 1:
        vector_copy = WittTrace(vector_copy)
        #for i in range(len(vector_copy)):
            #vector_copy[i] = vector_copy[i].trace()
    
    Z_p = Zp(p, ell, 'fixed-mod', 'terse') #Make p-adic integers for Teichmuller representatives
    teich_rep = Z_p.teichmuller_system() #Get Teichmuller representatives for isomorphism
    teich_rep.insert(0,0) #Add the Teichmuller representative for 0
    
    #Compute the isomorphism
    sum = 0;
    for i in range(ell):
        sum = sum + teich_rep[Integer(vector_copy[i])] * p^i;
    
    sum = Integers(p^ell)(sum); #Coerce to an element of Z/p^l Z
    
    return sum