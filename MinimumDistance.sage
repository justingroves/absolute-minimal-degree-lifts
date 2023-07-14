def HomogeneousWeight(x, p=None, ell=None, check=True):
    """
        Return the Homogeneous weight of an element from Z/p^l Z for some prime p.
        
        Arguments:
        
            -x: an element of Z/p^l Z
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x is in Z/p^l Z
            
        Returns:
            
            an integer
    """
    
    if x == 0:
        return 0
    
    if check:
        R = x.parent()
        factor_list = list(factor(R.order()))

        if len(factor_list) > 1:
            raise ValueError('Parent ring must be Z/ p^l Z for fixed prime p.')
    
    if (p is None) or (ell is None):
        if not check:
            factor_list = list(factor(R.order()))
        p = factor_list[0][0]
        ell = factor_list[0][1]
    
    if (x % p^(ell - 1)) == 0:
        return p^(ell - 1)
    else:
        return (p-1)*(p^(ell-2))
  
@parallel(8)
def HomogeneousWeightVec(x_vec, p=None, ell=None, check=True):
    """
        Return the Homogeneous weight of a vector over (Z/p^l Z)^n
        
        Arguments:
        
            -x_vec: a vector from (Z/p^l Z)^n
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x_vec is in (Z/p^l Z)^n
            
        Returns:
            
            an integer
    """
    
    return sum([ HomogeneousWeight(x, p=p, ell=ell, check=check) for x in x_vec ])
        
        
def HomogeneousDistanceElm(x,y,p=None, ell=None, check=True):
    """
        Return the Homogeneous distance of two elements in (Z/p^l Z)
        
        Arguments:
        
            -x: an element of (Z/p^l Z)
            
            -y: an element of (Z/p^l Z)
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x,y is in (Z/p^l Z)
            
        Returns:
            
            an integer
    """

    if check:
        R = x.parent()
        S = y.parent()

        if R != S:
            raise ValueError('{} and {} must be from the same ring'.format(x,y))
    
    return HomogeneousWeight(x-y, p=p, ell=ell, check=check)
    
    

def HomogeneousDistanceVec(x_vec, y_vec, p=None, ell=None, check=True):
    """
        Return the Homogeneous distance between two vectors from (Z/p^l Z)^n
        
        Arguments:
        
            -x_vec: a vector from (Z/p^l Z)^n
            
            -y_vec: a vector from (Z/p^l Z)^n
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x_vec,y_vec is in (Z/p^l Z)^n
            
        Returns:
            
            an integer
    """
    
    if len(x_vec) != len(y_vec):
        raise ValueError('Vectors must have the same length')
    
    return sum([ HomogeneousDistanceElm(x, y, p=p, ell=ell,check=check) for x,y in zip(x_vec,y_vec) ])
    

def ExtLeeWeight(x, p=None, ell=None, check=True):
    """
        Return the Extended Lee weight of an element from Z/p^l Z for some prime p.
        
        Arguments:
        
            -x: an element of Z/p^l Z
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x is in Z/p^l Z
            
        Returns:
            
            an integer
    """
    
    if x == 0:
        return 0
    
    if check:
        R = x.parent()
        factor_list = list(factor(R.order()))

        if len(factor_list) > 1:
            raise ValueError('Parent ring must be Z/ p^l Z for fixed prime p.')
    
    if (p is None) or (ell is None):
        if not check:
            factor_list = list(factor(R.order()))
        p = factor_list[0][0]
        ell = factor_list[0][1]
        
    if x <= p^(ell-1):
        return Integer(x)
    elif x <= (p-1)*(p^(ell-1)):
        return p^(ell-1)
    else:
        return p^ell - Integer(x)
        
@parallel(8)
def ExtLeeWeightVec(x_vec, p=None, ell=None, check=True):
    """
        Return the Extended Lee weight of a vector over (Z/p^l Z)^n
        
        Arguments:
        
            -x_vec: a vector from (Z/p^l Z)^n
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x_vec is in (Z/p^l Z)^n
            
        Returns:
            
            an integer
    """
    
    return sum([ ExtLeeWeight(x, p=p, ell=ell, check=check) for x in x_vec ])
        
        
def ExtLeeDistanceElm(x,y, p=None, ell=None, check=True):
    """
        Return the Extended Lee distance of two elements in (Z/p^l Z)
        
        Arguments:
        
            -x: an element of (Z/p^l Z)
            
            -y: an element of (Z/p^l Z)
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x,y is in (Z/p^l Z)
            
        Returns:
            
            an integer
    """
    
    R = x.parent()
    S = y.parent()
    
    if R != S:
        raise ValueError('{} and {} must be from the same ring'.format(x,y))
    
    return ExtLeeWeight(x-y , p=p, ell=ell, check=check)
    

def ExtLeeDistanceVec(x_vec, y_vec, p=None, ell=None, check=True):
    """
        Return the Extended Lee distance between two vectors from (Z/p^l Z)^n
        
        Arguments:
        
            -x_vec: a vector from (Z/p^l Z)^n
            
            -y_vec: a vector from (Z/p^l Z)^n
            
            -p: prime p (from Z/p^l Z)
            
            -ell: integer l (from Z/p^l Z)
            
            -check: Boolean, check if x_vec,y_vec is in (Z/p^l Z)^n
            
        Returns:
            
            an integer
    """
    
    if len(x_vec) != len(y_vec):
        raise ValueError('Vectors must have the same length')
    
    return sum([ ExtLeeDistanceElm(x, y, p=p, ell=ell, check=check) for x,y in zip(x_vec,y_vec) ])
    
@parallel(8)
def MatrixSpan(matrix,vec):
    """
        Return the product of a vector times a matrix on the right hand side. 
        Note: this function is for use where the vector is from a free module over
        a ring so the operation may be parallelized. 
        
        Arguments:
        
            -matrix: a matrix over a ring where the number of rows is equal to the 
                length of vec
                
            -vec: a vector from a free module
            
        Returns:
        
            A vector from the free module over the ring of the matrix with length the number of
                columns of matrix
    """
    return vec*matrix


def MinimumDistance(rdcd_mtrx, metric=0, verbose=False, check=True):
    """
        Given a (reduced) matrix over a free module, compute the minimum distance of
        the code obtained from the matrix using a brute force algorithm. Note 
        that this is the same as the minimum weight of a non-zero codeword.
        
        Arguments:
        
            -rdcd_mtrx : a matrix over a free module over Z/p^l Z in reduced form (see ReducedForm.sage)
            
            -metric : 0 - Homogeneous Distance, 1 - Extended Lee Distance, 2 - Both
            
            -verbose : a Boolean, suppress output
            
            -check : a Boolean, check if entries of matrix are over Z/p^l Z
        
        Returns:
        
            -minimum_distance : an integer
    """    

    if metric < 2:
        if metric == 0:
            metric_str = "Homogeneous"
        else:
            metric_str = "Extended Lee"
            
    R = (rdcd_mtrx.parent()).base_ring()
    factor_list = list(factor(R.order()))

    if len(factor_list) > 1 :
        raise ValueError("Matrix must be over a ring of the form Z/p^l Z for some prime p")

    p = factor_list[0][0]
    ell = factor_list[0][1]

    fm_dim = rdcd_mtrx.dimensions()
    
    #Check if reduced form is the identity matrix
    #if so, return the code that is the whole space
    if fm_dim[0] == fm_dim[1]:
        if rdcd_mtrx == rdcd_mtrx.parent().one():
            if verbose:
                print("--> Reduced matrix is the identity matrix, returning full code...")
                
            if metric == 0:
                return (p^ell)^(fm_dim[0]), (p-1)*(p^(ell-2))
            elif metric == 1:
                return (p^ell)^(fm_dim[0]), 1
            elif metric == 2:
                return (p^ell)^(fm_dim[0]), (p-1)*(p^(ell-2)), 1
            else:
                raise ValueError('metric must be either 0, 1, or 2')
    
    FM_domain = FreeModule(R,fm_dim[0])
    FM_range = FreeModule(R,fm_dim[1])

    codewords = []

    if verbose:
        print("--> Generating codewords from free module...")
        print("---> Generating domain of free module...", end='\r')

    list_of_vec = [ (rdcd_mtrx, vec) for vec in FM_domain]
    
    if verbose:
        print("---> Generating domain of free module...DONE!")
        print("---> Generating list of codewords...", end='\r')

    codewords = MatrixSpan(list_of_vec)
    codewords = [ x[1] for x in codewords ]

    if verbose:
        print("---> Generating list of codewords...DONE!")
        print("---> Isolating unique codewords...", end='\r')

    unique_cdwds = {tuple(word) for word in codewords}
    codewords = [ FM_range(word) for word in unique_cdwds ]
    del unique_cdwds

    codewords.remove(FreeModule(R,rdcd_mtrx.dimensions()[1]).zero()) #Remove the zero codeword for generating weight of nonzero codewords

    if verbose:
        print("---> Isolating unique codewords...DONE!")
        
    if metric == 0:
        
        if verbose:
            print("--> Codewords generated (n={}), computing {} weights...".format(len(codewords)+1,metric_str))
            
        weights = HomogeneousWeightVec([ (word, p, ell, check) for word in codewords ])
        weights = [ x[1] for x in weights ]
        
        min_weight = min(weights)
        
        if verbose:
            print("--> Minimum {} distance found: ".format(metric_str) + str(min_weight))
        
        return len(codewords)+1, min_weight
        
    elif metric == 1:
        
        if verbose:
            print("--> Codewords generated (n={}), computing {} weights...".format(len(codewords)+1,metric_str))
            
        weights = ExtLeeWeightVec([ (word, p, ell, check) for word in codewords ])
        weights = [ x[1] for x in weights ]
        
        min_weight = min(weights)
        
        if verbose:
            print("--> Minimum {} distance found: ".format(metric_str) + str(min_weight))
            
        return len(codewords)+1, min_weight
    
    elif metric == 2:
        
        if verbose:
            print("--> Codewords generated (n={}), computing weights...".format(len(codewords)+1))
            print("---> Computing Homogeneous weights...", end='\r')
            
        hom_weights = HomogeneousWeightVec([ (word, p, ell, check) for word in codewords ])
        hom_weights = [ x[1] for x in hom_weights ]
        
        if verbose:
            print("---> Computing Homogeneous weights...DONE!")
            print("---> Computing Extended Lee weights...", end='\r')
            
        lee_weights = ExtLeeWeightVec([ (word, p, ell, check) for word in codewords ])
        lee_weights = [ x[1] for x in lee_weights ]
        
        if verbose:
            print("---> Computing Extended Lee weights...DONE!")
            
        min_hom_weight = min(hom_weights)
        min_lee_weight = min(lee_weights)
        
        if verbose:
            print("--> Minimum Homogeneous distance: " + str(min_hom_weight) + ".")
            print("--> Minimum Extended Lee distance: " + str(min_lee_weight) + ".")
            
        return len(codewords)+1, min_hom_weight, min_lee_weight
            
    else:
        raise ValueError('metric must be either 0, 1, or 2')
        