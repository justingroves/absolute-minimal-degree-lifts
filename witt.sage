def WittSum(v1, v2, pols = []):
    
    """
    Calculate the sum of two Witt vectors.
    
    Arguments:
    
        v1 - a Witt vector (namely, a list of elements from the base field the 
            ring is over)
            
        v2 - a Witt vector
        
        pols - See Finotti's notes https://github.com/lrfinotti/witt/blob/master/README.md
        
    Returns:
    
        a Witt vector
    
    """
    F = v1[0].parent()
    p = F.characteristic()
    n = len(v1) - 1
    
    res = [ vRemoveZeros( [v1[i], v2[i] ]) for i in range(n+1)]
    
    if len(pols) == 0:
        pols = etapols(p,n)
        
    if len(pols) > n:
        pols = pols[0:(n-1)]
    
    for i in range(1,n+1):
        l = len(res[i-1])
        
        for j in range(1,l):
            
            temp = vetav(p, n+2-(i+1), [res[i-1][j], sum(res[i-1][k] for k in range(j))],
                         pols = pols)
            
            for t in range(i-1,n):
                if temp[t-i+1] != 0:
                    res[t+1].append(temp[t-i+1])

        pols.pop(n-i)
        
    sum_result = [ sum(x) for x in res ]
    
    #For bug returning Integers instead of elements of the field
    for i in range(len(sum_result)):
        if isinstance(sum_result[i], Integer):
            sum_result[i] = F(sum_result[i])
    
    return sum_result #[ F(sum(x)) for x in res ]
    
def WittProd(v1, v2, pols=[]):
    """
    Calculate the product of two Witt vectors from a ring of Witt vectors 
    over a finite field.
    
    Arguments:
    
        v1 - a Witt vector (namely, a list of elements from the base field the
            ring is over)
            
        v2 - a Witt vector
        
        pols - See Finotti's notes https://github.com/lrfinotti/witt/blob/master/README.md
        
    Returns:
    
        a Witt vector
        
    """
    P = v1[0].parent()
    p = P.characteristic()
    n = len(v1) - 1
    res = [ vRemoveZeros([ v1[j]^(p^(i-j))*v2[i-j]^(p^j) for j in range(0, i+1)]) for i in range(0,n + 1)]

    if n == 1:
        prod_result = [ sum(x for x in res[i]) for i in range(0, len(res))]
        
        #For bug in returning elements as Integers and not field elements
        for i in range(len(prod_result)):
            if isinstance(prod_result[i], Integer):
                prod_result[i] = P(prod_result[i])
                
        return prod_result #[ P(sum(x for x in res[i])) for i in range(0, len(res))]
    if pols == []:
        pols = etapols(p, n-1)
    if len(pols) > (n-1):
        pols = pols[0:n]
    
    for i in range(2, n+1):
        l = len(res[i-1])
        for j in range(1, l):
            temp = vetav(p, n+1-i, [res[i-1][j], sum(x for x in res[i-1][0:j])])
            for t in range(i, n+1):
                if temp[t-i] != 0:
                    res[t].append(temp[t-i])

        del pols[n-i] 
    
    prod_result = [ sum(x for x in res[i]) for i in range(0, len(res))]
    
    #For bug in returning elements as Integers and not field elements
    for i in range(len(prod_result)):
        if isinstance(prod_result[i], Integer):
            prod_result[i] = P(prod_result[i])
    
    return prod_result #[ P(sum(x for x in res[i])) for i in range(0, len(res))]

def WittNeg(v, pols = []):
    
    """
    Calculate the additive inverse of a Witt vector from a ring of Witt vectors 
    over a finite field.
    
    Arguments:
    
        v - a Witt vector (namely, a list of elements from the base field the
            ring is over)
        
        pols - See Finotti's notes https://github.com/lrfinotti/witt/blob/master/README.md
        
    Returns:
    
        a Witt vector
        
    """
    
    p = v[0].parent().characteristic()
    if p != 2:
        return [ -x for x in v]
    
    n = len(v)
    F = v[0].parent()
    vone = [ F(1) for _ in range(n) ]
    
    return WittProd(vone, v, pols = pols)
    
def WittInv(v, pols = []):
    
    """
    Calculate the multiplicative inverse of a Witt vector from a ring of Witt vectors 
    over a finite field.
    
    Arguments:
    
        v - a Witt vector (namely, a list of elements from the base field the
            ring is over)
        
        pols - See Finotti's notes https://github.com/lrfinotti/witt/blob/master/README.md
        
    Returns:
    
        a Witt vector
        
    """
    
    F = v[0].parent()
    p = F.characteristic()
    n = len(v) - 1
    
    if v[0].is_unit() == False:
        raise Exception("Witt vector is not a unit, inverse does not exist")
    
    if len(pols) == 0:
        pols = etapols(p,n-1)
        
    res = [ 1/v[0] ]
    P.<x> = PolynomialRing(F)
    
    for i in range(n):
        w1 = [ P(v[j]) for j in range(i+2) ]
        w2 = [ P(res[j]) for j in range(i+1) ] + [x]
        coord = WittProd(w1,w2, pols = pols)[i+1]
        
        res.append(-coord(0)/(v[0]^(p^(i+1))))
    
    return res
    
def WittToSeries(v, ZqAdic = None):
    
    FF = v[0].parent()
    p = FF.characteristic()
    n = len(v) - 1
    
    #Note the the resulting series will have entries with the generator aa
    if ZqAdic is None:
        k = FF.degree()
        ZqAdic.<aa> = Zq(p^k, prec = n+1)
    
    #Coerce the entries of v into the residue field for the Teichmuller lift
    GG = ZqAdic.residue_field()
    FF_GG_hom = FF.hom(GG.gen(),GG)
    v_mapped = [ FF_GG_hom(v[i]) for i in range(n+1) ]
    
    v_root = [v_mapped[i].nth_root(p^i) for i in range(n+1)]
    v_lift = [ ZqAdic.teichmuller(v_root[i]).mod(p^(n+1-i)) for i in range(n+1)]
    
    return sum( [ v_lift[i]*p^(i) for i in range(n+1)])
    
def SeriesToWitt(s):
    ZqAdic = s.parent()
    n = ZqAdic.precision_cap() - 1
    k = ZqAdic.degree()
    p = ZqAdic.prime()
    
    F = ZqAdic.residue_class_field()
        
    v = []
    
    a = s
    
    for i in range(n+1):
        t = F(a)
        v.append(t)
        
        if t == ZqAdic(0):
            continue
        else:
            a = (a - ZqAdic.teichmuller(t)) // p
            
    return [ v[i]^(p^i) for i in range(n+1) ]
    
def WittPower(v,k, pols=[], bintab=[]):
    F = v[0].parent()
    
    if F.is_finite() and F.is_field():
        p = F.characteristic()
        n = len(v)-1
        d = F.degree()
        v1 = WittToSeries(v)
        tmp = SeriesToWitt(v1^k)
        
        if d > 1:
            G = tmp[0].parent()
            F_G_hom = G.hom(F.gen(),F)
            tmp = [F_G_hom(x) for x in tmp]
            
    else:
        raise Exception("Not implemented for Witt Vectors not over a finite field.")
    
    return [xx for xx in tmp]

def WittFrobenius(v,ell=1):
    F = v[0].parent()
    p = F.characteristic()
    
    if p != 0:
        if F.degree() == 1:
            return v
        else:
            return [ F(xx^(p^ell)) for xx in v]
    else:
        raise ValueError("Base ring must have characteristic non-zero")
        
def WittTrace(v):
    F = v[0].parent()
    
    FF = F.base_ring()
    
    if not FF.is_field():
        raise ValueError("Base ring is not a field")
    
    p = F.characteristic()
    
    if F.degree() == 1:
        return v
    else:
        sum_result = list(v).copy()
        for i in range(1,F.degree()):
            sum_result = WittSum(sum_result,WittFrobenius(v,ell=i))
    
    return [ FF(x) for x in sum_result ]
    