#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OptExp (c) University of Manchester 2018

OptExp is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>

Created on Tue Nov 27 16:01:46 2018

@author: pablo
"""
import numpy as np
import pandas as pd
import itertools, re 
from scipy.stats import f as FDist, ncf as ncFDist
from .doebase import doeTemplate, promoterList, plasmidList, read_excel

def evaldes( steps, variants, npromoters, nplasmids, libsize, positional, 
             outfile=None, random=False ):
    """ Generate and evaluate an optimal design of a pathway circuit following the template: 
        1. Vector: 1 to nplasmids
        2. Promoter: 1 to npromoters
        3. Gene: 1 to variants
        4. Terminator + Promoter: None to npromoters, prob(None)= 0.5
        5. Gene: 1 to variants
        ...
        Parameters:
        - steps: number of steps in the pathway
        - variants: number of variants in each step (currently fixed)
        - npromoters: number of promoters in each step
        - nplasmids: number of plasmids in each step
        - libsize: desired library size
        - positional: gene rearrangement is allowed
        - outfile: output the doe into outfile if given
        - random: random DoE instead of optimal
    """

    plasmids = plasmidList(nplasmids)
    promoters = promoterList(npromoters)
    
    tree = []
    genes = {}
    for i in np.arange(steps):
        rid = "r%0d" % (i,)
        tree.append(rid)
        genes[rid] = []
        for j in np.arange(variants):
            gid = "g%0d_%0d" % (i,j)
            genes[rid].append(gid)
            
    doe = doeTemplate( tree, plasmids, promoters, genes, positional )
    if outfile is not None:
        doe.to_excel( outfile, index=False )
        fact, partinfo = read_excel( outfile )
    else:
        fact, partinfo = read_excel( None, doedf=doe )
    try:
        seed = np.random.randint(10000)
        starts = 1
        RMSE = 10
        alpha = 0.05
        factors, fnames, diagnostics = makeDoeOptDes(fact, size=libsize, 
                                                     seed=seed, starts=starts,
                                                     RMSE= RMSE, alpha=alpha,
                                                     random=random )
    except:
        raise Exception("No solution")
    diagnostics['steps'] = steps
    diagnostics['variants'] = variants
    diagnostics['npromoters'] = npromoters
    diagnostics['nplasmids'] = nplasmids
    diagnostics['libsize'] = libsize
    return diagnostics


def makeDoeOptDes(fact, size, seed=None, starts=1040, makeFullFactorial=False, RMSE=1, alpha=0.05, verbose=False, random=False):
    """ Full DoE script: 
        - fact: a dictionary contained the desired design    
    """
    # To Do: full factorial
    import pdb
    pdb.set_trace()
    factors = []
    fnames = []
    npos = 0
    nfact = 0
    for pos in sorted(fact):
        name = fact[pos].component+str(pos)
        if len(fact[pos].levels) > 1:
            nfact += 1
            # Currently only working with categorical
 #            if fact[pos]['component'] != 'gene' and '-' not in fact[pos]['levels']:
#                varType = 'Discrete Numeric'
#                theLevels = [ x for x in range(1, len(fact[pos]['levels'])+1 ) ]
#                factors.append( theLevels ) 
#                fnames.append(name)
#            else:
            # varType = 'Categorical'
            theLevels = [ '"L{}"'.format(x) for x in range(1, len(fact[pos].levels)+1 ) ]
            factors.append(set(theLevels))
            fnames.append(name)
        if fact[pos].positional is not None:
            npos += 1

    if npos >  1:
        # Total possible arrangements in orthogonal latin squares
        # varType = 'Categorical'
        theLevels = ['"L{}"'.format(x) for x in range(1, npos*(npos-1)+1)]
        factors.append( set( theLevels ) )
        fnames.append('pos')
        nfact += 1
    if seed is not None:
        np.random.seed( seed )
    else:
        seed = np.random.randint(100000, size=1)
        np.random.seed( seed )
    initGrid(factors)
    if random:
        # If set, perform a random design instead of a D-optimal design
        M = randExp( factors, n=int(size) )
        J = Deff2(M, factors)
    else:
        if np.product( [len(x) for x in factors] ) < size:
            # raise Exception('Library size is too large!')
            # TO DO: make a full factorial
            M = fullFactorial( factors )
            J = Deff2(M, factors)
            size = M.shape[0]
            ix = np.arange(size)
            np.random.shuffle( ix )
            M = M[ix,:]
        else:
            M, J = CoordExch(factors, n=int(size), runs=2, verb=verbose, mode='coordexch', seed=seed)
    M1 = MapDesign2(factors, M)
    X = mapFactors2( M, factors )
    df = pd.DataFrame(M1, columns=fnames)
    pows = CatPower(X , factors, RMSE=RMSE, alpha=alpha)
    rpvs = RPV(X)
    diagnostics = {'J': J, 'pow': pows, 'rpv': rpvs, 'X': X, 
                   'M': M, 'factors': factors,
                   'M1': M1, 'df': df, 'names': fnames, 'seed': seed}
    return factors, fnames, diagnostics


def Deff(X):
    # D-efficiency
    return (100.0/X.shape[0]) * ( np.linalg.det( np.dot( np.transpose( X ), X ) )**(1.0/X.shape[1]))

def Deff2(M, factors):
    X = mapFactors2(M, factors)
    return (100.0/X.shape[0]) * ( np.linalg.det( np.dot( np.transpose( X ), X ) )**(1.0/X.shape[1]))

def Dopt(X):
    # D-optimality
    return np.linalg.det( np.dot( np.transpose( X ), X ) )

def Dopt2(M, factors):
    # D-optimality
    X = mapFactors2(M, factors)
    return np.linalg.det( np.dot( np.transpose( X ), X ) )

def SE(X):
    # Estimation efficiency
    return np.diag( np.linalg.inv( np.dot( np.transpose( X ), X ) ) )

def RPV(X):
    # Relative prediction variance
    XXi = np.linalg.inv( np.dot( np.transpose( X ), X ) )
    return [np.dot( np.dot( np.transpose( X[i,:] ), XXi), X[i,:]) for i in np.arange(X.shape[0])]

def Contrib(X):
    cn = []
    for i in range(0, X.shape[0]):
        cn.append( Dopt( np.vstack( [X[:i,:], X[(i+1):,:]]  )  ) )
    return cn
def VarAdd(X,xj):
    # Variance of adding/removing one experiment
    return np.dot( np.dot( np.transpose(xj) , np.linalg.inv( np.dot( np.transpose( X ), X) ) ), xj )

def randExp( factors, n ):
    # Generate n random experiments
    V = None
    for levels in factors:
        vnew = np.random.randint(0, len(levels), n)
        if V is None:
            V = vnew
        else:
            V = np.vstack( [V, vnew] )
    if len(V.shape) == 1:
        V = np.expand_dims(V, axis=0)
    return np.transpose( V )


#%%
    


def grid(n, weighted=True):
    """ Provide normalized vectors of n-1 dummy variables 
    Useful for computing the model matrix (X) to 
    use pseudo-orthogonal terms in the n-1 hypercube.
    (Experimental)
    In JMP, grid(3) is multiplied by sqrt(2), grid(4) by
    sqrt(3), which brings back the weight of the number of 
    factors
    """
    from sklearn.preprocessing import normalize
    from sklearn.decomposition import PCA
    
    base = np.eye(n)*2 - 1
    pc = PCA(n-1, whiten=True, random_state=0)
    bt = pc.fit_transform(base)
    W = normalize(bt)
    if weighted:
        W = W*np.sqrt(n-1)
    return W


# Precompute the hypercube grids
gridList = {}
def initGrid(factors):
    global gridList
    vmax = set( [len(x) for x in factors] ) 
    for i in vmax:
        try:
            if i < 2:
                continue
        except:
            continue
        gridList[i] = grid(i)


#%%

def mapFactors( factors, M ):
    # Map a numerical factor into [-1,1] range, create dummy variables for categorical factors
    Mn = np.transpose( [np.ones( M.shape[0] )] )
    for i in np.arange( len(factors) ):
        v = factors[i]
        if type(v) == list:
            if len(set(v)) > 1:
                # Normalize between [-1,+1]
                Vn = (2*(M[:,i] - M[:,i].min())/(M[:,i].max()-M[:,i].min()) - 1)
                Vn = np.transpose( [Vn] )
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        else:
            if len(v) > 1:
                Vn = -np.ones( (M.shape[0],len(v)) )
                j = np.arange(M.shape[0])
                Vn[j,M[j,i]] = 1
                Vn = Vn[:,:-1]
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        if Mn is None:
            Mn = Vn
        else:
            Mn = np.hstack( [Mn, Vn])
    return Mn

def mapFactors2( M, factors ):
    # Map a numerical factor into [-1,1] range, 
    # create orthogonal coordinates for dummy variables for categorical factors
    global gridList
    # Add column for intercept
    Mn = np.transpose( [np.ones( M.shape[0] )] )
    for i in np.arange( len(factors) ):
        v = factors[i]
        Vn = None
        if type(v) == list:
            if len(set(v)) > 1:
                # Normalize between [-1,+1]
                Vn = (2*(M[:,i] - M[:,i].min())/(M[:,i].max()-M[:,i].min()) - 1)
                Vn = np.transpose( [Vn] )
            else:
                Vn = np.transpose( [np.ones( M.shape[0] )] )
        else:
            if len(v) > 1:
                # Use grid 
                j = np.arange(M.shape[0])
                Vn = gridList[len(v)][M[j,i],:]
            else: # 19/02/13: Constant factor, already contained in the intercept
                #Vn = np.transpose( [np.ones( M.shape[0] )] )
                pass
        if Vn is not None:
            if Mn is None:
                # Maybe not needed, we alrady have the intercept
                Mn = Vn
            else:
                Mn = np.hstack( [Mn, Vn])
    return Mn


def MapExp( E ):
    """ Read a design, transform into X matrix """
    # Define factors in the same way as for the library
    factors = [ set(np.unique(E[:,i])) for i in np.arange(E.shape[1])]
    initGrid(factors)
    EE = np.transpose( np.array([ list(np.unique(E[:,i], return_inverse=True)[1]) for i in np.arange(E.shape[1])] ) )
    M = mapFactors( factors, EE )
    return M, factors, EE

def MapDesign(factors, X):
    """ Map back from X to the factors """
    M = []
    for i in np.arange(X.shape[0]):
        row = []
        # skip intercept
        j = 1
        for fa in factors:
            levels = sorted(fa)
            # If none is set
            level = levels[-1]
            for l in levels[0:-1]:
                if X[i,j] == 1:
                    level = l
                j += 1
            row.append(level)
        M.append( row )
    return np.array( M )
    
def MapDesign2(factors, M):
    """ Map back from M to the factors """
    N = []
    for i in np.arange(M.shape[0]):
        row = []
        for j in np.arange(M.shape[1]):
            levels = sorted(factors[j])
            row.append( levels[ M[i,j] ])
        N.append( row )
    return np.array( N )


def JMPRead(jmpfile):
    """ This is a JMP example: """
    # Design Evaluation
    # Design Diagnostics
    # D Optimal Design	
    # D Efficiency	 87.98414
    # G Efficiency	 64.62616
    # A Efficiency	 76.00696
    # Average Variance of Prediction 1.229865
    # Design Creation Time (seconds) 11
    # Read design
    if re.search(jmpfile, 'xlsx$'):
        E = pd.read_excel(jmpfile)
    else:
        E = pd.read_csv(jmpfile)
    # Check which columns to keep
    # We assume that categorical variables are of format 'Li' (Exclude pattern or index columns)
    # Numeric nan factors are not considerd (extra Y column in JMP)
    cols = []
    for i in np.arange(0,E.shape[1]):
        if i == 0:
            if E.columns[i] == 'Pattern':
                continue
        if i == E.shape[1]-1:
            if E. columns[i] == 'Y':
                continue
        if type(E.iloc[0,i]) == str:
            try:    
                assert E.iloc[0,i].startswith('L')
                cols.append(i)
            except:
                continue
        else:
            try:
                assert np.isnan(E.iloc[0,i])
            except:
                cols.append(i)
    labels = E.columns[cols]
    D = np.array(E.iloc[:,cols])
    # Map into a binary matrix 
    DD, fac, EE = MapExp(D)
    # Compute D-efficiency (wrong use of categorical factors)
 #   print( Deff( DD ) )
    # 38.66
    # Compute D-efficiency (correct)
    Deff =  Deff2( EE, fac) 
    # D Efficiency	 87.98414
    return fac, DD, EE, Deff, labels

    
    
#%%

def fullFactorial( factors ):
    # Here we generate a full factorial but this is not possible for large designs
    # Replace by random sampling, ideally some descent algorithm (TBD)
    # For categorical variables, we randomize the levels that are then mapped into the dummy variables
    val = []
    # Add a constant for the intercept
#    val.append( [1] )
    for v in factors:
        val.append( np.arange(len(v)) )
    ffact = []
    for m in itertools.product(*val):
        ffact.append(m)
    ffact = np.array(ffact)
    return ffact

#%%


def DetMax( factors, n, m, it=1000, th=99.5, k=1 ):
    # n: Number of runs
    # m: sampled design space per iteration
    # w: maximum iterations
    # k: number of removal/additions per iteration
    # th: stop criterium for D-efficiency (<= 100.0)
    # Here I have implemented a simple DETMAX algorithm. At each iteration: 
    #  - remove the design with the lowest variance
    # - add the design with the highest variance
    # Many more efficent variants exist (kl-exchange, genetic algorithms, etc..)
    # See Mandal et al, Algorithmic searches for optimal design
    
    # Initial design: could we do something better than a purely random start?   
    M = randExp( factors, n )
    # Map factors into [-1,1] and convert categorical
    X = mapFactors( factors, M )
    # D-Efficiency of the initial design
    J = Deff(X)
    print(J)
    w = 0
    while ((J<100.0) and (w < it)):
        # X1 is the design space sample in the iteration. 
        # Here we loop through the full factorial, which is computationally expensive
        # First thing to fdo is to change it to random generation of a subset of the full library
        # It would be better to move across some surface like gradient descent...
        M1 = randExp( factors, m )
        X1 = mapFactors( factors, M1 ) 
        # Calcualte delta of removing an experiment
        sub = []
        for i in np.arange(X.shape[0]):
            sub.append( VarAdd(X, X[i,:]) )
        w += 1
        # Remove the experiments with the lowest delta
        Xsub = None
        dList = np.argsort( sub )[0:(k+1)]
        for i in np.arange(X.shape[0]):
            if i in dList:
                continue
            else:
                if Xsub is None:
                    Xsub = X[i,:]
                else:
                    Xsub = np.vstack( [Xsub, X[i,:]] )
        # Calculate the delta of adding an experiment from the sample
        add = []
        for j in np.arange(X1.shape[0]):
            add.append( VarAdd( Xsub, X1[j,:] ) )
        # Add the experiments with the highest delta
        aList = np.flip( np.argsort( add ), axis=0 )[0:(k+1)]
        Xn = Xsub
        for j in aList:
            Xn = np.vstack( [Xn, X1[j,:] ] )
        # Make the update if the resulting design increases the objective
        if w % 10 == 0:
            print(w,J,i,j, Dopt(X), Dopt(Xsub), Dopt(Xn))
        if Dopt(Xn) > Dopt(X):
            X = Xn
            J = Deff(X)
        elif Dopt(Xn) == Dopt(X):
            break
    print(w,J)
    return X



def CoordExch1( factors, n, mode='cordexch', verb=False, obj=Dopt, seed=None ): # Deff2
    # Start with an already sub-optimized design by DetMax
    # (it does not make too mauch difference)
    if verb:
        print('Init design')
 #   M = DetMax2( factors, n, 100 )
    M = randExp( factors, n )
    # No optimization (useful for debugging)
    if mode == 'random':
        if verb:
            print('Random')
        eff =  Deff2(M, factors)
        return M, eff
    elif mode == 'detmax':
        if verb:
            print('DetMax')
        M = DetMax2( factors, n, 100 )
        eff =  Deff2(M, factors)
        return M, eff      
    if verb:
        print('Start Coord Exchange Optimization')
    # D-Efficiency of the initial design
    J = 0
    Jn = Dopt2(M, factors)
    q = M.shape[0]
    while J < Jn:
        if verb:
            print(J,Jn) 
        J = Jn
        X = mapFactors2( M, factors ) 
        # Calcualte delta of removing an experiment
        if verb:
            print('Computing variances')
        sub = []
        for i in np.arange(X.shape[0]):
            sub.append( VarAdd(X, X[i,:]) )
        dList = np.argsort( sub )
#        for i in np.arange( M.shape[0] ):
        if verb:
            print('Start exchanging',q,'rows')
        for i in dList[0:q]:
            for j in np.arange( M.shape[1] ):
                Js = []
                for k in np.arange( len(factors[j]) ):
                    M[i,j] = k
                    Jk = Dopt2(M, factors)
                    if np.isnan( Jk ):
                        Jk = 0
                    Js.append( Jk )
                M[i,j] = np.argmax( Js )
        if verb:
            print('End mapping')
        Jn = Dopt2( M, factors )
        if q > 0.2*M.shape[0]:
            q = q-1
    eff =  Deff2(M, factors)
    if verb:        
        print( Deff2(M, factors) )
    return M, eff




def CoordExch( factors, n, mode='cordexch', runs=10, verb=True, seed=None ):
    M = None
    J = 0
    for i in np.arange( runs ):
        Mn, Jn = CoordExch1(factors,n,mode=mode, verb=verb, seed=seed)   
        if Jn > J:
            M = Mn
            J = Jn
    return M, J
    
def DetMax2( factors, n, m, it=1000, th=99.5, k=1, verb=False ):
    # n: Number of runs
    # m: sampled design space per iteration
    # w: maximum iterations
    # k: number of removal/additions per iteration
    # th: stop criterium for D-efficiency (<= 100.0)
    # Here I have implemented a simple DETMAX algorithm. At each iteration: 
    #  - remove the design with the lowest variance
    # - add the design with the highest variance
    # Many more efficent variants exist (kl-exchange, genetic algorithms, etc..)
    # See Mandal et al, Algorithmic searches for optimal design
    
    # Initial design: could we do something better than a purely random start?   
    M = randExp( factors, n )
    X = mapFactors2( M, factors ) 
    # D-Efficiency of the initial design
    J = Deff2(M, factors)
    if verb:
        print(J)
    w = 0
    while ((J<100.0) and (w < it)):
        # X1 is the design space sample in the iteration. 
        # Here we loop through the full factorial, which is computationally expensive
        # First thing to fdo is to change it to random generation of a subset of the full library
        # It would be better to move across some surface like gradient descent...
        M1 = randExp( factors, m )
        X1 = mapFactors2( M1, factors ) 
        # Calcualte delta of removing an experiment
        sub = []
        for i in np.arange(X.shape[0]):
            sub.append( VarAdd(X, X[i,:]) )
        w += 1
        # Remove the experiments with the lowest delta
        Xsub = None
        dList = np.argsort( sub )[0:(k+1)]
        for i in np.arange(X.shape[0]):
            if i in dList:
                continue
            else:
                if Xsub is None:
                    Xsub = X[i,:]
                    Msub = M[i,:]
                else:
                    Xsub = np.vstack( [Xsub, X[i,:]] )
                    Msub = np.vstack( [Msub, M[i,:]] )
        # Calculate the delta of adding an experiment from the sample
        add = []
        for j in np.arange(X1.shape[0]):
            add.append( VarAdd( Xsub, X1[j,:] ) )
        # Add the experiments with the highest delta
        aList = np.flip( np.argsort( add ), axis=0 )[0:(k+1)]
        Xn = Xsub
        Mn = Msub
        for j in aList:
            Xn = np.vstack( [Xn, X1[j,:] ] )
            Mn = np.vstack( [Mn, M1[j,:] ] )
        # Make the update if the resulting design increases the objective
        if w % 10 == 0 and verb:
            print(w,J,i,j, Deff2(M,factors))
        if Dopt(Xn) > Dopt(X):
            X = Xn
            M = Mn
            J = Deff2(M,factors)
        elif Dopt(Xn) == Dopt(X):
            break
    if verb:
        print(w,J)
    return M

#%% Power Analysis test
    

def SimpleCase():
    """ A simple example taken from JMP
    Power should be = 0.299 """
    ft = [{'"L1"', '"L2"', '"L3"', '"L4"'}]
    initGrid(ft)
#    M = randExp(ft, 8 )
#    M, J = CoordExch( ft, 8 )
    M = np.array( [ [0], [1],[2],[3],[0],[1],[2],[3] ])
    M = np.array([[0],[2],[1],[1],[0],[2],[3],[3]])
    X = mapFactors2( M, ft )
    return X, ft

def BigCase(slib,nf, RMSE=1, alpha=0.05):
    ft = []
    for i in np.arange(nf):
        ft.append( set(["L{}".format(x) for x in np.arange(np.random.randint(3,8))]) )
    initGrid(ft)
    M = randExp( ft,slib )
#    M, J = CoordExch( ft, 8 )
    X = mapFactors2( M, ft )
    return CatPower(X , ft, RMSE, alpha)

def Lexc(XX_inv, beta_A, i, varis, nvar, MSE=1):
    """ Compute L excluding the current categorical factor
    Probably this is not the correct one. """
    p_i = int( np.sum(varis[0:i]))
    n_i = varis[i]
    L = np.zeros( (beta_A.shape[0]-n_i, beta_A.shape[0]) )
    offset = 0
    for j in np.arange(0,nvar):
         p_j = int( np.sum(varis[0:j]))
         n_j = varis[j]
         if j > i:
             offset = n_i
         if j != i:
            L[p_j-offset:p_j+n_j-offset,p_j:p_j+n_j] = np.eye(n_j)
    left = np.transpose( np.dot(L,beta_A) )
    try:
        mid = np.linalg.inv( np.dot( np.dot(L, XX_inv),np.transpose(L) ) )/MSE**2
    except:
        import pdb
        pdb.set_trace()
    right = np.dot(L,beta_A) 
    lambda_i = np.dot( left, np.dot(mid, right) )
    return lambda_i

def Linc(XX_inv, beta_A, i, n_i, varis, nvar, RMSE=1):
    """ Compute L by keeping only the current categorical effect of a whole factor.
    Probably this is the correct one. """
    L = np.zeros( (n_i, beta_A.shape[0]) )
    for j in np.arange(0,nvar):
         p_j = int( np.sum(varis[0:j]))
         n_j = varis[j]
         if j == i:
            L[:,p_j:p_j+n_j] = np.eye(n_j)
    left = np.transpose( np.dot(L,beta_A) )
    mid = np.linalg.inv( np.dot( np.dot(L, XX_inv),np.transpose(L) ) )/(n_i*RMSE**2)
    right = np.dot(L,beta_A) 
    lambda_i = np.dot( left, np.dot(mid, right) )
    return lambda_i
   
def CatPower(X, factors, RMSE=1, alpha=0.05):
    """ Power analysis assuming that all
    variables are categorical.
    The returned value is different from the one returned by JMP 12 for SimpleCase().
    However, the result is the same with both the model matrices calculated here 
    or in JMP, so I assume that X is correct.
    The result is also identical to the example in Appendix C of 
    "Power Analysis Tutorial for Experimental Design Software"
     https://apps.dtic.mil/docs/citations/ADA619843
    For SimpleCase() change MSE=1.73
    """

    varis = [1] + [len(x) - 1 for x in factors]
    if np.sum(varis) != X.shape[1]:
        raise Exception("Number of variables do not match")
    n = X.shape[0]
    p = X.shape[1]-1
    nvar = len(varis)
    vec = []
    for i in np.arange(nvar):
        n_i = varis[i]
        vec.append( np.tile( [1,-1], int( np.ceil( n_i/2+1) ) )[0:n_i] )
    beta_A = np.concatenate( vec )
    XX = np.dot( np.transpose(X), X) 
    XX_inv = np.linalg.inv(XX)
    # skip intercept for the time being. To do: add power analysis for continuous parameters
    pows = []
    for i in np.arange(1,nvar):
        p_i = int( np.sum(varis[0:i]))
        n_i = varis[i]
        lambda_i = Linc(XX_inv, beta_A, i, n_i, varis, nvar, RMSE)
        fc_i = FDist.ppf(1-alpha, n_i, n-p-1)
        pow_i = 1-ncFDist.cdf(fc_i, n_i, n-p-1,lambda_i)
        pows.append(pow_i)
    return pows
        
#%% Routines for optimization based on genetic algorithms 
### The optimization so far is not very efficient  (coordinate exchange works better)   

def blending(A,B):
    blend = np.random.randint(2, size=A.shape[0])
    ib = np.argwhere( blend == 1 )
    offspring = A.copy()
    offspring[ib] = B[ib]
    return offspring

def crossover(A,B):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    x1 = A[i,np.arange(0,j)]
    x2 = A[i,np.arange(j,A.shape[1])]
    y1 = B[i,np.arange(0,j)]
    y2 = B[i,np.arange(j,A.shape[1])]
    A[i] = np.append(x1,y2)
    B[i] = np.append(y1,x2)

def mutation(A, th=2.0):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    epsilon = np.random.normal()
    if epsilon > th:
        if A[i,j] > 0:
            A[i,j] = -1
        else:
            A[i,j] = 1
        

def reproduction(A,B):
    A = A.copy()
    B = B.copy()
    for i in np.arange(100):
        if np.random.randint(100) > 95:
            mutation(A)
        if np.random.randint(100) > 95:
            mutation(B)
    # Crossover for multilevel numerical factors won't work
    for i in np.arange(10):
       if np.random.randint(100) > 80:
            crossover(A,B)
    C = blending(A,B)
    return C
    

    
def GenAlg( factors, n, m, nPop=100, it=10, th=99.5, k=1, func=Dopt):
    # Using genetic algorithms
    # Based on Using a Genetic Algorithm to Generate Doptimal Designs for Mixture Experiments
    # To be revised: mutations need to be done in the multilevel design space
    # 1. Generate an initial population of random designs
    # Take the best ones
    population = []
    for i in np.arange(nPop*20):
        M = randExp( factors, n )
        X = mapFactors( factors, M )
        population.append( X )
    population = np.array(population)
    eff = []
    for j in np.arange(nPop):
        eff = np.append( eff, func(population[j]) )
    effi = np.flip( np.argsort(eff), axis=0)
    population = population[effi[0:nPop]]
    # 2. Calculate score for each member of the population
    # and select elite chromosome
    for i in np.arange(it+1):
        eff = np.array( [] )
        for j in np.arange(nPop):
            eff = np.append( eff, func(population[j]) )
        effi = np.flip( np.argsort(eff), axis=0)
        elite = effi[0]
        print(i,Deff(population[elite]))
        if i == it:
            return population[elite]
        population = population[effi[0:nPop]]
        eff = eff[effi[0:nPop]]
#        population = np
        # 3. Random pairing
#        w = np.arange(1, len(eff)-1)
#        np.random.shuffle( eff[w] )
        pairs = []
  #      for i in np.arange(len(w), step=2):
        for i in np.arange(50):
            a = np.random.randint(50)
            b = np.random.randint(20)
            if a != b:
                pairs.append( (population[a], population[b]) )
        for p in pairs:
            offspring = reproduction(p[0],p[1])
            population = np.insert(population,-1,offspring, axis=0)
    
def blending2(A,B):
    blend = np.random.randint(2, size=A.shape[0])
    ib = np.argwhere( blend == 1 )
    offspring = A.copy()
    offspring[ib] = B[ib]
    return offspring

def crossover2(A,B):
    i = np.random.randint(A.shape[0])
    j = np.random.randint(A.shape[1])
    x1 = A[i,np.arange(0,j)]
    x2 = A[i,np.arange(j,A.shape[1])]
    y1 = B[i,np.arange(0,j)]
    y2 = B[i,np.arange(j,A.shape[1])]
    if np.random.randint(2)>0:
        A[i] = np.append(x1,y2)
        B[i] = np.append(y1,x2)
    else:
        A[i] = np.append(y1,x2)
        B[i] = np.append(x1,y2)

def mutation2(A, candidates, th=1.0):
    """ Mutation consists on replacing a run
    with one from the set of candidates """
    epsilon = np.random.normal()
    if epsilon > th:
        i = np.random.randint(A.shape[0])
        j = np.random.randint(len(candidates))
        A[i,:] = candidates[j][i,:] 
        
def reproduction2old(A,B):
    A = A.copy()
    B = B.copy()
    for i in np.arange(100):
        mutation2(A)
        mutation2(B)
    # Crossover for multilevel numerical factors won't work
    for i in np.arange(5):
       if np.random.randint(100) > 50:
            crossover2(A,B)
    C = blending2(A,B)
    return C

def reproduction2(population):
    """ Ranked candidates """
    pairs = []
    for i in np.arange(10):
        pairs.append( (population[0].copy(), population[i+1].copy()) )
    offsprings = []
    for p in pairs:
        A = p[0]
        B = p[1]
        # Crossover for multilevel numerical factors won't work
        for i in np.arange(5):
           if np.random.randint(100) > 50:
                crossover2(A,B)
        C = blending2(A,B)
        for i in np.arange(100):
            mutation2(C, population[0:10])
        offsprings.append( C )
    return offsprings
  
    
def GenAlg2( factors, n, nPop=100, it=10, th=99.5 ):
    # Using genetic algorithms
    # Based on Using a Genetic Algorithm to Generate Doptimal Designs for Mixture Experiments
    # To be revised: mutations need to be done in the multilevel design space
    # 1. Generate an initial population of random designs
    # Take the best ones
    population = []
    for i in np.arange(nPop):
#        M = randExp( factors, n )
        M = DetMax2(factors, n, m=100, it=10, k=2)
        population.append( M )
    population = np.array(population)
    eff = []
    for j in np.arange(nPop):
        eff = np.append( eff, Deff2(population[j], factors) )
    eff[ np.isnan(eff) ] = 0
    effi = np.flip( np.argsort(eff), axis=0)
    population = population[effi[0:nPop]]
    # 2. Calculate score for each member of the population
    # and select elite chromosome
    for i in np.arange(it+1):
        eff = np.array( [] )
        for j in np.arange(len(population)):
            eff = np.append( eff, Deff2(population[j], factors) )
        eff[ np.isnan(eff) ] = 0
        effi = np.flip( np.argsort(eff), axis=0)
        elite = effi[0]
        print(i,Deff2(population[elite],factors),Deff2(population[effi[1]],factors))
        if i == it:
            return population[elite]
        population = population[effi[0:nPop]]
        eff = eff[effi[0:nPop]]
        offsprings = reproduction2(population)         
        population = np.insert(population,-1,offsprings, axis=0)
#%%    

if __name__ == '__main__':
        
    n = 46 # Number of runs
    m = 100 # Sampled design space per iteration
    
    factors = [ [0,1,2,3,4], 
               [123,53,345],#[0,1], [0,1], [0,1]]
       {'Red', 'Green', 'Blue'}, 
       {'prom1', 'prom2', 'prom3', 'prom4'} ]
    
    factors = [
       {'Red', 'Green', 'Blue'}, 
       {'prom1', 'prom2', 'prom3', 'prom4','prom1', 'prom2', 'prom3', 'prom4'},
       {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
       {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
       {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
       {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
       {'Red', 'Green', 'Blue' ,'prom1', 'prom2', 'prom3', 'prom4'}, 
       {'prom1', 'prom2', 'prom3', 'prom4'} ]
    
    factors,DD, EE, Deff, labels = JMPRead('/mnt/SBC1/data/OptimalDesign/data/test3.csv')
    
    
    M , J = CoordExch(factors, n, runs=10)

   # M = GenAlg2(factors, n=46, it=10, nPop=5)
    
    #X = DetMax(factors, n, m, it=1000, k=2)
    
    
    #M = DetMax2(factors, n, m, it=100, k=2)
    
    #best = 0
    #while best < 100:
    ##    M = DetMax2(factors, n, m=100, it=10, k=2)
    #    M = GenAlg2(factors, n=46, it=10)
    #    de = Deff2(M, factors)
    #    
    #    if de > best:
    #        best= de
    #        print(best)
    #
    
    #X = GenAlg(factors, n, m=20, it=100, func=Dopt)
    
    #X = GenAlg2(factors, n=46, it=100)
    
