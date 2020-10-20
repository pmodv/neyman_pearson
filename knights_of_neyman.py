# neyman pearson excursion

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations, chain
from scipy.special import binom
import scipy.stats as ss


# code stolen from web - basis for excursion
# powerset recipe from itertools page
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

 
def isWeaklyDominated(a,x):
    """Returns True if a is weakly dominated by a member of x.
    a : list; x : list of lists."""
    if weaklyDominates(x[0],a):
        return True
    elif len(x) == 1:
        return False
    else:
        return isWeaklyDominated(a,x[1:])
    
def findDiscreteRegions(null,alt):
    """Finds size and powers of rejection regions: rR[i] has size sizes[i] and power powers[i].
    Parameters:
      null, alt : list of probabilities summing to 1, e.g. null = [.1,.7,.2].
    Returns:
      rR, list of tuples : [(), (0), .., (0,4,7), .., (0,..,n-1)]
      sizes, list of floats: [size0, size1, ..]
      powers, list of floats: [power0, power1, ..]"""
    
    # CHECK INPUTS MAKE SENSE
    assert len(null) == len(alt), "Inputs have to be the same length!"
    
    # FORM REGIONS AND WORK OUT SIZES AND POWERS
    listOutcomes = list(range(len(null))) # [0,..,n-1] where n outcomes
    rR = list(powerset(listOutcomes)) # [(), (0,), .., (0,4,7), .., (0,..,n-1)]
    sizes = []; powers = []
    for region in rR:
        # rR[i] has size sizes[i] and power powers[i]
        sizes.append(np.sum(list(map(lambda x: null[x], region)))) # (0,1) --> [.1,.7] --> .8
        powers.append(np.sum(list(map(lambda x: alt[x], region))))

    return rR, sizes, powers

def printDiscreteRegions(null,alt):
    """Prints sizes and powers of rejection regions."""
    
    # get data
    rR, sizes, powers = findDiscreteRegions(null,alt)
    
    # make dataframe and print
    foo = np.array([rR,sizes,powers]).T
    df = pd.DataFrame(foo,columns=['region', 'size', 'power'])
    print(df)
    
def plotDiscreteRegions(null, alt):
    """Plots size and powers of rejection regions"""
    
    # get data
    rR, sizes, powers = findDiscreteRegions(null,alt)
    
    # make the plot
    plt.scatter(sizes, powers, alpha=0.5)
    plt.xlabel('size'); plt.ylabel('power')
    plt.title("""Sizes and powers of all rejection regions.""")
    plt.show()

def weaklyDominates(a,b):
    """Returns True if a[i] >= b[i] for all i and a[j] > b[j] for some j.
    a,b array-like. Must support slicing, e.g. a[1:]"""
    if a[0] < b[0]:
        return False
    elif len(a) == 1:
        return a[0] > b[0]
    elif a[0] > b[0]:
        return veryWeaklyDominates(a[1:],b[1:])
    else: # a[0] == b[0] and length > 1
        return weaklyDominates(a[1:],b[1:])
    
def veryWeaklyDominates(a,b):
    """Returns True if a[i] >= b[i] for all i.
    a,b : lists of equal length."""
    if a[0] < b[0]:
        return False
    elif len(a) == 1:
        return True
    else:
        return veryWeaklyDominates(a[1:],b[1:])



def findDiscreteRegionsPlus(null,alt):
    """Finds data about rejection regions: sizes and powers, dominated regions, likelihood test regions.
    Indexes do book-keeping work: data about rR[i] is in sizes[i], powers[i], domRegions[i], lRegions[i].
    Parameters:
      null, alt : list of probabilities summing to 1, e.g. null = [.1,.7,.2].
    Returns:
      rR, list of tuples : [(), (0), .., (0,4,7), .., (0,..,n-1)]
      sizes, list of floats: [size0, size1, ..]
      powers, list of floats: [power0, power1, ..]
      domRegions, boolean list : jth entry True if jth region is dominated
      lRegions, boolean list: jth entry True if jth region is lrt region"""
    
    # CHECK INPUTS MAKE SENSE
    
    assert len(null) == len(alt), "Inputs have to be the same length!"
    
    # GENERATE REGIONS AND WORK OUT SIZES AND POWERS -- rR, sizes, powers
    
    listOutcomes = list(range(len(null))) # [0,..,n-1] where n outcomes
    rR = list(powerset(listOutcomes)) # [(), (0,), .., (0,4,7), .., (0,..,n-1)]
    sizes = []; powers = []
    
    for region in rR:
        # rR[i] has size sizes[i] and power powers[i]
        sizes.append(np.sum(list(map(lambda x: null[x], region)))) # (0,1) --> [.1,.7] --> .8
        powers.append(np.sum(list(map(lambda x: alt[x], region))))

    # FIND DOMINATED REJECTION REGIONS -- domRegions
    
    negSizesAndPowers = list(zip(list(map(lambda x: 1-x, sizes)),powers)) # [(1-size0,power0), (1-size1,power1),..]
    domRegions = list(map(lambda x: isWeaklyDominated(x, negSizesAndPowers), negSizesAndPowers))
    
    # FIND LIKELIHOOD RATIO TEST REGIONS -- lRegions
     
    likelihoodRatios = np.divide(alt,null)
    outcomesAndRatios = [[i,likelihoodRatios[i]] for i in range(len(likelihoodRatios))]
    outcomesAndRatios.sort(reverse=True,key = lambda x: x[1]) # sort by decreasing likelihood ratio
    
    # form lrt regions by adding in outcomes in turn, in order of decreasing likelihood ratio
    # careful: outcomes with the same likelihood ratio need to be added in jointly
    likelihoodRegions = [[]] # initialize with empty list, since empty region is always lrt region
    workingRegion = []
    for counter, (outcome,ratio) in enumerate(outcomesAndRatios):
        workingRegion.append(outcome)
        try:
            if outcomesAndRatios[counter+1][1] < ratio: # next ratio is strictly less than current ratio
                likelihoodRegions.append(workingRegion[:])
        except: # reached end of list
            likelihoodRegions.append(workingRegion[:])
    
    # generate boolean list from likelihoodRegions
    # to test equality, need to convert regions to sets
    regionsToSets = [set(x) for x in likelihoodRegions] 
    lRegions = [True if set(x) in regionsToSets else False for x in rR]
        
    return rR, sizes, powers, domRegions, lRegions
    
def printDiscreteRegionsPlus(null,alt):
    """Prints data about rejection regions: sizes and powers, dominated regions, likelihood ratio test regions."""
    
    # get data
    rR, sizes, powers, domRegions, lRegions = findDiscreteRegionsPlus(null,alt)
    
    # make dataframe and print
    foo = np.array([rR,sizes,powers,domRegions,lRegions]).T
    df = pd.DataFrame(foo,columns=['region', 'size', 'power', 'dominated?', 'lrt region?'])
    print(df)
        
def plotDiscreteRegionsPlus(null, alt):
    """Plots data about rejection regions: sizes and powers, dominated regions, likelihood ratio test regions."""
    
    # get data
    rR, sizes, powers, domRegions, lRegions = findDiscreteRegionsPlus(null,alt)
    
    # domination
    # generate list of colours to pass to plt.scatter
    dataColours = ['orange' if x else 'blue' for x in domRegions]
    
    # likelihood ratio tests
    # generate list of sizes to pass to plt.scatter
    dataSizes = [100 if x else 30 for x in lRegions]  
    
    # make the plot
    plt.scatter(sizes, powers, c=dataColours, s=dataSizes, alpha=0.5)
    plt.xlabel('size'); plt.ylabel('power')
    plt.title("""Rejection regions. Orange dots are dominated; blue ones aren't.
    Big dots are likelihood ratio tests; small ones aren't.""")
    plt.show()
    
def plotDiscretizedNormalsPlus():
    """Plots sizes and powers of rejection regions testing N(0,1) against N(2,1).
    Divides up R into 11 intervals. Rejection regions are arbitrary unions of the intervals."""
    y1Probs, y2Probs = findDiscretizedNormals(10)
    plotDiscreteRegionsPlus(y1Probs,y2Probs)


plotDiscreteRegionsPlus([.1,.05,.7,.15],[.3,.15,.4,.15])
printDiscreteRegionsPlus([.1,.05,.7,.15],[.3,.15,.4,.15])





