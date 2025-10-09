import numpy as np
import matplotlib.pyplot as plt
import warnings
from tqdm import tqdm
from numpy.typing import ArrayLike, NDArray
def BinaryH(a)->float:
    """
    Compute the binary entropy
    """
    a = np.asarray(a, dtype = float)
    return np.where(np.logical_or(a==0, a==1), 0, a*np.log2(1/a) + (1-a)*np.log2(1/(1-a)) )

def H(p:np.ndarray)->float:
    """

    """
    # We'll keep the negative elemets

    inds = np.logical_not(p==0)
    # print(inds.shape)
    # inds = inds.astype(np.bool_)
    # print(inds, '\n', p[inds])
    if np.any(p[inds]<0):
        # print('in', p[inds])
        return np.nan
    return -np.sum(p[inds]*np.log2(p[inds]))

def bool_conv(a,b):
    """ Returns the boolean convolution of a and b i.e a*(1-b) + b*(1-ab)"""
    return a*(1-b) + b*(1-a)

def rate_dist_weight(p):
    """Returns the rate distortion weight for unit string, p is the relative Hamming distance"""
    return 1-BinaryH(p)

def len_gen_ab( p, qhat, alpha, rhat, s):
    """ Returns the length of the genome given that the length of the fixed part n0 ==1 """
    def A(rhat, s):
        r = rhat
        probs = np.array([1-r-s/2, s/2, s/2, r -s/2])
        del r
        return 2-H(probs)
    l = (1/(1-alpha))* (A(rhat, s)*alpha + (1- 2*alpha)*rate_dist_weight(rhat)- rate_dist_weight(bool_conv(p, bool_conv(qhat,p))) )
    return l

def inv_bool_conv(b, c):
    """
        Returns a such that bool_conv(a,b) = c
    """
    return (1/2) * ( (c - b)/((1/2) - b) )


def compute_l_r_substr_dist( strs:NDArray[np.bool], alpha:float, disttype: str = 'relative'):
    """
    Compute the distance between the left alpha*n substring and the substring of the same length contigously next to it for each string in the array. Each row is assumed to be a boolean string.
    ## Inputs:
    - strs (NDArray[np.bool]): An array of bits with each row representing a string of length n, which is the number of columns.
    - alpha (float): The length fraction of the substring compared to the length of the total string. Cannot be `>0.5`.
    - disttype (str): Whether to compute the absolute hamming distance or the relative hamming distance. Options: `{'absolute', 'relative'}`.
    ## Returns:
    - s_arr (NDArray[np.bool]): An
    """
    assert alpha<=0.5, "Alpha cannot be greater than 0.5"
    #assert np.isdtype()
    assert strs.ndim ==2 
    n = strs.shape[1]
    STEP = int(alpha*n)
    cmp_left_right = np.logical_xor(strs[:, :STEP], strs[:, STEP: 2*STEP])
    assert cmp_left_right.shape == (strs.shape[0], STEP)
    if disttype == 'relative':
        s_arr = np.sum( cmp_left_right, axis = 1) /(alpha*n)
    elif disttype=='absolute':
        s_arr = np.sum( cmp_left_right, axis = 1)
    print(f's_arr shape {s_arr.shape}')
    return s_arr

def compute_dist_to_target(target_str:  NDArray[np.bool], strs:NDArray[np.bool], disttype: str = 'relative'):
    """
    Return a distance array for an array of strings to a given target string `target_str`.
    # Inputs:
    - target_str(NDArray[np.bool]): The target string
    - strs (NDArray[np.bool]): An array of bits with each row representing a string of length n, which is the number of columns.
    - alpha (float): The length fraction of the substring compared to the length of the total string. Cannot be `>0.5`.
    - disttype (str): Whether to compute the absolute hamming distance or the relative hamming distance. Options: `{'absolute', 'relative'}`.
    """
    assert target_str.ndim ==1
    assert strs.ndim ==2
    assert strs.shape[1] == target_str.shape[0]
    if disttype == 'relative':
        target_dist = np.sum( np.logical_xor(target_str, strs), axis= 1)/target_str.shape[0]
    elif disttype == 'absolute':
        target_dist = np.sum( np.logical_xor(target_str, strs), axis= 1)
    return target_dist


def compute_death_rate(
        mother:NDArray,
        fps:NDArray,
        alpha:float,
        trials:int,
        p:float =0,
        choice_p: float =0.5
): 
    """ 
    Compute the death rate of the given a mother string.
    """
    assert mother.ndim == 1
    assert fps.ndim == 2
    assert 0<alpha <=0.5
    assert isinstance(trials, int)
    n_ = mother.shape[0]
    STEP = int(np.ceil(alpha*n_) )
    substr_1 = mother[:STEP]
    substr_2 = mother[STEP:2*STEP]
    #toehold = np.where(
    #    np.random.random(np.shape(STEP, ))>np.random.random(np.shape(STEP, )),
    #    substr_1,
    #    substr_2
    #)
    toehold = np.hstack((substr_1[:STEP//2], substr_2[STEP//2:]))
    assert toehold.shape[0] == STEP
    death_rate = 0.0
    last_str = None
    dist_mother = np.zeros(trials)
    dist_minfp = np.zeros(trials)

    for i in range(trials):
        random_part = np.random.random(size = n_-STEP) <choice_p
        #print(f'Toehold shape {toehold.shape}')
        init_str = np.hstack( (toehold, random_part))
        parent_toehold_dist = np.logical_xor( init_str, mother)
        parent_toehold_dist = np.sum(parent_toehold_dist)/n_
        minm = np.inf
        for a_str in fps:
            assert a_str.shape == init_str.shape
            dist = np.sum(np.logical_xor(init_str, a_str))/n_
            minm = min(dist, minm)
        dist_mother[i] = parent_toehold_dist
        dist_minfp[i] = minm
        if minm<parent_toehold_dist:
            death_rate= death_rate + 1
        #if last_str is not None:
        #    print(f'islat = {np.all(last_str == init_str)}')
        #last_str = init_str
    death_rate = death_rate/trials
    return death_rate



class infinite_len_ab:

    def __init__(self) -> None:
        print('Initialised the infinite length class')

    def BinaryH(self, a:float):
        """ Compute the actual binary entropy"""
        # warnings.filterwarnings('default')
        return BinaryH(a)

    def H_n(self, p):
        """ Compute the actual entropy given a list of probabilities"""
        ps = p.T
        ## Number of columns should be 4
        # warnings.filterwarnings('error')
        assert ps.shape[1] == 4
        return np.array([H(ps[i,:]) for i in range(ps.shape[0])])

    def n(self, n0, alpha):
        return n0/(1-alpha)

    def len(self, q, n0, p, r, s, alpha):
        """ 
        Compute the length for the infinte length 
        string genome information.
        ## Inputs:
        - q(float): The size of the basin.
        - n0(int): The fixed part of the string.
        - p (float): The developmental noise.
        - r(float): The "functional" distance "plus" the developmental noise. bool_conv(p,rhat).
        - s(float): The distance between the replicated entities (replicated substrings).
        - alpha(float): The initial part of the string as the intitial condition. 
        """
        # print( f's=  {s}', f'a= {alpha}, r = {r}')
        n = self.n(n0 = n0, alpha=alpha)
        C= 1- self.BinaryH(q)
        rh =inv_bool_conv( p, r)
        B= 1- self.BinaryH(rh)
        probs = np.array([1-rh-s/2, s/2, s/2, rh -s/2])
        A= 2 - self.H_n(probs)
        return n*( alpha*A + (1- 2*alpha) * B - C)


class finite_len_ab:
    """ 
    Finite length approximation wrapper class for the duplication AB model
    has no data attributes just class for clubbing together the functions of the 
    same type.
    """
    def __init__(self) -> None:
        print('Initialised the finite length class')

    def BinaryH(self, a:float, n:int):
        """
        'Binary Entropy modification' for the finite length case

        ## Inputs:
        - a: Probability of one of the events such that p = [a, 1-a].
        - n: Length of the string.
        Note: This migh sometimes throw a DividebyZero warning but that's because \
            the np.where in `BinaryH` function does not perform boolean short circuiting so it will evaluate `a==0` or\
            `a==1` case leading to warnings, but it's still a good practice to know about the 
            warnings so I have left them in.
        """
        return BinaryH(a) + 1/n - (1/(2*n)) * np.log2(2*np.pi*n)

    def H_n(self, p:ArrayLike, n:int ,t:int=4):
        """
        General entropy modification for the case of finite strings.

        ## Inputs:
        - p (ArrayLike): The probability array such that `sum(p)==1` and `np.all(p>0) ==True`.
        - n (int): The length of the string.
        - t (int): Length of the `p`. By default it is 4.
        """
        
        if p.ndim ==1:
            ps = p[:, np.newaxis]
            ps = ps.T
        else:
            ps = p.T
        assert ps.shape[1] == t
        entrop_actual = np.array([H(ps[i,:]) for i in range(ps.shape[0])])
        additive =  ((1/2-t/2)*np.log2((2*np.pi*n)) + (t/2)*np.log2(t))/n
        return entrop_actual + additive

    def n(self, n0, alpha):
        """
        Full length of the string given the length of the fixed part
        and alpha.
        """
        return n0/(1-alpha)

    def len(self, q, n0, p, r, s, alpha):
        """ 
        Compute the length for the infinte length 
        string genome information.
        ## Inputs:
        - q(float): The size of the basin.
        - n0(int): The fixed part of the string.
        - p (float): The developmental noise.
        - r(float): The "functional" distance "plus" the developmental noise. bool_conv(p,rhat).
        - s(float): The distance between the replicated entities (replicated substrings).
        - alpha(float): The initial part of the string as the intitial condition. 
        """ 
        n = self.n(n0 = n0, alpha=alpha)
        assert n>0
        c_t= 1- self.H_n(np.array((q, 1-q)), n, t=2)
        rh =inv_bool_conv(p, r)
        b_t= 1- self.H_n(np.array((rh, 1-rh)), n*(1-2*alpha), t = 2) if n*(1-2*alpha) >0 else 0
        probs = np.array([1-rh-s/2, s/2, s/2, rh -s/2])
        A= 2-self.H_n(probs, alpha*n, t = 4) if n*alpha>0 else 0
        return n*( alpha*A + (1- 2*alpha) * b_t - c_t)

    def A(self, rh, s, n):
        probs = np.array([1-rh-s/2, s/2, s/2, rh -s/2])
        A= 2 - self.H_n(probs, n)
        return A

