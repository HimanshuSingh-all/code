import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import poisson_binom, binom
from len import BinaryH, finite_len_ab, bool_conv, inv_bool_conv
from len import compute_dist_to_target, compute_l_r_substr_dist, compute_death_rate



CHOICE_P = 0.5
ALPHA = 0.5
NPTS_1  = 20
NPTS_2  = 20
SPTS = 100

qhmin = 0.5*(ALPHA*0.5 + 0.5 *(1-ALPHA))
_n0 = 70 # 110
_n = int(np.ceil(_n0/(1- ALPHA) ))
assert np.ceil(_n*ALPHA)%2 == 0
_p = 0

qh = np.linspace(qhmin, 0.5, NPTS_1)
r  = np.linspace(0, 0.5, NPTS_2)

finiapprox = finite_len_ab()

lmin = np.zeros( (NPTS_2, NPTS_1) )
smin = np.zeros( (NPTS_2, NPTS_1) )

#Grid Search
for i, rr in enumerate(tqdm(r)):
    rh = inv_bool_conv(_p, rr)
    for j, qq in enumerate(qh):
        q_tmp = bool_conv(qq, _p )
        smax = min( 2*rh, 1 - (1- 2*qq)/ALPHA)
        s = np.linspace( 0, smax, SPTS)
        minml = np.inf
        minms = None
        for ss in s:
            l = finiapprox.len(
                q = q_tmp,
                n0= _n0,
                p = _p,
                r = rr,
                s = ss,
                alpha=ALPHA
            )
            if l<minml:
                minml = l
                minms = ss
        smin[i, j] = minms
        lmin[i, j] = np.squeeze(minml)

print(f'Last minml {minml}')
# Getting the masked array where l_min <=0
Q, R  = np.meshgrid(qh ,r)
# Masking the region with less than zero genome
mask = lmin <=0
Z = np.ma.masked_array(lmin, mask)

lowest_inds = np.any(Z.mask, axis=1)
rmin = r[lowest_inds][0]
qind = Z.mask[lowest_inds,:][0]
qpick = qh[qind][0]
spick = smin[lowest_inds,:][0]
spick = spick[qind][0]
PICKEDVALS: dict = {
        'rpick': rmin,
        'spick':spick,
        'qpick':qpick
}

print("The picked value of s: ",spick)
print("The picked value of q: ",qpick)
print("The picked value of r: ",rmin)

#plt.figure(1)
#Plotting the regions
fig, axs = plt.subplots(1, 2 , constrained_layout = True)
map = axs[1].imshow(smin, origin = 'lower', extent=[qh[0], qh[-1], r[0], r[-1]])#, aspect='square')
axs[1].set_title(r'$s_{min}$')
plt.colorbar(mappable= map, ax = axs[1],  shrink = 0.6)
map = axs[0].imshow(Z, origin = 'lower', extent=[qh[0], qh[-1], r[0], r[-1]] ) #, aspect='square')
plt.colorbar(mappable= map, ax = axs[0], shrink = 0.6)
# axs[0].imshow(mask, cmap = 'grey', origin = 'lower', extent=[q[0], q[-1], r[0], r[-1]], aspect = 'equal', alpha = 0.5)

axs[0].contour(Q,R, lmin, colors = ['white'], levels = 20)
axs[0].annotate(r'$\hat{q} = $' + f'{qpick:.3f}'
                +'\n'+ r'$r = $' + f'{rmin:.3f},'
                +'\n'+ r'$s= $' + f'{spick:.3f}',
                xy = (qpick, rmin),
                xytext = (0.01, 0.01),
                arrowprops=dict(arrowstyle='->,head_width=.15',color='red')
                )
axs[0].set_title(r'$l_{min}$')
axs[0].set_xlabel(r'$\hat{q}\to$')
axs[0].set_ylabel(r'$r\to$')

axs[1].set_xlabel(r'$\hat{q}\to$')
axs[1].set_ylabel(r'$r\to$')
#fig.suptitle(r"$p_{repl} = $"+ f"{p_repl}, p = {p}, "+ r"$\alpha= $" + f'{ALPHA}')



# Generate random strings
_k: int = int(np.ceil( _n * (1 - finiapprox.BinaryH(a = PICKEDVALS['qpick'], n = _n) )) )
assert isinstance(_k, int)
fig2, axs2 = plt.subplots(2, 2, constrained_layout = True)
print('Plotting the scatterplot')
for i , k in tqdm(enumerate([ 10 , 13, 16, _k])):
    NUMSTR = 2**k
    NUMTARGS = 100 # number of different target to be constructed
    NUMCONSTR = 100 

    target_strings = np.random.random( size = (1, _n-_n0) ) < CHOICE_P
    fp_strings = np.random.random( size = (NUMSTR, _n) ) < CHOICE_P
    death_rate = np.zeros(shape = (NUMTARGS, ) )
    assert target_strings.ndim==2 

    target = target_strings[0, :] 
    STEP = _n - _n0
    s_arr = compute_l_r_substr_dist(strs = fp_strings, alpha= ALPHA, disttype='relative')
    arr_1 = compute_dist_to_target(target_str= target, strs = fp_strings[:, :STEP])
    arr_2 = compute_dist_to_target(target_str= target, strs = fp_strings[:, STEP:2*STEP])
    r_arr = np.where(
            arr_1>arr_2,
            arr_1,
            arr_2
    )

    #plt.figure(2)
    axs2[i//2, i%2].scatter( s_arr, r_arr, s = 2)
    axs2[i//2, i%2].set_title(f'k = {k}, # strings = {2**k}')




NUMSTR = 2**_k
NUMTARGS = 1 # number of different target to be constructed
NUMCONSTR = 100 

target_strings = np.random.random( size = (1, _n-_n0) ) < CHOICE_P
fp_strings = np.random.random( size = (NUMSTR, _n) ) < CHOICE_P
death_rate = np.zeros(shape = (NUMTARGS, ) )
assert target_strings.ndim==2

target = target_strings[0, :] 
STEP = _n - _n0
s_arr = compute_l_r_substr_dist(strs = fp_strings, alpha= ALPHA, disttype='relative')
arr_1 = compute_dist_to_target(target_str= target, strs = fp_strings[:, :STEP])
arr_2 = compute_dist_to_target(target_str= target, strs = fp_strings[:, STEP:2*STEP])
r_arr = np.where(
        arr_1>arr_2,
        arr_1,
        arr_2
    )
#print('here')
#S, R = np.meshgrid(s_arr, r_arr)
#print('here')
# I do not need kdtrees as the _n is an integer which is fixed and the 
# relative error is e/_n or e/_n0, where e = {0, 1, ... _n} or e = {0, 1, ... _n0}
# so it is very likely the 2d plot will have a "lattice" like structure in the dense points
# or the bulk of the scatterplot
bound_inds = np.logical_and(s_arr<=0.5, r_arr<=0.5)
s_unique, counts_s= np.unique(s_arr[bound_inds], return_counts=True)
r_unique, counts_r= np.unique(r_arr[bound_inds], return_counts=True)

print(s_unique.shape)
print(r_unique.shape)


death_rates = np.zeros((r_unique.shape[0], s_unique.shape[0]))
death_ratesmap = np.zeros_like(death_rates)

s_unique = np.sort(s_unique)
r_unique = np.sort(r_unique)
r_ind_elems = list(range(r_unique.shape[0]))
s_ind_elems = list(range(s_unique.shape[0]))
#r_val_ind:dict = dict(zip(r_unique, r_ind_elems))
#s_val_ind:dict = dict(zip(s_unique, s_ind_elems))
#print(MOTHER_STRINGS.shape)
for i, r in tqdm(enumerate(r_unique)):
    inds_r = r_arr==r
    for j, s in enumerate(s_unique):
        inds_s = s_arr==s
        ind = np.logical_and(inds_r, inds_s)
        if not(np.any(ind)):
            death_rates[i,j] = np.nan
            continue
        del inds_s
        mother = fp_strings[ind][0]
        assert mother.shape[0] == _n
        death_rate = compute_death_rate(
            mother = mother,
            fps = fp_strings,
            alpha = ALPHA,
            trials = NUMCONSTR
        )
        death_rates[i, j] +=death_rate



#for index ,mother in tqdm(zip(MOTHER_STR_INDICES, MOTHER_STRINGS)):
#    rval = r_arr[index]
#    sval = s_arr[index]
#    store_ind = (r_val_ind[rval], s_val_ind[sval])
#    death_rate = compute_death_rate(
#        mother = mother,
#        fps = fp_strings,
#        alpha = ALPHA,
#        trials = NUMCONSTR
#    )
#    death_rates[store_ind] +=death_rate
#    death_ratesmap[store_ind]+=1
#
death_rates = death_rates/death_ratesmap
plt.figure(3)
plt.imshow(death_rates, origin = 'lower', extent=[s_unique[0], s_unique[-1], r_unique[0], r_unique[-1]]) 
plt.colorbar()
print('Last Elem')
plt.show()
