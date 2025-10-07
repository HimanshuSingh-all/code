import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.stats import poisson_binom, binom
from len import BinaryH, finite_len_ab, bool_conv, inv_bool_conv




ALPHA = 0.5
NPTS_1  = 100
NPTS_2  = 200
SPTS = 100

qhmin = 0.5*(ALPHA*0.5 + 0.5 *(1-ALPHA))
_n0 = 110
_n = np.ceil(_n0/(1- ALPHA) )
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

plt.figure(1)
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
_k: int = np.ceil( _n * (1 - finiapprox.BinaryH(a = PICKEDVALS['qpick'], n = _n) )) 
assert isinstance(_k, int)
NUMSTR = 2**_k


plt.show()
