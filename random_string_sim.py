import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
from len import bool_conv, inv_bool_conv, BinaryH, finite_len_ab, H

finiteapprox = finite_len_ab()
# Defining the constants
n0_ = 194//2
s_min = 0.402
r_min = 0.3679
qhat_pick_ = 0.350 # 0.33568
ALPHA = 0.5
p_ = 0.00

# Other things
q_ = bool_conv(p_, bool_conv(qhat_pick_,p_))
n = int(n0_/(1-ALPHA))
kp=  (1- BinaryH(q_))
print(f' kp is {kp};')
k = np.ceil(n*(1- finiteapprox.BinaryH(q_, n)))
print(f'k is {k}')
num_str = int(2**k) #* 16 
print(f'the number of strings are {num_str}')
choice_p = 0.5 # the probability of generating a 1 from the source if it is generated from the binomial




# Generating the random fixed points strings
fp_strings = np.random.choice(a = [False, True], size = (num_str, n), p = (1- choice_p, choice_p)) 
print(f' The shape of the string matrix {fp_strings.shape}')
print(f' The memory occupied of the string matrix {fp_strings.nbytes} bytes')
target_string = np.random.choice( a = [ False, True], size = (n0_,), p = [1 -choice_p, choice_p] )



STEP = int(ALPHA*n)
cmp_left_right = np.logical_xor(fp_strings[:, :int(ALPHA*n)], fp_strings[:, int(ALPHA*n): int(ALPHA*2*n)])
assert cmp_left_right.shape[0] == num_str
print(f'cmp left right shape {cmp_left_right.shape}')
s_arr = np.sum( cmp_left_right, axis = 1) /(ALPHA*n)
print(f's_arr shape {s_arr.shape}')

target_dist_1 = np.sum( np.logical_xor(target_string, fp_strings[:,:STEP]), axis= 1)
target_dist_2 = np.sum( np.logical_xor(target_string, fp_strings[:,STEP:2*STEP]), axis= 1)
print(target_dist_1.shape)
print(target_dist_2.shape)



rhat_arr = np.where( target_dist_1 < target_dist_2, target_dist_1, target_dist_2)/(STEP)
print(rhat_arr.shape)
del target_dist_1, target_dist_2, cmp_left_right



#
def root(rh, s, k, n):
    probs = np.array([1 - rh - s/2, s/2, s/2 , rh -s/2])
    # print('ps:',probs)
    return (k/n) + finiteapprox.H_n(probs,n, 4) - 2
    #return (k/n) + H(probs) - 2
print( 'The extra thing: ',  ( (1/2-4/2)*np.log2(2*np.pi*ALPHA*n) +1/(ALPHA*n) )/n )
s_plt = list()
r_plt = list()

smax = 1 - (1-2*qhat_pick_)/ALPHA
s_check = np.linspace(np.min(s_arr), smax , 100)
for ss in (s_check):
    #args = (ss, k, ALPHA, n)
    args = (ss, k, n0_)
    #print(f'Sol args: {args}')
    fun = lambda x:root(x, *args)
    try:
        sol =root_scalar(fun, bracket= (ss/2, 0.5))
        if sol.converged:
            s_plt.append(ss)
            r_plt.append(sol.root)
    except ValueError:
        continue
#arg_inds = np.argsort(s_plt)
s_plt = np.array(s_plt)#[arg_inds]
r_plt = np.array(r_plt)#[arg_inds]
#
# tolerance: Just doing a basic 1D grid search if the root(rh) -0 = tol 
# we say rh is solution
#r = np.linspace(rhat_arr.min(), rhat_arr.max(), 100)
#tol = 0.0001 
#smax = 1 - (1-2*qhat_pick_)/ALPHA
#r_check = np.linspace(0, 0.5, 100)
#s_check = np.linspace(np.min(s_arr), smax , 100)
#r_plt = list()
#s_plt = list()
#
#for ss in (s_check):
#    args = (ss, k, ALPHA, n)
#    rcheck  = list()
#    for rr in r_check:
#        root_check = root(rr, *args)
#        cond = np.abs(root_check - 0) <tol
#        if cond:
#            #s_plt.append(ss)
#            rcheck.append((rr, root_check))
#            #print(f'Broken, {ss} {rr}')
#    if len(rcheck) == 0:
#        continue
#    else:
#        s_plt.append(ss)
#        rmin = sorted(rcheck, key = lambda x: np.abs(x[1]) )[0]
#        r_plt.append(rmin[0])
#
#
##assert len(s_plt) <=3
#print(f' lens = {len(s_plt)}' )
#s_plt = np.array(s_plt)
#r_plt = np.array(r_plt) 
#
# Plot 
ind = s_arr <= s_min
size = 2
ind2 = np.logical_and(ind, rhat_arr<=r_min)

plt.figure(figsize=(6,6))
plt.scatter( s_arr, rhat_arr, s = size, marker='x')
plt.scatter( s_arr[ind], rhat_arr[ind], s = size, marker='x',color = 'orange')
plt.scatter( s_arr[ind2], rhat_arr[ind2], s = size, color = 'red', marker='x')
plt.plot(s_plt, r_plt, color = 'red')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()


#
#

#actualminr = np.zeros((num_str,))
#ind_min = np.zeros( (num_str,), dtype=int)

#for i in range(num_str):
#    ind = None 
#    mindist = rhat_arr[i]
#    for j in range(num_str):
#        if j==i:
#            continue
#    dist1 =  np.logical_xor( fp_strings[i, :STEP], fp_strings[j, :STEP]) 
#    dist2 =  np.logical_xor( fp_strings[i, :STEP], fp_strings[j, STEP:2*STEP]) 
#
