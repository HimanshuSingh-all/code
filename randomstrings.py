import numpy as np
import matplotlib.pyplot as plt
import argparse


def BinH(p):
    return 0 if p==0 or p==1 else p*np.log2(1/p) + (1-p)*np.log2(1/(1-p))
n0 = 64
alpha = 0.5
n = int(n0/(1-alpha))
q = 0.25
k = int(np.ceil(n*(1 - BinH(q))))
num_fp = 2**k
print(f'k = {k}')
print(f'n = {n}')
print(f'2**k = {2**k}')
draw = 0.5
X = np.random.choice(a= [0, 1], p=[draw, 1- draw], size = (2**k, n) )

Y_target = np.random.choice(a = [0, 1], size = (n0,) )

absolute_dist_lr = np.zeros((2**k,))
absolute_dist_target_leg1 = np.zeros((2**k,)) 
absolute_dist_target_leg2 = np.zeros((2**k,)) 

lr = np.logical_xor(X[:,:n//2], X[:, n//2:])
print(lr.shape)
lr = np.sum(lr, axis = 1)
print(lr.shape)
absolute_dist_lr = lr
ycpy = np.array([Y_target]*num_fp)
absolute_dist_target_leg1[:num_fp] = np.sum(np.logical_xor(ycpy, X[:, :n//2]), axis =1)
absolute_dist_target_leg2[:num_fp] = np.sum(np.logical_xor(ycpy, X[:, n//2:]), axis =1)

plt.scatter(absolute_dist_lr/(n*alpha), absolute_dist_target_leg1/(n*alpha), color = 'blue')
plt.scatter(absolute_dist_lr/(alpha*n), absolute_dist_target_leg2/(n*alpha), color = 'blue')
plt.xlabel(r"s")
plt.ylabel(r"$\hat{r}$")

plt.show()
