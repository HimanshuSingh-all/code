import numpy as np
import matplotlib.pyplot as plt

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


class infinite_len_ab:

    def __init__(self) -> None:
        print('Initialised the infinite length class')

    def BinaryH(self, a:float):
        """ Compute the actual binary entropy"""
        return BinaryH(a)

    def H_n(self, p):
        """ Compute the actual entropy given a list of probabilities"""
        ps = p.T
        ## Number of columns should be 4
        assert ps.shape[1] == 4
        return np.array([H(ps[i,:]) for i in range(ps.shape[0])])

    def n(self, n0, alpha):
        return n0/(1-alpha)

    def len(self, q, n0, p, r, s, alpha):
        """ 
        Compute the length for the infinte length 
        string genome information
        """
        # print( f's=  {s}', f'a= {alpha}, r = {r}')
        n = self.n(n0 = n0, alpha=alpha)
        C= 1- self.BinaryH(q)
        rh =inv_bool_conv(p, r)
        B= 1- self.BinaryH(rh)
        probs = np.array([1-rh-s/2, s/2, s/2, rh -s/2])
        A= 2 - self.H_n(probs)
        return n*( alpha*A + (1- 2*alpha) * B - C)


class finite_len_ab:

    def __init__(self) -> None:
        print('Initialised the finite length class')

    def BinaryH(self, a:float, n:int):
        return BinaryH(a) - (1/(2*n)) * np.log2(2*np.pi*n) 

    def H_n(self, p, n ,t=4):
        ps = p.T
        assert ps.shape[1] == 4
        entrop_actual = np.array([H(ps[i,:]) for i in range(ps.shape[0])])
        additive =  ((1/2-t/2)*np.log2((2*np.pi*n)) + np.log2(t/2))/n
        return entrop_actual - additive

    def n(self, n0, alpha):
        return n0/(1-alpha)

    def len(self, q, n0, p, r, s, alpha):
        n = self.n(n0 = n0, alpha=alpha)
        C= (1- self.BinaryH(q, n))
        rh =inv_bool_conv(p, r)
        B= (1- self.BinaryH(rh, n))
        probs = np.array([1-rh-s/2, s/2, s/2, rh -s/2])
        A= (2 - self.H_n(probs, n, t = 4) ) 
        return n*( alpha*A + (1- 2*alpha) * B - C) 




if __name__ == "__main__":
    p = 0.0
    NUMPTS = 1000
    
    alpha = 0.5 #:1 - 2*qh + 0.001
    qmin = 1/4
    qmax = 3*alpha/8
    qhat = [0.25, 0.3, 0.35]
    ls = list()

    s = np.linspace(0, 1/2, NUMPTS)
    rhat = np.linspace(0, 1/2, NUMPTS)
    for qh in qhat:
        l = np.zeros((s.shape[0], rhat.shape[0]))
        for i in range(s.shape[0]):
            for j in range(rhat.shape[0]):
                l[i,j] = len_gen_ab(p = p, qhat= qh, alpha= alpha, rhat= rhat[j], s =s[i] )
        ls.append(l)

    xlabel = r"$\hat{r}\to$"
    ylabel = r"$s\to$"
    # Colormap
    cmap = "viridis"
    # Create figure with 3 subplots (columns)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    x_vals = rhat
    y_vals = s
    for i, ax in enumerate(axes):
        im = ax.imshow(ls[i], cmap=cmap, origin='lower', aspect='auto' ,extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]])
        ax.set_aspect('equal')
        ax.set_title(r"$\hat{q}:$" + f"{qhat[i]}"+ f' p = {p}')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        X, Y = np.meshgrid(x_vals, y_vals)
        mask = Y==2*X 
        # Plot filled contours
        cont = ax.contour(X, Y, ls[i], levels=20, cmap="viridis")
        # Optional: plot contour lines
        ax.plot(X[mask], Y[mask], color='red', label = r'$s=2\hat{r}$')
        zero_mask = np.isclose(ls[i],0)
        print(zero_mask.shape)
        ax.scatter(X[zero_mask], Y[zero_mask], color='white') #, label = r'$s=2\hat{r}$')

        nan_mask = np.isnan(ls[i])
        ax.scatter(X[nan_mask], Y[nan_mask], color='orange') #, label = r'$s=2\hat{r}$')
        ax.plot( rhat, rhat*(1-rhat)*2, color='peru', label = r'$\hat{r} \bigoplus\hat{r}$' )
        print(alpha)
        smin = 1 - (1-2*qhat[i])/alpha
        ax.axhline(y = smin, color  = 'black', linestyle = 'dashed')
        lcpy = np.where(nan_mask, np.inf, ls[i])
        idx = np.argmin(lcpy, axis = 0)
        ax.scatter(rhat,s[idx], marker = 'x', color = 'slategrey', label = r'$s_{min}$ (no constraints)')
        
        ## Accountinf for the constrained s
        mask = Y <=smin
        lcpy = np.where(mask, lcpy, np.inf)
        idx = np.argmin(lcpy, axis = 0)
        ax.scatter(rhat,s[idx], marker = '.', color = 'magenta', label = r'$s_{min}$ (constrained)')
        # ax.scatter(rhat,s[idx], marker = 'o', edgecolor = 'magenta', facecolor = 'none', label = r'$s_{min}$ (constrained)')
        # lines = ax.contour(X, Y, ls[i], levels=10, colors="white", linewidths=0.7)
        # ax.clabel(lines, inline=True, fontsize=8)  # label contour lines
        # Add individual colorbar
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Value")

    # Tight layout for spacinghandles, labels = axes[0].get_legend_handles_labels()
    handles, labels = axes[0].get_legend_handles_labels()
    # One global legend outside the figure
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.01), ncols = len(labels))
    plt.tight_layout()
    plt.savefig(f'l_svr_alpha_{alpha:.2f}_p_{p}.png', dpi = 300, bbox_inches = 'tight')
    plt.show()
           #def appl_func(s,rhat):
           #    X,Y = np.meshgrid(s, rhat, indexing='ij')
           #    print('Grid \n', X)
           #    print(Y)
           #    return len_gen_ab(p, qh, alpha, Y,X )
           #ls.append(appl_func(s, rhat))
           #print(ls[-1])
