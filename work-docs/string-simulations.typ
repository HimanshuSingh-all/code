#import "@preview/cetz:0.4.2"
#show figure.caption: set align(left)
#let targ = `target`
/*#cetz.canvas({
  import cetz.draw: *
  ...
})
*/
// Align the content to the left side of the page.
//

#let rect_width = 6cm
#let rect_height = 17pt
/*
#rect(
  width: rect_width,
  height: rect_height,
  stroke: 1pt + black,
  radius: 0pt,
  inset: 0pt,
  grid(
    columns: (1fr, 1fr),
    column-gutter: 0pt,
    stroke: (y: 1pt + black),
    box(
      width: 100%,
      height: 100%,
      fill: none,
      stroke: 1pt + black,
      pad(x:15pt,y: 5pt, [
        $<dash(A) bar dash(A)'>$
      ])
    ),
    box(
      width: 100%,
      height: 100%,
      fill: luma(220),
      stroke: 1pt + black,
      pad(x:24pt,y: 5pt, [
        .................
      ])
    ),
  )
)
*/
#set heading(numbering: "1.")
#show raw: set text(size: 7pt)

= Simulating Random Strings
Currently everything was done by assuming zero developmental noise $p=0$. So I
might use $r$ and $hat(r)$, $q$ and $hat(q)$ interchangeably.
#figure(
  image("gridsearch.png", width: 100%, height: 7cm,fit: "contain"),
  caption: [ Minimum length of genome required for finite length strings, 
  the white region marks $l_(m i n) lt.eq 0$, and the lowest $r$ for zero length
genome is found for the string simulations. $r_m=0.401$, $hat(q) = 0.35$, $s_m=
0.397$ and $n_0 = 110$ with $alpha= 0.5$.]
)<param>
Using the parameters obtained above as shown in @param, we use the appropriate
$k = ceil(n(1-H(q,n, t=2)))=
19$ to generate $2^(k)$ random strings(fixed points) of length $n_0/(1- alpha)$. Another
string called the _target string_, $X^(targ)$ is also generated, with length $alpha n_0$.
Since $alpha = 0.5$ , we compare the relative hamming distance between the
right and left half of the fixed point strings $s$. The distance $r$ is computed for each fixed
point by taking the max relative distance out of  $X^targ$ and left substring
or the right substring of the fixed point (@randstr). 

#figure(
  image("random_strings_rvs.png", width: 60%, height:8.5cm, fit: "contain"),
  caption: [ $hat(r)$ vs s plot for randomly generated strings of total length
  220 against a randomly generated target. The curve in red is the solution of
  the equation $k/(alpha n ) + H_(t = 4)({1 -hat(r) - s/2, s/2 , s/2, hat(r)-s/2}; n
= alpha n ) - 2 = 0$, which represents the points where the expected number of
strings is 1.]
)<randstr>

#let wo = `worst`
#let ini = $X^(t)_(i n i)$
The scatterplot can be divided into 3 regions by taking two lines $s<s_m$ and
$r<r_m$. For our next subset of simulation we pick the string with maximum $s
= s_wo$ such that $s_wo lt.eq s_m$ and with minimum $hat(r)$. Usually there
are approximately $1$ or $2$ strings. From one of these strings, say $X^(t-1)$ we construct
multiple strings by taking the left half of the left substring of $X^(t-1)$ and right half
of the right substring of $X^(t-1)$ to construct the left half of the initial
condition string, $ini$ and the right half is generated randomly (@figconstr).

#figure(
image("constr-str.png", width: 80%, height: 8cm, fit: "contain"),
  caption: [Generating multiple strings from the selected string with maximum
possible $s=s_wo$, and minimum $r$. The grey part is randomly generated for
  each new constructed string]
)<figconstr>

#let constr = $X_(c o n s t r)$
Now for each of the constructed strings, $constr^i$ ($i in {1, ..,1000}$), we now compute the
minimum relative distance out of all the $2^k$ fixed points strings to it
:$op("min", limits:#true)_j {d^(r e l)(X^j_(F P), constr^i)}$. Additionally we
also compute the distance between each of the constructed strings $constr^i$
and the string used to construct them $X^(t-1)$: $d^(r e l)(constr^i,
X^(t-1))$. We plot the both of these distance for each of the constructed
strings against each other in @constgivenother.

#figure(
image("death_repl.png", fit: "contain", width: 65%, height: 7cm)
)<constgivenother>


The number of constructed strings which are closer to the one of the
random fixed point strings than the string they were constructed from is $456$
out of a 1000 constructed strings, which is denoted by $y>x$ region.
/*
== Helper Classes and Functions <helpers>
```python

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
        Note: This might sometimes throw a DivideByZero warning but that's because the np.where in `BinaryH` function does not perform boolean short circuiting so it will evaluate `a==0` or `a==1` case leading to warnings, but it's still a good practice to know about warnings so I have left them in.
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
```

*/
