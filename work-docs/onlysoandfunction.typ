
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
= S.O. + functional (only component A)

Let us go back to the information given by the initial condition. Given that the _complexity_ is given by $~1/q$, the equation is:

$
k = n(1-H(q))
$

and we may rewrite $n$ in the form of $alpha$ and $n_0$ as:

$
n = n_0/(1-alpha)
$

Now going back to the image of the Self Organised plus the functional Component.

#figure(
image("so-funct.png"),
caption:[The constraint diagram for self organised copying with a functional constraint, note that the $r^*$ in the figure is just $hat(r)$. I didn't have the functionality to add hat on top of a variable in the paint software that was used to create the image.]
)

We define $r$ to be the folllowing (note that $r^*=hat(r)$).
$
r = hat(r) plus.circle p\
hat(r) = (1/2) (r - p)/(1/2 - p)
$

if the strings are generated iid with bernoulli trials, of length $l_s$ each then the probability of finding a string with some empirical probaility distribution of the element as ${q}$ when the probability of the underlying thing is given by ${p}$ is:
$
cal(P) = 2^(-l_s D_(K L) (tilde(q)||p) )
$

Thus the probability of finding two strings with the relative distance $s$ with each other while both mutually being at a distance $hat(r)$ from some other string (calculation can be done for a string of all zeros).


from the normalisation condition,

$
tilde(q)_1+ tilde(q)_2+ tilde(q)_3+ tilde(q)_4 = 1
$

and the conditions on relative distances will give three more equations,

$
hat(r) = tilde(q)_2 + tilde(q)_3
$

$
hat(r) = tilde(q)_3 + tilde(q)_4
$

$
s = tilde(q)_2 + tilde(q)_3
$
 one can solve for $tilde(q)$'s and then get the following,

$
tilde(q)_1 = 1 - 1/2 ( 2 hat(r) +s),\
tilde(q)_2 = tilde(q)_3 = s/2,\
tilde(q)_4 = hat(r) - s/2
$

for a randomly drawn string with bernulli trial and $p_0 = p_1 = 1/2$ then, we can compute the $D_(K L) (q||p)$ easily,
$
D(q||p) = sum q_i log(q_i/p_i)\
D(q||p) = sum q_i log(4q_i)\
D(q||p) = 2+ sum q_i log(q_i)\
D(q||p) = 2 - H({q})
$

thus the probability of finding these strings is:
$
cal(P) = 2^(-alpha n (2 - H({tilde(q)})) )
$

$
log(1/cal(P)) = alpha n (2- H(tilde(q)))
$

and the number of initial condition and genomic information bits:
$
k+l = alpha n (2 - H({tilde(q_i)}))\
l = alpha n (2 - H({tilde(q_i)})) - k \
l = n( alpha (2 - H({tilde(q_i)})) - (1-H(q) ) )\
l = (n_0/(1-alpha)) [ alpha [2 - H({tilde(q_i)})] - [1-H(q) ] ]
$

Minimising l with respect to $alpha$, we have now:


