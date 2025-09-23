= Finite Length Stirling's Approximation

#set math.equation(numbering:"(1)")
*Stirling's Approximation*

$
n! tilde n^n e^(-n) sqrt(2 pi n)
$<stir>

*Log Sum Inequality*
For two sets of nonnegative numbers ${a_1, ... a_t}$ and ${b_1, ... b_t}$. The log sum inequalty states that:
$
sum_(i) a_i log(a_i/b_i)  gt.eq (sum_i a_i) times  log((sum_i a_i)/(sum_i b_i)  )
$<logsum>
the inequality becomes an equality when $b_i = c a_i$ for all $i$.
== Approximating the multinomial coefficient
To approximate the multinomial coefficient $vec( n, q_1 t ..q_t t)$, we use the Stirling's appriximation (@stir):
$
binom( n, q_1 t,...,q_t t) &= n^n e^(-n) sqrt(2 pi n)/ (product_i ( q_i^(q_i n) n^(q_i n) e^(-n q_i) sqrt(2 pi q_i n) ) ) \
&tilde cancel(n^n e^(-n)) sqrt(2 pi n) /  ( cancel(n^n e^(-n)) (2 pi n)^(t/2) product_i (q_i)^(n q_i) (q_i)^(1/2) ) \
&= (2 pi n) ^(1/2 - t/2) 2^(n H(q_i) )times (1/ (product_i q_i) )^(1/2)  
$<plug>


Using the log-sum inequality, we can say that 
$
log(1/ (product_i q_i) )  = sum_(i=0)^t 1 log&(1/q_i) gt.eq (sum_(i=0)^t 1)times log( (sum_(i=0)^t 1)/(sum_(i=0)^t q_i))\
log(1/ (product_i q_i) ) &gt.eq t  log( t/(1)) = log(t^t)\
$

$
1/ (product_i q_i)  &gt.eq t^t
$<logsumres>

going back to @plug, we then have 

$
binom( n, q_1 t,...,q_t t) = (2 pi n) ^(1/2 - t/2) 2^(n H(q_i) ) (1/ (product_i q_i) )^(1/2)  gt.eq (2 pi n) ^(1/2 - t/2) 2^(n H(q_i) ) t^(t/2)
$



The inequality becomes equaltity when $q_i^* = c times 1$ where c is some constant and from the normalisation condition on $sum_i^t q^*_i=1 $ implies that $c = q_i^* =  1/t$.

$
binom( n, q_1 t,...,q_t t) approx (2 pi n) ^(1/2 - t/2) 2^(n H(q_i) ) t^(t/2)
$


