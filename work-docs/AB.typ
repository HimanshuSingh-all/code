#let numbered_eq(content) = math.equation(
    block: true,
    numbering: "(1)",
    content,
)

#let geq = $gt.eq$
#let leq = $lt.eq$
#let hr = $hat(r)$
#let cpl = $ plus.circle$

#let ti(var) = $ tilde(var) $
#set par(justify:true)
#show heading.where(level: 1): set align(center)
#set page(header: [
  Himanshu Singh, NCBS
  #h(1fr)
  #datetime.today().display()
], 
numbering: "(1)")
=  Self Organised and Function + Pure function

We assume the final string representing the stable FP as a combination of substring of the self organised and functional substring (denoted by A) plus the purely functional substring (denoted by B).


#set heading(numbering: "1.1")
== Some diagrams and introduction to the problem

#figure(
  grid(
      columns: (auto, auto),
      rows:(auto, auto),
      gutter: 1em,
      [#image("sofunctional.jpeg", width: 80%)],
      [#image("functional.jpeg", width: 80%)],
      [#text("(a)")],
      [#text("(b)")],
      [#image("complexityparam.jpeg", width: 80%)],
      [#image("abstringrep.jpeg", width: 100%, height:15%)],
      [#text("(c)")],
      [#text("(d)")]
), 
caption: []
)<strc>

*These are the following fixed parameters: *
- $p$, the developmental noise, specifically the relative distance due to the developmental noise. 
- $q$, the complexity parameter related to the number of fixed points.
- $hr$, which the relative distance between some ideal _copy_ to the noisefree copy.
Where both $0 leq p,q leq 1/2$. We see a case (@rboundsection) where $hr$ and hence $r= hr cpl p$, be constrained in the following region $0 leq hr <1/4$, in order for the inqualities obtained by the positivity of the probabilities to be consistent. 

The information (in bits) given by the initial condition is given as:
#numbered_eq(
$
k = n(alpha)(1-H(q))
$
)<k>

Now going back to @strc(c) and @strc(d) 
, we may easily write the following  inquality on the initial string:
#numbered_eq(
$
(1-alpha)1/2 + alpha(p cpl s cpl p) leq q
$
)
#numbered_eq(
$
p cpl (s alpha + 1/2 (1-alpha)) cpl p leq q\
p cpl overline(s) cpl p leq q
$
)<init>
where,
$
overline(s) = (s alpha + 1/2 (1 - alpha))
$
Going back to @init, we may now write:
#numbered_eq(
$
&overline(s) leq 1/2 (( q - p cpl p)/(1/2 - p cpl p) )
$
)

#numbered_eq(
$&s leq 1/(2alpha)[ (( q - p cpl p)/(1/2 - p cpl p) ) - 1 ] + 1/2 = s_(m a x)$
)<smaxfirst>

for ease we can write 
$ 1 - (( q - p cpl p)/(1/2 - p cpl p) ) =  D(q,p) $

since $ 0 leq q leq 1/2$, $D(q,p) geq 0$. From @smaxfirst we have the following:
#let sm = $s_(m a x)$
#numbered_eq(
$ s leq-D(q,p)/(2 alpha) + 1/2  = sm$
)<smax>
Let us now take a moment to  rewrite all the constraints and inqualities that we have so far:


 $ 0 leq p leq 1/2, $
 $ 0 leq q leq 1/2, $
 $ s leq -D(q,p)/(2 alpha) + 1/2  = sm, $
 $ 0 leq alpha leq 1, $
 $ beta = 1 - 2 alpha, $
 $ r = hr cpl p $

== Contribution of the A term
/*
#figure(
image("aprob.jpeg"),
  caption:[Computing the probability of observing a two strings drawn randomly with a distance $s$ while having a mutual distance of $hr$ from some other fixed string ]
)
*/
#let ts = $tilde(s)$
Without going into details we just write the probability of the strings is:
#numbered_eq(
$
P = 2^(-l_s D_(K L)({ts_i||p_i}))
$
)
in our case $l_s = alpha n$. We have the observed probabilities as:
#numbered_eq(
$
ts_1 = 1 - 1/2(2 hr +s)
$
)
#numbered_eq(
$
ts_2 = ts_3 = s/2
$
)
#numbered_eq(
$
ts_4 =  (hr -s/2)
$ 
)

And assuming that the generation of the strings is random ($p_i=1/4$), we have:
#numbered_eq(
$
D_(K L)({ts_i||p_i}) = 2 - H({ts_i})
$
)
and we can use that to get,
#numbered_eq(
$
P = 2^(-alpha n(2 - H(ts_i)))
$
)

and then we can compute $log(1/P)$ as,
#numbered_eq(
$ 
log(1/P) = alpha n(2 - H(ts_i))
$
)<abitcontri>

== Contribution of the B part
The contribution of the B part is straightforward to write:

#numbered_eq(
$
log(1/P_B) = beta n (1- H(r))
$
)<bcontrib>

== A and B string together 

From @k, @abitcontri and @bcontrib, we write:

#numbered_eq(
$
k+l = alpha n(alpha) (2 - H(ts_i)) + beta n(alpha)(1-H(r))
$)
Plugging in the value of $k$ from @k and subtracting it from the right-hand side:

$
l = alpha n(alpha) (2 - H({ts_i})) + beta n(alpha)(1-H(r)) - n(alpha) (1 - H(q))
$
#numbered_eq(
$
l = n(alpha)[ alpha  (2 - H({ts_i})) + beta (1-H(r)) - (1 - H(q))]
$
)
remembering that $beta = 1-2 alpha$,
#numbered_eq(
$
l = n(alpha)[ alpha  (2 - H({ts_i})) + (1-2alpha) (1-H(r)) - (1 - H(q))]
$
)<lequation>

where,
#numbered_eq(
  $
  n(alpha) = n_0/(1- alpha)
  $
)<strlen>
Let us rename the following:
#numbered_eq(
  $
  A({ts_i}) = (2- H({ts_i}))
  $
)<A>

#numbered_eq(
  $
  B(r) = (1- H(r))
  $
)<B>

#numbered_eq(
  $
  C(q) = (1- H(q))
  $
)<C>
Now we may rewrite @lequation as
#numbered_eq(

  $
  l = n(alpha)[ alpha A({ts_i}) + (1-2alpha)B(r) - C(q)]
  $
)<leqnabc>

Now again since $n(alpha)geq 0$, the positivity constraint on $l$ implies that
#numbered_eq(
  $
  [ alpha A({ts_i}) + (1-2alpha)B(r) - C(q)] geq 0
  $
)<lpostivty>

for now I have not looked much into this constraint but it should be kept in mind moving forward, since it is a much general constraint that is true irrespective of the domains we shall talk about. (In Mathematica we simply assert $l geq 0 $). If this is positivity constraint is not kept in mind the optimizer fails, citing imaginary values of the function.

#text(blue)[*TODO*: Implement a numerical solver to find minimum l given the bounds on the $alpha$ and $hr$. Keeping $q$ and $p$ fixed. Also need to make sure all these bounds that we pass would be consistent.]

On taking the derivative of @leqnabc with respect to $alpha$, we have
#numbered_eq(
$
(d l)/(d alpha) = (n(alpha)/(1-alpha)) [ A({ts_i}) - B(r) - C(q) + alpha(1-alpha) sum_i (d A)/(d ts_i) (d ts_i)/(d alpha) ] 
$
)<dldalp>

_I might omit the function arguments of $A$, $B$ and $C$ sometimes, it is implied that they are functions of ${ts_i}$
, $r$ and $q$ respectively._
== The Two Regimes $hr cpl hr< sm$ and $hr cpl hr geq sm$

There are two regimes,
- if $hr cpl hr< sm$, then $s = hr cpl hr$ is a _natural_ distance.
- if  $hr cpl hr geq sm$, we use $s = sm$, and in that case $s=s(alpha)$ and $ts_i = ts_i (alpha)$

#let der(u, D) = $(d #u)/(d #D)$
=== Regime 1: $hr cpl hr < sm$ 

In this regime, $der(l, alpha) = 0$, then we have:

#numbered_eq(
$
der(l, alpha) = (n(alpha)/(1-alpha)) [ A - B - C]
$
)

#numbered_eq(
$
A - B - C = H(r) + H(q) - H({ts_i})
$
)

#numbered_eq(
$
 der(l, alpha) = n(alpha)/(1-alpha)[H(r) + H(q) - H({ts_i})]
$
)

=== Regime 2: $hr cpl hr geq sm$ 

In the unexpanded form, (leaving all the $ts_i$'s as is) the derivative of $l$ in this regime is given by

#numbered_eq(
$
der(l, alpha) = n(alpha)/(1-alpha)[A(alpha) - B - C + alpha(1-alpha)D(q,p)/(4alpha^2)log_2((ts_1(alpha) ts_4(alpha))/(s_2^2(alpha)))]
$
)
==== Positivity constraint on probabilities will give a bound on $hr$ <rboundsection>
Now in this regime we put $ s= sm$ and then the probabilies ${ts_i}$ are a function of $alpha$. There is a positivity constraint on the probability.

$
ts_1 = 1 - 1/2(2hr + s(alpha)) geq 0\
ts_2 = ts_3 = s(alpha) geq 0 \
ts_4 = hr - s(alpha)/2 geq 0
$

from these constraints and plugging $ s= sm$ from @smax, we get the following equations
$
sm leq 2(1-hr)
$


$
sm leq 2hr
$
and since $hr leq 1/2$, these will be redundant and we can just write:

#numbered_eq(
$
sm leq 2hr
$
)<smubound>

and then we have additional constraint from positivity of $ts_2 = ts_3= s/2$

#numbered_eq(
$
sm geq 0 
$
)

After plugging in the value of $sm$, we get the two inequalities:

#numbered_eq(
$
D(q,p) leq alpha
$
)<lbalpha>
and an upper bound from @smubound


#numbered_eq(
$
  alpha leq D(q,p)/(1-4hr)
$
)<ubalpha>

therefore, for the bounds given by @ubalpha and @lbalpha to be consistent, recalling that $D(q,p) geq 0$, we require:
#numbered_eq([#rect(
$
hr < 1/4
$
)])

#text(blue)[*so we cannot set $hr$ to be any value in 0 to $1/2$ like we do with $p$ and $q$.*]

