$l$ \~ Genome Length

$2^k$ fixed points

$k/n = 1 - H(r)$, where $r$ is the complexity parameter. It is the relative distance to end up in a basin.\


Initial condition yields us, $log(2^k) = k$ bits of information. If $l$ is the genome length, then the total information is $k+l$ (verify if this is true since it shall only be true if $I(I n i t; G e n o m e) = 0$). Then the _required information_ is 
$
  k/n + l/n = 1 - H(r) + l/n
$

== Distance between two strings is the linear combination of the Relative Hamming Distances
Consider two strings of length $n$ divided into $t$ substrings of length ${alpha_1n, ...alpha_t n}$, and the relative distance between the $i$th substrings is given by ${q_i}$, then the total relative distance between the two strings is:
$
  q_(t o t) = sum alpha_i q_i
$


== Cell + Random substring

Let us assume the state of the beaker described by a $n$ length string with initially the left fraction $alpha$ describing the cell and the right $(1-alpha)$ being a random substring. The left substring is the _foothold_ and after a time the string should evolve to 3 substrings with the two left alpha fractions being close to each other in distance and the rightmost substring being random.