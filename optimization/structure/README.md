Structure update submodule
====================================

Updates the value of z, 
in order to minimize the sum of Q(z) + g(z), 
where g(z) is the structure objective function.


Description
-----------

The structure submodule allows for four update schemes
with which to update p, where g(z) = I_0(z - m(p)) + w(p):

*   'continuous' which bounds the values of p to a continuous allowable range,
*   'continuous-linear' which additionally signals that both m(p) and w(p)
    are linear functions,
*   'discrete' which sets the range of p 
    to a finite set of descrete values, and
*   'discrete-diagonal' which additionally signals that m(p), w(p), and Q(z)
    are linear and diagonal functions.
