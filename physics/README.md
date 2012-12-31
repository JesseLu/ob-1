Physics module
==============

Produces a specification of the design/optimization problem in the language
of electromagnetics. 


Description
-----------

The specification produced by the physics module
consists first of the specification of z, given by
  
*   the structure parameterization function m(p),
*   the parameter weighting function w(p), 
*   the allowable values of p, and
*   the structure update scheme to be used.

In addition, the specification outlines the desired response of the structure
  to various input excitations.
These are called modes, and the specification for each mode consists of

*   the operational frequency,
*   the input excitation,
*   the desired output responses and 
    the corresponding minimum and maximum allowable amplitudes of the response,
*   the s-parameters which determine the spacing of the simulation grid,
*   the values of the permeability on the simulation space,
*   the background values of the permittivity (epsilon), and
*   a selection matrix (S) which describes
    the change in epsilon as a function of z.



