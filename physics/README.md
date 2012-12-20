Physics module
==============

Produces a specification of the design/optimization problem in the language
of electromagnetics. 


Description
-----------

The physical specification which the physics module produces consists 
  first of the allowable values of the parameters (z), 
  which are used to describe the structure to be designed.
The range of z can either be one continuous and bounded range,
  or else a finite set of discrete allowable values.

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
*   a function handle that returns the change in epsilon with respect to 
    various values of the structure parameters (z).



