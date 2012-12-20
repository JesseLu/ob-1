Optimization module
===================

Attempts to find a structure parameterization (z) that satisfies the constraints
of an optimization problem of a very specific form.


Description
-----------

The optimization submodule is divided into two separate submodules.

*   The paradigm submodule allows for the choice between a local or global
    optimization strategy. 
    Both these strategies search for optimal z 
    by successively producing quadratic functions Q(z), 
    which are used by the structure submodule to update the value of z.
*   The structure submodule updates z by minimizing 
    the sum of Q(z) and the structure objective function g(z).
    The updated value of z can then be passed to the paradigm submodule
    in order to continue iterating toward an optimal design.
