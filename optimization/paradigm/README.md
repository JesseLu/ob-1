Optimization paradigm submodule
===============================

Based on the optimization state, the paradigm submodule produces 
  a quadratic problem, Q(z), 
  which can be passed to the structure submodule
  in order to update z.
The paradigm submodule can user either 
  a local (adjoint) or a global (objective-first) optimization paradigm.


Description
-----------

The paradigm submodule produces Q(z) given

*   the current value of z,
*   the quantities for the field design objectives: alpha, beta, and C,
*   the quantities for the physics residual: A, b, B, and d,
*   function handles to solve for A as well as its conjugate transpose, 
*   the previous state-variable of the paradigm submodule and
*   various user-defined settings.

The specific paradigm can be chosen by using 
  either the prodQ_local() or prodQ_global() functions.

The quadratic problem, Q(z), is returned in the form of matrix P and vector q,
 where Q(z) = 1/2 * || Pz - q ||^2.  
