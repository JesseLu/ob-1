Structure parameterization submodule
====================================

Updates the value of z, the structure parameterization
in order to minimize the sum of Q(z) + g(z), 
where g(z) is the structure objective function.


Description
-----------

If the allowable range of z is continous, 
then the structure submodule uses a steepest-descent strategy
to minimize Q(z) + g(z).

However, if the allowable range of z is discrete,
then the structure submodule updates z by trial-and-error
utilizing the assumption that the elements of z are independent of one another.
