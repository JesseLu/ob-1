Description
-----------

Objective-first design for nanophotonic structures.


Outline
-------

Ob-1 is divided into four modules:

*   physics module:
    Enables the user to describe the desired functionality of the structure 
    to be designed.
    Returns a design specification in the language of electromagnetics.

*   simulation module:
    Provides the functionality needed to solve the electromagnetic wave equation
    for a variety of scenarios.
    Critically, allows for batch solves of multiple three-dimensional simulations
    through the Maxwell package (maxwell.lightlabs.co).

*   translation module:
    Combines the design specification from the physics module
    with the capabilities of the simulation module to produce
    an optimization problem that is entirely in the language of linear algebra.

*   optimization module:
    Given an optimization problem, attempts to produce an optimal design for it.


Example usage
-------------



Do-now
------

*   Implement the physics and simulation modules.


Do-later
--------


