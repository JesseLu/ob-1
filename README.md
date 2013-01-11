Description
-----------

Objective-first design for nanophotonic structures.


Do-now
------

*   Figure out why there is no "wiggle room" in prodQ_global.
*   Incorporate the progress tracking in restarts.
*   When alpha = 0, remove lower limit for prodQ_global.
*   Generalize track_progress to be able to select the plane and component
    on which to visualize E and epsilon.


Do-later
--------



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


Application ideas
-----------------

*   Solar cell: normal incidence plane wave on to a lossy membrane (Si) with mirror.
    Objective would be no reflection.
*   CCD camera cell: normal incidence plane wave at red, green, and blue.
    The underlying structure would be three point-like detectors.
    Objective would be appropriate field amplitude at each detector.
*   Lens
*   Nanophotonic interconnect
*   LED: active qw layer that must radiate surface normal. 
    Would need to try lots of positions.
*   Photonic crystal laser: central point source that must radiate into the normal
    direction, with Gaussian profile, and high build-up.
    Alternatively, can input Gaussian, and output amplified Gaussian,
    with only a small amplifying volume at the center.

Example usage
-------------


