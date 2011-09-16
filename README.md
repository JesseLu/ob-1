Purpose
=======

The purpose of this package is to implement objective-first optimization.


Motivation
==========

When compared to most approaches, the objective-first scheme

1.  requires less work per iteration, which may lead to a working solution
    in a significantly reduced computation time, and
1.  may be less sensitive to starting conditions, especially where non-intuitive
    designs are needed.


What is objective-first optimization?
=====================================

A common optimization problem is to minimize a design objective, while
satisfying physical laws,

    minimize f(x)           [design objective]
    subject to g(x, p) = 0  [underlying physics].

Here, the x variable is typically a field of some kind, while p describes a 
structure. For example, x can be an electric field and p can be the dielectric
structure supporting it.

Most common approach
--------------------

This problem is most often approached by choosing an initial p and then 
repeating the following steps:

1.  solve g(x, p) = 0 for x,
1.  evaluating f(x),
1.  updating p to decrease f(x).

Objective-first approach
------------------------

The objective-first approach reformulates the optimization problem as

    minimize ||g(x, p)||^2  [physics residual]
    subject to f(x) = 0     [design objective].

x is now forced to optimally satisfy the design objective, although x and p
most likely no longer satisfy the underlying physics.

The optimization, then, simply proceeds by updating x and p while ensuring that
f(x) = 0 is always satisfied.


Progress
========

The current goal is to produce a simulation module that can be adapted to 
various physical systems. However, because of time constraints, the simulation
module will only be tested on the sourceless Maxwell equations for
electromagnetics.

Issues currently being worked on
--------------------------------

1.  Change plotting to gnuplot.
1.  Document with Sphinx.


Dependencies
============

This package uses petsc4py and the packages which it subsequently depends on.
Also, the numpy and scipy packages are required.
The Sphinx (python package) is used to generate documentation.
Use of mpi4py and the CUDA environment is anticipated in the near future. 
