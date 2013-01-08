%% run_optimization
% Obtain an optimized structure.

%% Description

function [z, p] = run_optimization(opt_prob, g, p0)
    p = p0;
    z = g.m(p);
