%% update_structure
% Updates z by solving argmin Q(z) + g(z).

%% Description
% Q(z) is described by matrix P and vector q, via Q(z) = 1/2 || Pz - q ||^2.
%
% Typically, the optimal value for z is found by solving the problem in terms
% of its parameterization variable p.
% For this reason, the problem is typically converted into "p-space".

function [z] = update_structure(P, q, g)
    


    switch(g.scheme)
        case 'continuous'

        %% Continuous-linear case
        %
        % We transform $Q(z) + g(z)$ into p-space via
        % $$ 1/2\|P(A_m p - b_m)-q\|^2 + \mbox{real}(c_w^\dagger p) = 
        %    1/2\|P A_m p - (q + b_m + (P A_m)^{-\dagger} c_w)\|^2, $$
        % where we have used the linearizations of m(p) and w(p).
        case 'continuous-linear'
            path(path, genpath(strrep(mfilename('fullpath'), ...
                                            'update_structure', 'cvx')));
            p0 = g.p_range(:,1);

            % Linearization of m(p).
            A_m = get_gradient(g.m, p0);
            b_m = g.m(p0);

            % Linearization of w(p).
            c_w = get_gradient(g.w, p0)';

            % Transform into p-space.
            A = P * A_m;
            b = q + b_m + (A' \ c_w);

            cvx_quiet(true)
            cvx_begin
                variable p(length(p0))
                minimize norm(A*p - b)
                subject to
                    p <= g.p_range(:,2)
                    p >= g.p_range(:,1)
            cvx_end

            % Test
            real(A'*(A * p - b))
            p
            z = g.m(p);
        case 'discrete'
        case 'discrete-diagonal'
        otherwise
            error('Invalid scheme for structure design objective.');
    end

       
    
