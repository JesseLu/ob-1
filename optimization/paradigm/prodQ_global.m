%% prodQ_global
% Produce Q(z) using a global (ADMM) optimization paradigm.

%% Description
%

function [P, q, status] = prodQ_global(z, opt_prob, state, varargin)

%% Input parameters

%% Output parameters

    %% Parse inputs
    t = 1;
    phi = 0;
    xv = varargin{1};

    % Short-hand variables from the values of the opt_prob parameter.
    fobj = [opt_prob.field_obj];
    pres = [opt_prob.phys_res];
    invA = {opt_prob.solve_A};
    invAd = {opt_prob.solve_A_dagger};

    N = length(fobj); % Number of fields/modes.
    n = length(z);

    %% Form the relaxed form of the field design objective
    f = @(alpha, beta, C, t, phi, x) - (1/t) * ...
                (ln(beta.^2 - norm(C'*x)^2) + ...
                ln(real(exp(i*phi)*C'*x) - alpha));
    df = @(alpha, beta, C, t, phi, x) (1/t) * ...
                diag((C'*x)./(beta.^2 - norm(C'*x)^2) + ...
                    exp(i*phi)./(real(exp(i*phi)*C'*x) - alpha)) * C';
    f_indiv = @(alpha, beta, C, t, phi, x) - (1/t) * ...
                (beta.^2 - abs(C'*x).^2);

    f_indiv = @(alpha, beta, C, t, phi, x) - (1/t) * ...
                (beta.^2 - norm(C'*x).^2);
    df = @(alpha, beta, C, t, phi, x) (1/t) * ...
                diag((conj(C'*x)./abs(conj(C'*x)) ))* C';
    
    f = @(alpha, beta, C, t, phi, x) (-1/(2*t)) * ln(beta.^2 - abs(C'*x).^2);
    df = @(alpha, beta, C, t, phi, x) diag(conj(C'*x)./(beta.^2 - abs(C'*x).^2)) * C';

    f = @(alpha, beta, C, t, phi, x) (1/2) * sum(beta.^2 - abs(C'*x).^2);
    grad_f = @(alpha, beta, C, t, phi, x) C * (C'*x);
    Hess_f = @(alpha, beta, C, t, phi, x) (C) * (C');

    f = @(alpha, beta, C, t, phi, x) -(1/2) * sum(ln(beta.^2 - abs(C'*x).^2));
    grad_f = @(alpha, beta, C, t, phi, x) sum((diag((C'*x)./(beta.^2 - abs(C'*x).^2)) * C')',2);
    Hess_f = @(alpha, beta, C, t, phi, x) C * diag(1./(beta.^2-abs(C'*x).^2) + ...
                                        (C'*x)./(beta.^2-abs(C'*x).^2).^2) * (C');
%     f = @(alpha, beta, C, t, phi, x) 0.5 * norm(beta - x)^2;
%     df = @(alpha, beta, C, t, phi, x) x';
    for k = 1 : N
        alpha = fobj(k).alpha;
        beta = fobj(k).beta;
        C = fobj(k).C;

        fi = @(x) f(alpha, beta, C, t, phi, x);
        grad_fi = @(x) grad_f(alpha, beta, C, t, phi, x);
        Hess_fi = @(x) Hess_f(alpha, beta, C, t, phi, x);
        x = xv{k};
        derivative_tester(fi, grad_fi(x)', fi(x), x, 1e-6, @real);
        derivative_tester(grad_fi, Hess_fi(x)', grad_fi(x), x, 1e-6);
    end

%% Other private functions
function [z] = ln(x)
% Custom log function which returns nan for non-zero imaginary part and 
% real part <= 0.

    if ~all(isreal(x)) || any(real(x) <= 0)
        z = nan;
    else
        z = log(x);
    end

