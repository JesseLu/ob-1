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
                                        2*abs(C'*x).^2./(beta.^2-abs(C'*x).^2).^2) * (C');

    f = @(alpha, beta, C, t, phi, x) ln(real(exp(i*phi)*C'*x) - alpha);
    grad_f = @(alpha, beta, C, t, phi, x) (diag(exp(i*phi)./(real(exp(i*phi)*C'*x) - alpha)) * C')';

    f = @(alpha, beta, C, t, phi, x) sum(ln(real(C'*x) - alpha));
    grad_f = @(alpha, beta, C, t, phi, x) sum(...
                diag(1./(real(C'*x) - alpha)) * C', 1)';
    Hess_f = @(alpha, beta, C, t, phi, x) (C * ...
                diag(1./(real(C'*x) - alpha).^2) * C');

    grad_f = @(alpha, beta, C, t, phi, x) ...
                C * C' * x;
    Hess_f = @(alpha, beta, C, t, phi, x) ...
                C * C';
    % This is basically the working set for alpha (lower limit).
    f = @(alpha, beta, C, t, phi, x) sum(ln(real(C'*x)));
    grad_f = @(alpha, beta, C, t, phi, x) ...
                sum(C * diag(1./real(C'*x)), 2);
    Hess_f = @(alpha, beta, C, t, phi, x, dx) ...
                C * diag(-1./real(C'*x).^2) * real(C' * dx);

    grad_f = @(alpha, beta, C, t, phi, x) ...
                1./(beta^2 - abs(C'*x).^2);
    Hess_f = @(alpha, beta, C, t, phi, x) ...
                1./(beta^2 - abs(C'*x).^2).^2 * (2*C*C'*x);
%     f = @(alpha, beta, C, t, phi, x) 0.5 * norm(beta - x)^2;
%     df = @(alpha, beta, C, t, phi, x) x';
    for k = 1 : N
        alpha = fobj(k).alpha;
        beta = fobj(k).beta;
        C = fobj(k).C;

        C = C(:,1);
        beta = beta(1);
        fi = @(x) f(alpha, beta, C, t, phi, x);
        grad_fi = @(x) grad_f(alpha, beta, C, t, phi, x);
        % Hess_fi = @(x, dx) Hess_f(alpha, beta, C, t, phi, x, dx);
        Hess_fi = @(x) Hess_f(alpha, beta, C, t, phi, x);
        x = xv{k};
        derivative_tester(fi, grad_fi(x)', fi(x), x, 1e-6, @real);
        derivative_tester(grad_fi, Hess_fi(x)', grad_fi(x), x, 1e-6, @real);
%         % For alpha.
%         dt2(grad_fi, @(dx) Hess_fi(x, dx), grad_fi(x), x, 1e-6);
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

