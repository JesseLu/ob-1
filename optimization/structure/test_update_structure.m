%% test_update_structure
% Test the update_structure function.

function test_update_structure(l_max, test_time)

    start_time = tic;
    
    while toc(start_time) < test_time 
        l = randi(l_max) + 1;
        n = randi(5) * l;
        %% Test the continuous-linear scheme
        A_p = randn(l) + 1i * randn(l);
        b_p = randn(l, 1) + 1i * randn(l, 1);
        c_p = randn(l, 1) + 1i * randn(l, 1);
        g = struct( 'm', @(p) A_p * p - b_p, ...
                    'w', @(p) real(c_p' * p), ...
                    'p_range', ones(l, 1) * [0 1], ...
                    'scheme', 'continuous-linear');

        P = randn(n, l) + 1i * randn(n, l);
        q = randn(n, 1) + 1i * randn(n, 1);

        p0 = rand(l, 1);
        [z, p1] = update_structure(P, q, g, p0);

        %% Test the continuous scheme.
        g.scheme = 'continuous';
        [z, p2] = update_structure(P, q, g, p0);

        f = @(p) 1/2 * norm(P * g.m(p) - q)^2 + g.w(p);
        if (abs(f(p1) - f(p2)) / abs(f(p1))) > 0.1
            error('Continuous and continuous-linear results do not match.');
        end

        %% Test the discrete scheme
        A_p = diag(diag(A_p));
        P = repmat(diag(diag(P)), n/l, 1);
        n
        l
        size(P)
        size(A_p)
        g = struct( 'm', @(p) A_p * p - b_p, ...
                    'w', @(p) real(c_p' * p), ...
                    'p_range', ones(l, 1) * [0 1], ...
                    'scheme', 'discrete');

        p0 = randi(2, l, 1) - 1;
        [z, p1] = update_structure(P, q, g, p0);


        %% Test the discrete-diagonal sheme
        g.scheme = 'discrete-diagonal';
        p0 = randi(2, l, 1) - 1;
        [z, p2] = update_structure(P, q, g, p0);

        f = @(p) 1/2 * norm(P * g.m(p) - q)^2 + g.w(p);
        if any(p1 ~= p2)
            f(p1)
            f(p2)
            error('Discrete and discrete-diagonal results differ for linear case.');
        end
        fprintf('.');
    end
    fprintf('\n');


