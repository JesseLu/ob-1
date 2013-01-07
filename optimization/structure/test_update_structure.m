%% test_update_structure
% Test the update_structure function.

function test_update_structure(l)

    
    %% Test the continuous-linear scheme
    A_p = randn(l) + 1i * randn(l);
    b_p = randn(l, 1) + 1i * randn(l, 1);
    c_p = randn(l, 1) + 1i * randn(l, 1);
    g = struct( 'm', @(p) A_p * p - b_p, ...
                'w', @(p) real(c_p' * p), ...
                'p_range', ones(l, 1) * [0 1], ...
                'scheme', 'continuous-linear');

    P = randn(l) + 1i * randn(l);
    q = randn(l, 1) + 1i * randn(l, 1);

    p0 = rand(l, 1);
    [z, p] = update_structure(P, q, g, p0);


    %% Test the discrete-diagonal sheme
    A_p = diag(diag(A_p));
    P = diag(diag(P));
    g = struct( 'm', @(p) A_p * p - b_p, ...
                'w', @(p) real(c_p' * p), ...
                'p_range', ones(l, 1) * [0 1], ...
                'scheme', 'discrete');

    p0 = randi(2, l, 1) - 1;
    [z, p1] = update_structure(P, q, g, p0);


    %% Test the discrete-diagonal sheme
    g = struct( 'm', @(p) A_p * p - b_p, ...
                'w', @(p) real(c_p' * p), ...
                'p_range', ones(l, 1) * [0 1], ...
                'scheme', 'discrete-diagonal');

    p0 = randi(2, l, 1) - 1;
    [z, p2] = update_structure(P, q, g, p0);

    if any(p1 ~= p2)
        error('Discrete and discrete-diagonal results differ for linear case.');
    end

