%% test_update_structure
% Test the update_structure function.

function test_update_structure(l)

    A_p = randn(l) + 1i * randn(l);
    b_p = randn(l, 1) + 1i * randn(l, 1);
    c_p = randn(l, 1) + 1i * randn(l, 1);
    g = struct( 'm', @(p) A_p * p - b_p, ...
                'w', @(p) real(c_p' * p), ...
                'p_range', ones(l, 1) * [0 1], ...
                'scheme', 'continuous-linear');

    %% Test the continuous-linear scheme
    P = randn(l) + 1i * randn(l);
    q = randn(l, 1) + 1i * randn(l, 1);

    z = update_structure(P, q, g);

    norm(P * z - q)
