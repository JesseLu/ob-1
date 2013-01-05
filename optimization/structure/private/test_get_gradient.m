%% test_get_gradient
% Just a simple test of the get_gradient function.

    n = 100;
    A = randn(n) + 1i * randn(n);
    b = randn(n, 1) + 1i * randn(n, 1);
    x0 = randn(n, 1) + 1i * randn(n, 1);
    f = @(x) A*x - b;
    
    grad = get_gradient(f, x0);

    err = norm(A(:) - grad(:)) / norm(A(:))

