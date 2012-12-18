%% test_line_search_convex
% Test the |line_search_convex| function.
    n = 10;
    A = 10 * randn(10) + 1i * randn(10);
    b = randn(10, 1) + 1i * randn(10, 1);
    x0 = randn(10, 1) + 1i * randn(10, 1);

    f = @(x) 1/2 * norm(A*x - b)^2;
    grad_f = @(x) A' * (A*x - b);

    line_search_convex(f, grad_f, -grad_f(x0)*1e-0, x0, 1e-6)
