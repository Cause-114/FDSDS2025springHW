f = @(x) 1 + cos(x);
df = @(x) -sin(x);

function root = newton_solve(initial_guess, f, df, tolerance)
    root_list = zeros(1, 100); index = 1;
    root = initial_guess; root_list(1) = root;
    while index == 1 || abs(root_list(index - 1) - root) > tolerance
        root = root - f(root) / df(root);
        index = index + 1;
        root_list(index) = root;
        if index > 100
            fprintf('Maximum iteration reached.\n'); break;
        end
    end
    plot(1:index - 1, log10(abs(root_list(1:index - 1) - root)));
end

r = newton_solve(3, f, df, 1e-9);
fprintf('Root is %.8f\n', r);
title("Root convergence history")
ylabel('$log_{10}|error|$', 'Interpreter', 'latex');
xlabel('Iteration time');
