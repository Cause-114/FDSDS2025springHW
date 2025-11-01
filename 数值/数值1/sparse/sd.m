A=[20,0;0,1];
b=[0;0];
x0=[1;5];
x=SD(A, b, x0, 1e-6, 100);
disp(x)

function x = SD(A, b, x0, tol, maxiter)
    % This function solves the system Ax=b using the Steepest Descent method
    % with initial guess x0 and tolerance tol. The maximum number of iterations
    % is set to maxiter.
    % Initialize variables
    f=@(x) dot(A*x,x)-2*dot(b,x);
    fprintf("Initial f = %f\n", f(x0));
    x = x0;
    r = b - A * x0;
    x1_=[x(1);zeros(10,1)];
    x2_=[x(2);zeros(10,1)];
    for i = 1:maxiter
        % Compute the step size alpha
        alpha = dot(r, r) / dot(r,A*r);
        % Update the solution
        x = x + alpha * r;
        if(i<=10)
            x1_(i+1)=x(1);
            x2_(i+1)=x(2);
        end
        % Compute the residual
        r = r - A *r* alpha;
        % Check for convergence
        if norm(r) < tol
            fprintf("Converged in %d iterations\n", i);
            break
        end
    end
    figure;
    plot(x1_,x2_);
    title('variations of x');
    xlabel('x_1');
    ylabel('x_2');
    fprintf("Final f = %f\n", f(x));
    fprintf("Residual norm: \n");
    disp(norm(r));
end
