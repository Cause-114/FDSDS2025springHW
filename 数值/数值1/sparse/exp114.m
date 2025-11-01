log_absolute_err=zeros(1,199);
log_relative_err=zeros(1,199);
for n=2:200
    dia = rand(1, n);
    U = rand(n, n);
    A = U \ diag(dia) * U;
    anss = U \ diag(exp(dia)) * U;
    res = my_exp(A);
    log_absolute_err(n-1) = log10(norm(res - anss));
    log_relative_err(n-1) = log10(norm(res - anss) / norm(anss));
end
plot(2:200, log_relative_err,"DisplayName", "Relative Error");
hold on;
plot(2:200, log_absolute_err, "DisplayName", "Absolute Error");
legend show;
xlabel("n");
ylabel("log10(norm)");
title("Relative and Absolute Error of my expm Function");
grid on;

function R = my_exp(A)
    k = ceil(max(log2(norm(A, inf)), 0));
    A = A / 2 ^ k;
    R = eye(size(A, 1));
    A_ = A; x = 1;
    i = 1;
    while norm(A_) > 1e-20 * norm(R)
        R = R + x * A_;
        i = i + 1;
        x = x / i;
        A_ = A_ * A;
    end
    for i = 1:k
        R = R * R;
    end
end
