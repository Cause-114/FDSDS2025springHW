n = 15;
A = randi([1, 100], n, n);
A = A * A';
disp(A);
B = A;
for i = 1:n
    A(i, i) = sqrt(A(i, i));
    A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
    for j = i + 1:n
        A(j:n, j) = A(j:n, j) - A(j:n, i) * A(j, i);
    end
end

L = zeros(n, n);
for i = 1:n
    L(i:n, i) = A(i:n, i);
end
disp(L * L');
