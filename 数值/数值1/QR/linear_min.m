function delta = f(y1, y2)
    delta = 0;
    n = size(y1, 1);
    for i = 1:n
        delta = delta + (y1(i) - y2(i)) ^ 2;
    end
end

x_list = [2; 3; 4; 5; 6; 7]; x_copy = x_list;
y_list = log(x_list); one = ones(6, 1);

R = zeros(2, 2);
R(1, 1) = norm(x_list);
x_list = x_list / R(1, 1);
R(1, 2) = dot(x_list, one);
one = one - R(1, 2) * x_list;
R(2, 2) = norm(one);
one = one / R(2, 2);
a = R \ [dot(x_list, y_list); dot(one, y_list)];

one = ones(6, 1);
y_heta = [x_copy, one] * a;

plot(x_copy, y_list, 'o', 'DisplayName', 'ln(x)');
hold on
plot(x_copy, y_heta, 'DisplayName', 'linear fitting');
legend show
fprintf("Min delta located at a=%f , b=%f\n", a(1), a(2));

delta_y2 = f(y_heta, y_list);
min_delta_y2 = 65535;
for i = 0:0.0001:0.5
    for j = 0:0.0001:0.5
        y_heta1 = x_copy * i + one * j;
        min_delta_y2 = min(min_delta_y2, f(y_heta1, y_list));
    end
end
fprintf("The min delta is %f,and %f if find precisely\n", delta_y2, min_delta_y2);
% What should be mentioned is that in order to get more precise ans,
% you need to make the step of i,j smaller! 