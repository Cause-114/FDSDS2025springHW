% n = 6;
% all = randi([1, 100], 2 * n, 1);
% all = sort(all, "descend");
% dia = all(2:2:end);
% alpha = all(1:2:end);
% z = ones(n, 1);

% for i = 1:n
%     z(i) = alpha(i) - dia(i);

%     for j = 1:n

%         if i ~= j
%             z(i) = z(i) * (alpha(j) - dia(i)) / (dia(j) - dia(i));
%         end

%     end

%     z(i) = sqrt(z(i));
% end

% A = diag(dia) + z * z';

% % here we generate distinct and sorted diagonal elements for the matrix A.
% % And they are quite far from each other, which makes the pucture more beautiful.
% [~, V] = eig(A);
% v = sort(diag(V));
% disp(v);
% disp(alpha);
% % disp(v);
% % disp(dia);

% % In the following code, we plot the result part by part,in order to avoid the Inf value.
% figure;
% x = linspace(dia(1) - 1, dia(1) - 0.0001, 1000);
% y = f(x, z, dia);
% plot(x, y, 'b');
% hold on;

% for i = 2:n
%     l = v(i - 1) - 0.999 * (v(i - 1) - dia(i - 1));
%     r = v(i - 1) + 0.999 * (dia(i) - v(i - 1));
%     x = [linspace(l, v(i - 1), 500), linspace(v(i - 1), r, 500)];
%     y = f(x, z, dia);
%     plot(x, y, 'b');
%     hold on;
% end

% l = v(n) - 0.999 * (v(n) - dia(n));
% r = v(n) + 1;
% x = [linspace(l, v(n), 500), linspace(v(n), r, 500)];
% y = f(x, z, dia);
% plot(x, y, 'b', 'DisplayName', 'function');
% hold on;

% plot(v, zeros(n, 1), '.r', 'DisplayName', 'Real Eigenvalues');
% legend ("fuction", "Real Eigenvalues");

% function y = f(lambda, z, dia)
%     y = zeros(length(lambda), 1);

%     for k = 1:length(lambda)
%         y(k) = 1 + z' * (z ./ (dia - lambda(k)));
%     end

% end

function an = ff(d)
    l = length(d);
    an = 0;
    for i = 1:l
        temp = 1;
        for j = 1:l
            if i ~= j
                temp = temp /(d(j) - d(i));
            end
        end
        % for j = 1:i - 1
        %     temp = temp * (a(j) - d(i)) / (d(j) - d(i));
        % end
        % for j = i + 1:l
        %     temp = temp * (a(j - 1) - d(i)) / (d(j) - d(i));
        % end
        an = an + temp;
    end
end

n = 10;
% all=randi([1,100],1,n);
% all=sort(all);
% d=all(1:2:end);
% a=all(2:2:end);
for j = 1:n
    d=randperm(100);
    d=d(1:n);
    % for i = 1:n
    %     a = randi([1, 100], 1, n - 1);
    %     an = ff(d, a);
    %     disp(an);
    % end
    a=ff(d);
    disp(a);
end
