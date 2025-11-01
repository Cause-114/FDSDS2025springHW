n = 20;
A = eye(n);
for i = 1:n - 1
    A(i, i + 1) = 1;
end
A(n, 1) = 1;A = A / 2;
x = randi([-n, n], n - 1, 1);
y = randi([-n, n], n - 1, 1);
x = [x; -sum(x)];
y = [y; -sum(y)];
x = x / norm(x);
y = y / norm(y);

phase1(A, x, y);
phase2(A, x, y);
phase3(A, 0, pi*37/180);

function phase1(A, x, y)
    figure;
    subplot(2, 2, 1);
    plot([x;x(1)], [y;y(1)],'b-*');  
    title('Original');
    for i = 1:100
        x = A * x;
        y = A * y;
        if(i==5)
            subplot(2, 2, 2);
            plot([x;x(1)], [y;y(1)], 'b-*');
            title('After 5 iterations');
        elseif(i==20)
            subplot(2, 2, 3);
            plot([x;x(1)], [y;y(1)],'b-*');  
            title('After 20 iterations');            
        end
    end
    subplot(2, 2, 4);
    plot([x;x(1)], [y;y(1)],'b-*');  
    title('After 100 iterations');
end

function phase2(A, x, y)
    figure;
    subplot(2, 2, 1);
    plot([x;x(1)], [y;y(1)],'b-*');
    title('Original');
    for i = 1:200
        x = A * x;
        y = A * y;
        x = x / norm(x);
        y = y / norm(y);
        if(i==5)
            subplot(2, 2, 2);
            plot([x;x(1)], [y;y(1)],'b-*');
            title('After 5 iterations');
        elseif(i==20)
            subplot(2, 2, 3);
            plot([x;x(1)], [y;y(1)],'b-*');  
            title('After 20 iterations');
        end
    end
    subplot(2, 2, 4);
    plot([x;x(1)], [y;y(1)],'b-*');
    title('After 200 iterations');
end

function phase3(A, theta1, theta2)
    figure;
    n=size(A,1);
    c=cos(0:2*pi/n:2*(n-1)*pi/n);
    s=sin(0:2*pi/n:2*(n-1)*pi/n);
    x=c.*cos(theta1)+s.*sin(theta1);
    y=c.*cos(theta2)+s.*sin(theta2);
    x=x(:);y=y(:);
    for i = 1:100
        x = A * x;
        y = A * y;
        x = x / norm(x);
        y = y / norm(y);
        % if i < 50
        %     continue;
        % end
        if mod(i, 2)
            plot(x, y, 'b*');
            hold on;
        else
            plot(x, y, 'r*');
            hold on;
        end
    end
end
