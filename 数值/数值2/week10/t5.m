Richardson(10);
function Richardson(n)
    T = zeros(1, n + 1); m = 3;
    T(1) = m*sin(pi/m);
    for i = 1:n
        m = m * 2; d = 4;
        tmp = m * sin(pi/m);
        for j = 2:i + 1
            tmp2 = (tmp * d - T(j - 1)) / (d - 1);
            T(j - 1) = tmp; tmp = tmp2; d = d * 4;
        end
        if(abs(pi-tmp)<5e-8)
            break;
        end
        T(i + 1) = tmp;
    end
    fprintf("minimum value of n to get 7 digits: %d\n",m);
end