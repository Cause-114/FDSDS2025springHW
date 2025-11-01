f1 = @(x) x;
f2 = @(x)3/2 * x ^ 2 -1/2;
f3 = @(x)5/2 * x ^ 3 -3/2 * x;
df1 = @(x) 1; df2 = @(x) 3 * x;
df3 = @(x)15/2 * x ^ 2 -3/2;
ddf1 = @(x) 0; ddf2 = @(x) 3; ddf3 = @(x) 15 * x;
r_list = root_solve(f1, df1, ddf1, 1);
disp(r_list);
r_list = root_solve(f2, df2, ddf2, 2);
disp(r_list);
r_list = root_solve(f3, df3, ddf3, 3);
disp(r_list);

function r_list = root_solve(f, df, ddf, n)
    % make sure your input has no repeated roots
    % and all the roots loated on the interval [-1,1]
    r_list = zeros(n, 1); id = 0;
    r0 = -1;
    while (true)
        while (true)
            if (df(r0) == 0)
                r0 = r0 + 1e-2;
            end
            r0 = r0 - f(r0) / df(r0);
            if (abs(f(r0))<1e-6)
                break;
            end
        end
        r_list(n - id) = -r0; id = id + 1;
        r_list(id) = r0;
        if(id>=n/2)
            break;
        end
        while (true)
            if(ddf(r0)==0)
                r0=r0+1e-6;
            end
            r0 = r0 - df(r0) / ddf(r0);
            if (abs(df(r0)) < 1e-6)
                break;
            end
        end
        r0 = r0 + 1e-2;
    end
    if (bitand(n, 1))
        r_list(ceil(n / 2)) = 0;
    end
end
xx=linspace(-1,1,100);
plot(xx,arrayfun(f1,xx));hold on;
plot(xx,arrayfun(f2,xx));hold on;
plot(xx,arrayfun(f3,xx));hold on;
grid on;