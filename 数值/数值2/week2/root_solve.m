clc; clear;
f = @(x) x ^ 64 - 0.1;
r1 = bise_regula_falsi(f, 0, 1, 1e-10, 1);
fprintf("Root 1:%.10f\n", r1);
r2 = bise_regula_falsi(f, 0, 1, 1e-10, 2);
fprintf("Root 2:%.10f\n", r2);
legend show;
title('Root convergence history using bisection and false position method');
ylabel('$log_{10}|error|$', 'Interpreter', 'latex');
xlabel('Iteration time');

function m = bise_regula_falsi(f, l, r, tol, choice)
    m_list = zeros(1, 100); index = 1;
    fl = f(l); fr = f(r);
    while (r - l > tol)
        if choice == 1
            m = (l + r) / 2;
        else
            m = (fl * l - fr * r) / (fl - fr);
        end
        fm = f(m); m_list(index) = m; index = index + 1;
        if (abs(fm) < 1e-12)
            fprintf("break for tiny f(r)!\n"); break;
        elseif (fm * fl < 0)
            r = m; fr = fm;
        else
            l = m; fl = fm;
        end
        if (index > 100)
            fprintf("Maximum iteration reached\n");
            break;
        end
    end
    if (choice == 1)
        plot(1:index - 2, log10(abs(m_list(1:index - 2) - m)), 'DisplayName', 'Bisection method');
    else
        plot(1:index - 2, log10(abs(m_list(1:index - 2) - m)), 'DisplayName', 'False position method');
    end
    hold on;
end
