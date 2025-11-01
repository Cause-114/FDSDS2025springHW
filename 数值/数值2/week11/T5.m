% This function is good for this given "good-quality" DTMF tones,
% i.e. there are alomst no noise in the signal.
% But it may be sensitive to the noise, which needs to be improved.

% Ideas: divide the signal into groups by larget intervals, and find begining
% and ending time of each sound in each group. Then use fft to find the frequency
% spectrum of each sound, and use the frequency to find the corresponding digit.
[y, fs] = audioread('DTMF_dialing.ogg');
solvewave(y, fs);

function solvewave(y, fs)
% This function is used to figure out the actual digits of the DTMF tones.
    [lid, rid] = get_stage(y);
    [r, c] = size(lid);
    for i = 1:r
% Reach the non-used buffer space, the same below.
        if (lid(i, 1) == 0)
            break;
        end
        fprintf("\nThe group %d signals: ", i);
        for j = 1:c
            if (lid(i, j) == 0)
                break;
            end
% Here, I want to use my fft_cyy, so I should make sure the number of samples is 2^N,
% if you'd like to use MATLAB-builtin fft, just type "freq = fft(y(lid(i,j):rid(i,j)));"
            center = round((lid(i, j) + rid(i, j)) / 2);
            freq = fft_cyy(y(center - 256:center + 255));
            freq = freq(2:257);
            get_digit(freq, fs);
            fprintf(get_digit(freq, fs) + " ");
        end
    end
end

function [lid, rid] = get_stage(y)
% This fucntion is used to figure out the begin and end time of each sound.
% Typically if we conider all sound between two "large" intervals as one group, then
% each row of lid and rid represent a group; and there are "small" intervals between
% sounds in a group, then each element of lid and rid represent a sound.
    res = 0; stage = false; tmp = 0;
% You can set this buffer size larger according to your need.
    lid = zeros(8, 10);
    rid = zeros(8, 10);
    id1 = 0; id2 = 0;
    for i = 51:length(y)
        res = (res * 50 + abs(y(i)) - abs(y(i - 50))) / 50;
% If the average amplitude of the nearby 50 points is greater than 0.2, we consider it
% in the sound period (stage=true). Otherwise it is considered as interval (stage=false).
        if res > 0.2
            if ~stage
% Begin of a new sound.
                lid(id1, id2) = i - 25;
                stage = true;
            end
        else
            if stage
% End of a sound.
                rid(id1, id2) = i - 25;
                tmp = rid(id1, id2);
                stage = false;
                id2 = id2 + 1;
            elseif (i - tmp == 1000)
                % Then this interval is considered as "large", so we need to start a new group.
                id1 = id1 + 1; id2 = 1;
            end
        end
    end
end

function x = fft_cyy(x)
    % Only for N=2^m length, i.e. N = 1,2,4...
    N = length(x);
    if (bitand(N, N - 1))
        error('N must be a power of 2');
    end
    x = x(:).';
    omega = exp(-2i * pi / N * (0:N / 2 - 1));
    step = bitshift(N, -1);
    while (step >= 1)
        for j = 1:step
            x(j + step:2 * step:N) = x(j + step:2 * step:N) .* omega(1:step:N / 2);
            tmp = x(j:2 * step:N) + x(j + step:2 * step:N);
            x(j + N / 2:step:N) = x(j:2 * step:N) - x(j + step:2 * step:N);
            x(j:step:j + N / 2 - 1) = tmp;
        end
        step = bitshift(step, -1);
    end
end

function res = get_digit(freq, fs)
% This function is used to figure out the actual digit in a given frequency spectrum.
    answer = ["1", "2", "3", "A"; "4", "5", "6", "B"; "7", "8", "9", "C"; "*", "0", "#", "D"];
    lower_freq = [697, 770, 852, 941];
    upper_freq = [1209, 1336, 1477, 1633];
    freq = abs(freq);
    basefreq = fs / length(freq) / 2;
    lf = 0; uf = 0;
% find the lower and upper frequency of the given frequency spectrum.
    for i = 1:length(freq)
        if (freq(i) > 30 && freq(i) > freq(i + 1) && freq(i) > freq(i - 1))
            if (lf == 0)
                lf = i * basefreq;
            elseif i * basefreq >= upper_freq(1) - 100
                uf = i * basefreq;
                break;
            end
        end
    end
% find the closest frequency in the lower_freq and upper_freq.
    r1 = abs(lower_freq(1) - lf);
    r2 = abs(upper_freq(1) - uf);
    id1 = 1; id2 = 1;
    for i = 2:4
        if (abs(lower_freq(i) - lf) < r1)
            id1 = i; r1 = abs(lower_freq(i) - lf);
        end
        if (abs(upper_freq(i) - uf) < r2)
            id2 = i; r2 = abs(upper_freq(i) - uf);
        end
    end
    res = answer(id1, id2);
end
