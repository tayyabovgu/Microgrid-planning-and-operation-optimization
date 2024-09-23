function WT1 = wtarrifs(PLoadProf, Horizon)
    kz = 0.18;
    el = sum(PLoadProf)'; % Sum of load profile
    time = (1:Horizon);
    a = reshape(el, 24, 365); % Reshape load profile for 365 days, 24 hours each
    
    % Separate weekdays and weekends
    W_work = a; 
    W_work(:, 7:7:end) = []; % Exclude Sundays from weekdays
    W_weekend = a(:, 7:7:end); % Only Sundays (weekends)

    % Initialize demand vectors
    Din = zeros(24, 1); 
    Dop = zeros(24, 1);
    DP = zeros(24, 1);

    % Calculate weekday tariffs
    for i = 1:313
        for j = 1:24
            if j == 18 || j == 19 || j == 20 || j == 21
                DP(j) = W_work(j, i);
                Din(j) = 0;
                Dop(j) = 0;
            elseif j == 17 || j == 22
                Din(j) = W_work(j, i);
                DP(j) = 0;
                Dop(j) = 0;
            else
                Dop(j) = W_work(j, i);
                DP(j) = 0;
                Din(j) = 0;
            end
        end
        WTwd(:, i) = kz * (5 * DP + 3 * Din + Dop);
        WT.WT_work(:, i) = 3 * WTwd(:, i);
    end

    % Set weekend tariffs
    WT.WT_weekend = 30 * ones(24, length(W_weekend(1, :)));

    % Reshape weekday tariffs
    A = WT.WT_work(:, 1:312);
    AAA = WT.WT_work(:, 313);
    AA = WT.WT_weekend(:, 1:52);

    kkk = 1;
    kk = 1;
    while kk <= 312
        B = A(:, kk:kk+5);
        kk = kk + 6;
        C(:, kkk) = reshape(B, [], 1);
        kkk = kkk + 1;
    end

    kkkk = 1;
    while kkkk <= 52  
        r(:, kkkk) = [C(:, kkkk); AA(:, kkkk)];
        kkkk = kkkk + 1;
    end

    WT1 = reshape(r, [], 1);
    WT1 = ([WT1; AAA])';
end
