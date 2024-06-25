function A = siddon_alg(xi, xf, gri)   
    startpos = xi; % Initial position of ray (both coordinates)
    endpos = xf; % Final position of ray (both coordinates)
    e = 1e-3; % Error margin
    d = 8 / gri; % Pixel size (length of side)
    % Slope of ray
    k = ((endpos(2) - startpos(2)) / (endpos(1) - startpos(1))); 
    X = []; % Array of lengths of paths inside pixels
    C = []; % Vector including indexed pixels 
    prev = 0; % Accrued lengths of paths inside pixels
    x = d - rem(startpos(1), d); % Indicating amount that x needs to be moved when changing pixels
    y = d - rem(startpos(2), d); % Indicating amount that y needs to be moved when changing pixels
        while x <= endpos(1) - startpos(1) + e || y <= endpos(2) - startpos(2) + e
            hx = sqrt(x ^ 2 + (x * k) ^ 2);
            hy = sqrt(y ^ 2 + (y / k) ^ 2);
            re = rem(1 / d * gri * (x + startpos(1) - d) + e, gri);
            if startpos(2) > endpos(2)
                h = -1 / d * y + 1;
            else
                h = 1 / d * y;
            end
            C = [C, round(1 / d * gri * (x + startpos(1) - d) - re + floor(h + 1 / d * startpos(2) + e))];
                if hy > hx 
                    X = [X, hx - prev];
                    prev = hx;
                    x = x + d;
                elseif hx > hy
                    X = [X, hy - prev];
                    prev = hy;
                    y = y + d;
                else
                    X = [X, hx - prev];
                    prev = hx;
                    y = y + d;
                    x = x + d;
                end
        end
    Q=[C; X];
    A = [];
    for j = 1:(gri * gri)
        o = 0;
            for k = 1:length(C)
                if C(k) == j
                    A = [A, Q(2,k)];
                    o = 1;
                end
            end
            if o == 0
                A = [A, 0];
            end
    end
end