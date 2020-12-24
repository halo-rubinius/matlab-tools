% reference: https://rdrr.io/rforge/labeling/src/R/labeling.R
function ticks = fineticks(dmin, dmax, nticks, varargin)
% ticks = fineticks(dmin, dmax, nticks, varargin)

tmin = dmin;
tmax = dmax;

algorithm = 'extended';

switch algorithm
    case 'heckbert'
        ticks = heckbert(tmin, tmax, nticks);
    case 'wilkinson'
        if isempty(varargin)
            ticks = wilkinson(tmin, tmax, nticks);
        else
            ticks = wilkinson(tmin, tmax, nticks, varargin{:});
        end
    otherwise
        if isempty(varargin)
            ticks = extended(tmin, tmax, nticks);
        else
            ticks = extended(tmin, tmax, nticks, varargin{:});
        end
end

end


function d = heckbert(dmin, dmax, m)

range = heckbert_nicenum((dmax-dmin), false);
lstep = heckbert_nicenum(range/(m - 1), true);
lmin = floor(dmin/lstep) * lstep;
lmax = ceil(dmax/lstep) * lstep;
d = lmin:lstep:lmax;

end


function num = heckbert_nicenum(x, round)

e = floor(log10(x));
f = x / 10^e;

if round
    if f < 1.5
        nf = 1;
    elseif f < 3
        nf = 2;
    elseif f < 7
        nf = 5;
    else
        nf = 10;
    end
else
    if f <= 1
        nf = 1;
    elseif f <= 2
        nf = 2;
    elseif f <= 5
        nf = 5;
    else
        nf = 10;
    end
end

num = nf * 10^e;

end


function d = wilkinson(dmin, dmax, m, varargin)

if nargin > 3
    Q = varargin{1};
else
    Q = [1, 5, 2, 2.5, 3, 4, 1.5, 7, 6, 8, 9];
end

if nargin > 4
    mincoverage = varargin{2};
else
    mincoverage = 0.8;
end

if nargin > 5
    mrange = varargin{3};
else
    mrange = max(floor(m / 2), 2):ceil(6*m);
end

best = [];
for i = 1:length(mrange)
    k = mrange(i);
    result = wilkinson_nicescale(dmin, dmax, k, Q, mincoverage, mrange, m);
    if ~isempty(result) && (isempty(best) ...
            || result.score > best.score)
        best = result;
    end
end

d = best.lmin:best.lstep:best.lmax;

end


function best = wilkinson_nicescale(dmin, dmax, k, varargin)

if nargin > 3
    Q = varargin{1};
else
    Q = [1, 5, 2, 2.5, 3, 4, 1.5, 7, 6, 8, 9];
end

if nargin > 4
    mincoverage = varargin{2};
else
    mincoverage = 0.8;
end

if nargin > 5
    mrange = varargin{3};
else
    mrange = [];
end

if nargin > 6
    m = varargin{4};
else
    m = k;
end

Q = [10, Q];

range = dmax - dmin;
intervals = k - 1;
granularity = 1 - abs(k-m) / m;

delta = range / intervals;
base = floor(log10(delta));
dbase = 10^base;

best = [];
for i = 1:length(Q)
    tdelta = Q(i) * dbase;
    tmin = floor(dmin/tdelta) * tdelta;
    tmax = tmin + intervals * tdelta;

    if tmin <= dmin && tmax >= dmax
        if tmin <= 0 && tmax >= 0
            roundness = 1 - ((i - 1) - 1) / length(Q);
        else
            roundness = 1 - ((i - 1) - 0) / length(Q);
        end
        coverage = (dmax - dmin) / (tmax - tmin);
        if coverage > mincoverage
            tnice = granularity + roundness + coverage;
            if isempty(best) || tnice > best.score
                best.lmin = tmin;
                best.lmax = tmax;
                best.lstep = tdelta;
                best.score = tnice;
            end
        end
    end
end

end


function r = floored_mod(a, n)

r = a - n * floor(a/n);

end


function d = simplicity(q, Q, j, lmin, lmax, lstep)

eps = 1.e-10;

n = length(Q);
i = find(Q == q);
i = i(1);
if (floored_mod(lmin, lstep) < eps || lstep - floored_mod(lmin, lstep) < eps) ...
        && lmin <= 0 && lmax >= 0
    v = 1;
else
    v = 0;
end

d = 1 - (i - 1) / (n - 1) - j + v;

end


function r = simplicity_max(q, Q, j)

n = length(Q);
i = find(Q == q);
i = i(1);
v = 1;
r = 1 - (i - 1) / (n - 1) - j + v;

end


function r = coverage(dmin, dmax, lmin, lmax)

range = dmax - dmin;
r = 1 - 0.5 * ((dmax - lmax)^2 + (dmin - lmin)^2) / ((0.1 * range)^2);

end


function r = coverage_max(dmin, dmax, span)

range = dmax - dmin;
if span > range
    half = (span - range) / 2;
    r = 1 - 0.5 * (half^2 + half^2) / ((0.1 * range)^2);
else
    r = 1;
end

end


function d = density(k, m, dmin, dmax, lmin, lmax)

r = (k - 1) / (lmax - lmin);
rt = (m - 1) / (max([lmax, dmax]) - min([dmin, lmin]));
d = 2 - max([r / rt, rt / r]);

end


function d = density_max(k, m)

if k >= m
    d = 2 - (k - 1) / (m - 1);
else
    d = 1;
end

end


function r = legibility(lmin, lmax, lstep)
r = 1;
end


function dd = extended(dmin, dmax, m, varargin)

if nargin > 3
    loose = varargin{1};
else
    loose = false;
end

if nargin > 4
    Q = varargin{2};
else
    Q = [1, 5, 2, 2.5, 4, 3];
end

if nargin > 5
    w = varargin{3};
else
    w = [0.25, 0.2, 0.5, 0.05];
end

n = length(Q);
best.score = -2;

j = 1;
while (j < Inf)
    for iq = 1:length(Q)
        q = Q(iq);
        sm = simplicity_max(q, Q, j);
        if w(1) * sm + w(2) + w(3) + w(4) < best.score
            j = Inf;
            break;
        end

        k = 2;
        while (k < Inf)
            dm = density_max(k, m);
            if w(1) * sm + w(2) + w(3) * dm + w(4) < best.score
                break;
            end
            delta = (dmax - dmin) / (k + 1) / j / q;
            z = ceil(log10(delta));
            while (z < Inf)
                step = j * q * 10^z;
                cm = coverage_max(dmin, dmax, step*(k - 1));

                if w(1) * sm + +w(2) * cm + w(3) * dm + w(4) < best.score
                    break;
                end

                min_start = floor(dmax/(step)) * j - (k - 1) * j;
                max_start = ceil(dmin/step) * j;

                if min_start > max_start
                    z = z + 1;
                    continue;
                end

                for start = min_start:max_start
                    lmin = start * (step / j);
                    lmax = lmin + step * (k - 1);
                    lstep = step;

                    s = simplicity(q, Q, j, lmin, lmax, lstep);
                    c = coverage(dmin, dmax, lmin, lmax);
                    g = density(k, m, dmin, dmax, lmin, lmax);
                    l = legibility(lmin, lmax, lstep);

                    score = w(1) * s + w(2) * c + w(3) * g + w(4) * l;
                    if score > best.score && (~loose || (lmin <= dmin && lmax >= dmax))
                        best.lmin = lmin;
                        best.lmax = lmax;
                        best.lstep = lstep;
                        best.score = score;
                    end

                end
                z = z + 1;
            end
            k = k + 1;
        end
    end
    j = j + 1;
end

dd = best.lmin:best.lstep:best.lmax;

end
