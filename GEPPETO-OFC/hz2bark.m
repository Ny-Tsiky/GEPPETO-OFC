function res = hz2bark(f)
res = (26.81*f)./(1960+f) - 0.53;
lowmask = res < 2;
highmask = res > 20.1;
res(lowmask) = res(lowmask) + 0.15*(2-res(lowmask));
res(highmask) = res(highmask) + 0.22*(res(highmask)-20.1);
end
