function res = hz2mel(hz)
    res = 2595*log10(1 + hz/700);
end