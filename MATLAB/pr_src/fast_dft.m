function [corrmap, maxval] = fast_dft(Mq,Mi)
    Fq = fft(Mq); % fft along theta axis
    Fn = fft(Mi);
    corrmap_2d = ifft(Fq.*conj(Fn));

    corrmap = sum(corrmap_2d,2);
    maxval = max(corrmap);
end
