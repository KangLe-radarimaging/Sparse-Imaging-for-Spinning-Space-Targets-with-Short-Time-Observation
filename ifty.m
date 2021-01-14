function s=ifty(fs)
%% ifft for row
s=fftshift(ifft(fftshift(fs.'))).';
end

