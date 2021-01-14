function fs=fty(s);
%% fft for row
fs=fftshift(fft(fftshift(s.'))).';
end

