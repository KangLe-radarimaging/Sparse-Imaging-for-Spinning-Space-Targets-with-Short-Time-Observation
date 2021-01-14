function s=iftx(fs)
%% ifft for column
s=fftshift(ifft(fftshift(fs)));
end

