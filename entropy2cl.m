function entropy2cl(sWqM127,sWqM255,sWqM511,sWqM1023,sWqM2047)
% variation of the image entropy with the code length
L = [127,255,511,1023,2047];
E = zeros([1,length(L)]);
E(1) = image_entropy(abs(sWqM127));
E(2) = image_entropy(abs(sWqM255));
E(3) = image_entropy(abs(sWqM511));
E(4) = image_entropy(abs(sWqM1023));
E(5) = image_entropy(abs(sWqM2047));
figure;
plot(L,E,'o');
hold on
plot(linspace(127,2047,2047-127+1),interp1(L,E,linspace(127,2047,2047-127+1),'makima'));
xlabel('Code length');
ylabel('Image entropy');

end