function [H_x] = image_entropy(I)
% Calculating image entropy by 256 quantification
I = I/max(max(I));
I = floor(I*255);
imagesc(I)
[C,L]=size(I); 
Img_size=C*L; 
G=256; 
H_x=0;
nk=zeros(G,1);
for i=1:C
for j=1:L
Img_level=I(i,j)+1; 
nk(Img_level)=nk(Img_level)+1; 
end
end
for k=1:G  
Ps(k)=nk(k)/Img_size; 
if Ps(k)~=0; 
H_x=-Ps(k)*log2(Ps(k))+H_x; 
end
end

end

