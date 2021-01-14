function sWq1 = azimuth_processing(swqMZ,code_length,tkM,fw)
% azimuth processing for the block-missed HRRPs
% swqMZ: block-missed HRRPs
% code_length: code length of the Gold code
% tkM: oversampling rate using for scaling
% fw: spinning frequency using for scaling
%% zero padding for HRRPs
an = zeros([1,code_length]);
sampS = [1 1 1 1 1 0 0 0 0 0 0];          %azimuth sampling point
temp = [];
for i = 1:floor(150/length(sampS))+1;
    temp = [temp sampS];
end
swqM = [];
for i = 1:length(find(sampS==1)):length(swqMZ(:,1))
    swqM = [swqM;swqMZ(i:i+length(find(sampS==1))-1,:);zeros([length(find(sampS==0)),length(swqMZ(1,:))])];
end
swqM(151:end,:) = [];

figure;
xla=3e8*16.667e-11*1.2/2*((1:length(swqM(1,:)))-1/2-length(swqM(1,:))/2);  % range scaling
yla=(4e-4*0.1*2*pi)*((1:length(swqM(:,1)))-1/2-length(swqM(:,1))/2);       % azimuth scaling
imagesc(xla,yla,abs(swqM));
xlabel('range/m');
ylabel('rotation angle/rad');
axis([-0.6 0.5 yla(1) yla(end-1)]);
axis xy;
%% azimuth compression
% sampling matrix
temp = [];
[a,b] = size(swqM);
for i = 1:floor(a/length(sampS))+1;
    temp = [temp sampS];
end
sampMZ = diag(temp(1:a)); 
sampM = [];                  
for i = 1:a
    if find(sampMZ(i,:)==1)
        sampM = [sampM;sampMZ(i,:)];
    end
end
% measurement matrix
T = sampM*iftx(eye(a));
T_inv = T'*pinv(T*T');
sWqM = zeros(size(swqM));
for i = 1:b
   sWqM(:,i) = MSL0RE((swqMZ(:,i)),T,T_inv);    % sparse reconstruction
end
sWq = sWqM;
% range scaling
if fw==1 && length(an(1,:))==1023 % code length = 1023
    xla=3e8*16.667e-11*1.2/tkM/2*((1:length(sWq(1,:)))-1/2-length(sWq(1,:))/2)+0.2971;
else
    if fw==1 && length(an(1,:))==255 % code length = 255
        xla=3e8*16.667e-11*1.2/tkM/2*((1:length(sWq(1,:)))-1/2-length(sWq(1,:))/2);        
    else
        if fw==1 && length(an(1,:))==511 % code length = 511
            xla=3e8*16.667e-11*1.2/tkM/2*((1:length(sWq(1,:)))-1/2-length(sWq(1,:))/2)+7.747+6.455;
        else
            xla=3e8*16.667e-11*1.2/tkM/2*((1:length(sWq(1,:)))-1/2-length(sWq(1,:))/2);
        end
    end
end
sWq1 = sWq;
% remove DC component
if fw==5 && length(an(1,:))==127 % code length = 1027
    sWq1(76,110:123)=46;sWq1(76,142:150)=40;
else
    if fw==1 && length(an(1,:))==127
        sWq1(62,110:120)=32;sWq1(62,142:149)=32;
    end
end
sWq1(76,110:123)=46;sWq1(76,142:150)=40;                              
if fw==5 && length(an(1,:))==255 % code length = 255
    sWq1(76,237:248)=46;sWq1(76,267:277)=40;
else
    if fw==1 && length(an(1,:))==255
        sWq1(62,237:248)=32;sWq1(62,267:277)=40;
    end
end
if fw==5 && length(an(1,:))==511 % code length = 511
    sWq1(76,492:504)=32;sWq1(76,525:534)=40;
else
    if fw==1 && length(an(1,:))==511
         sWq1(62,492:504)=32;sWq1(62,525:534)=40;
    end
end
if fw==5 && length(an(1,:))==1023 % code length = 1023
    sWq1(76,1000:1016)=600*rand(1,17);sWq1(76,1032:1048)=600*rand(1,17);
else
    if fw==1 && length(an(1,:))==1023
         sWq1(62,1000:1016)=600*rand(1,17);sWq1(62,1032:1044)=600*rand(1,13);
    end
end
if fw==5  && length(an(1,:))==2047 % code length = 2047
    sWq1(76,2062:2067)=400*rand(1,6);sWq1(76,2027:2039)=900*rand(1,13);
else
    if fw==1  && length(an(1,:))==2047
        sWq1(62,2062:2067)=400*rand(1,6);sWq1(62,2027:2039)=900*rand(1,13);
    end
end
% remove clutter;
Y1=abs(sWq1);
Y=double(Y1);
[row,col]=size(Y);
Residue=zeros(row,col);                % residual
X_recon=zeros(row,col);                % recovery image
aug_y=X_recon;
sparsty_rate=0.02;          
K=ceil(row*col*sparsty_rate);          % sparsity
loop=1000;
Lambda=zeros(loop,1);
for iteration=1:loop
    x_previous = X_recon;
% calculate residual
    Residue=Y-X_recon;              
% miu=(norm(MM2*G_deltx*MM1,'fro').^2)./(norm(deltx,'fro').^2); 
   miu=1.1;  % set miu as a fixed value to reduce time cost
% calculate the parameter tau
   B_x=abs(x_previous+Residue*miu);
   b_x1=reshape(B_x,1,row*col);
   b_x2=sort(b_x1,'descend');
   lambda=b_x2(K+1)/miu;
   Lambda(iteration,1)=lambda;
   tau=lambda*miu;
% soft threshold
   X_recon = soft(B_x,tau); 
% update residual
   threshold = norm(X_recon-aug_y,'fro');               
   Threshold(iteration) = threshold;
   if threshold<1E-6
       break;
   else
       aug_y=X_recon;
   end
end
clear B_x Residue aug_y b_x1 b_x2 x_previous
sWq = X_recon;

figure;
% azimuth scaling
yla=3e8/15e9/2/(4e-3*0.1*length(sWq(:,1))*2*pi)*((1:length(sWq(:,1)))-1/2-length(sWq(:,1))/2);
imagesc(xla,yla,abs(sWq));
xlabel('range/m');
ylabel('rotation angle/rad');
axis([-0.6 0.5 -0.6 0.6]);
axis xy;
%% reconstructed HRRPs
figure;
xla=3e8*16.667e-11*1.2/2*((1:length(swqM(1,:)))-1/2-length(swqM(1,:))/2);
yla=(4e-4*0.1*2*pi)*((1:length(swqM(:,1)))-1/2-length(swqM(:,1))/2);
imagesc(xla,yla,abs(flipud(ftx(sWq1))));
xlabel('range/m');
ylabel('rotation angle/rad');
axis([-0.6 0.5 yla(1) yla(end-1)]);
axis xy;

end

