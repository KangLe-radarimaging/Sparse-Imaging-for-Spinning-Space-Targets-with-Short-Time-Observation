clear all;
close all;
clc;
%% parameters
ith=1;
for i=-0:5 %generate different Gold sequences in the same family  
    goldseq = comm.GoldSequence('FirstPolynomial',[1 0 1 1 1 0 0 0 1],...
        'SecondPolynomial',[1 1 1 0 0 1 1 1 1],...
        'FirstInitialConditions',[0 0 0 0 0 0 0 1],...
        'SecondInitialConditions',[0 0 0 0 0 0 0 1],...
        'Index',i,'SamplesPerFrame',255);
    an(ith,:) = (goldseq())';
    ith=ith+1;
end
fc=15e9;% carrier frequency
c0=3e8;% light speed
Tb=1.5e-9;% code width
tkM=1;% oversampling rate
xn=[0 0 0 -1 1];% scattering position，intensity is 1
yn=[-1 0 1 0 0];% scattering position
fw=5; % spinning frequency
dtm=5e-05;% pulse interval
thm=14.4;% azimuth angle
tm=0:dtm:0.0012;%slow time
%% -------------------------SIMULATION-------------------------------
%% echo generation
n_tmj=round(128/length(tm));% number of slow time splitting intervals
if(n_tmj<=1) % observation time is enough
    dtmj=1;
    tmj=0;
else
    dtmj=dtm/n_tmj;% sampling interval of code division signal
    tmj=0:dtmj:dtm;% sampling time series of code division signal 
end

Tl=length(an)*Tb;
Tspre=-Tl/2;Tfpos=Tl/2;
for ia=1:length(tm)
    Tv=0;%translation compensation is completed       
    for i=1:length(tmj)%calculate the echo of all orthogonal codes
        tmz=tm(ia)+tmj(i);%equivalent slow time
        thz=2*pi*fw*tmz;%equivalent rotation angle
        xns=xn*cos(thz)-yn*sin(thz);%distance in LOS
        Td=xns(1)*2/c0;% delay time
        %calculate the echo of the scattering centers
        sn(i,:)=PskLiang(an(i,:),fc,Tb,Tspre,Tfpos,Td,Tv,tkM);
        for j=2:length(xn)
            Td=xns(j)*2/c0;
            sn(i,:)=PskLiang(an(i,:),fc,Tb,Tspre,Tfpos,Td,Tv,tkM)+sn(i,:);       
        end
    end
    snz(ia,:)=ones(1,length(tmj))*sn;% Add the echo signals with same slow time together
end
snz = awgn(snz,5);
%% imaging processing
for i=1:length(tmj)
        sc(i,:)=PskLiang(an(i,:),fc,Tb,Tspre,Tfpos,0,Tv,tkM); % reference signal before rearrangement
end
snq=[];scq=[];

for i=1:length(tm)% rearrange the reference signal to tmq(1:length(tm)*length(tmj))
    snq=[snq;ones(length(tmj),1)*snz(i,:)];
    scq=[scq;sc];
end
%% conventional ISAR imaging method
% zero padding
snqCon = zeros(size(snq));
scqCon = zeros(size(scq));
snqCon(1:13:end,:)= snq(1:13:end,:);
scqCon(1:13:end,:)= scq(1:13:end,:);
% range compression
swqCon=ifty(fty(snqCon).*conj(fty(scqCon)));
figure(1);
imagesc(abs(swqCon)); 
xlabel('range');
ylabel('azimuth');
% title('HRRPs');
axis xy;
% azimuth compression
sWqCon=ftx(swqCon);
figure(2);
imagesc(abs(sWqCon));
xlabel('range');
ylabel('azimuth');
% title('Imaging result');
axis xy;

figure(3);% range profile
plot(abs(sWqCon(76,:)/max(abs(sWqCon(76,:)))))
xlabel('range');
% title('Range profile');
axis([1 511 0 1]);

figure(4)% azimuth profile
plot(abs(sWqCon(1:143,255)/max(abs(sWqCon(1:143,255)))))
xlabel('azimuth');
% title('Azimuth profile')
%% proposed method
% sampling matrix
sampS = [1 1 1 1 1 0 0 0 0 0 0 0 0];%azimuth sampling point
temp = [];
[a,b] = size(snq);
for i = 1:floor(a/length(sampS))+1;
    temp = [temp sampS];
end
FTr = ftx(eye(a));
sampMZ = diag(temp(1:a));% zero-padded sampling matrix
sampM = []; % sampling matrix
for i = 1:a
    if find(sampMZ(i,:)==1)
        sampM = [sampM;sampMZ(i,:)];
    end
end
% range compression
swq=ifty(fty(sampM*scq).*conj(fty(sampM*snq))); %no zero-padding
swqZ=ifty(fty(sampMZ*scq).*conj(fty(sampMZ*snq)));%zero-padding
figure(5);
imagesc(abs(swqZ(1:end,:)));
xlabel('range');
ylabel('azimuth');
% title('HRRPs');
axis xy;
% azimuth compression
T = sampM*iftx(eye(a));% measurement matrix
T_inv = T'*pinv(T*T');
sWq = zeros(size(snq));
for i = 1:b
   sWq(:,i) = MSL0RE((swq(:,i)),T,T_inv);
end
figure(6);
imagesc(abs(sWq));
xlabel('range');
ylabel('azimuth');
% title('Imaging result');
axis xy;%二维成像结果

figure(7)
plot(abs(sWq(76,:)/max(abs(sWq(76,:)))))
xlabel('range');
axis([1 511 0 1]);
figure(8)
plot(abs(sWq(:,256)/max(abs(sWq(:,256)))))
xlabel('azimuth');

%% -------------------------EXPERIMENT-------------------------------
%% conventional ISAR imaging method
%% load the experiment data by transmitting the LFM signal
load swqLFM 

%%  zero padding imaging
swqLFMZ = [];
for i = 1:length(swqLFM(:,1))
    swqLFMZ = [swqLFMZ;zeros([9,length(swqLFM(1,:))]);swqLFM(i,:);];% Zeros padded high-resolution range profiles
end
figure; 
xla=3e8/6.0e9/2*((1:length(swqLFMZ(1,:)))-1/2-length(swqLFMZ(1,:))/2)-0.05;%range scaling
yla=linspace(-1.1/180*pi,1.1/180*pi,length(swqLFMZ(:,1)));%azimuth scaling
imagesc(xla,yla,abs(swqLFMZ(:,end:-1:1)));
axis([[-0.6 0.5]  yla(1) yla(end)]);
xlabel('range/m');
ylabel('rotation angle/rad');
axis xy;

% Azimuth compression
sWqLFMZ = ftx(swqLFMZ);
figure; 
xla=3e8/6.0e9/2*((1:length(sWqLFMZ(1,:)))-1/2-length(sWqLFMZ(1,:))/2)-0.14;%range scaling
yla=3e8/15e9/2/(4e-3*0.1*length(sWqLFMZ(:,1))*2*pi)*((1:length(sWqLFMZ(:,1)))-1/2-length(sWqLFMZ(:,1))/2);%azimuth scaling
imagesc(xla,yla,abs(sWqLFMZ(end:-1:1,end:-1:1)));
axis([[-0.6 0.5] -0.6 0.6]);
xlabel('range/m');
ylabel('azimuth/m');
axis xy;

%% Interpolation imaging
swqLFMI = interp1(1:12,swqLFM,1:11/119:12);%Interpolated high-resolution range profiles;
figure;
xla=3e8/6.0e9/2*((1:length(swqLFMI(1,:)))-1/2-length(swqLFMI(1,:))/2)-0.05;%range scaling
yla=linspace(-1.1/180*pi,1.1/180*pi,length(swqLFMI(:,1)));%azimuth scaling
imagesc(xla,yla,abs(swqLFMI(:,end:-1:1)));
axis([[-0.6 0.5]  yla(1) yla(end)]);
xlabel('range/m');
ylabel('rotation angle/rad');
axis xy;

%Azimuth compression
sWqLFMI = ftx(swqLFMI);
figure;
xla=3e8/6.0e9/2*((1:length(sWqLFMI(1,:)))-1/2-length(sWqLFMI(1,:))/2)-0.14;%range scaling
yla=3e8/15e9/2/(4e-3*0.1*length(sWqLFMI(:,1))*2*pi)*((1:length(sWqLFMI(:,1)))-1/2-length(sWqLFMI(:,1))/2);%azimuth scaling
imagesc(xla,yla,abs(sWqLFMI(end:-1:1,end:-1:1)));
axis([[-0.6 0.5] -0.6 0.6]);
xlabel('range/m');
ylabel('azimuth/m');
axis xy;
%% proposed method
% load the experiment data by transmitting the orthogonal coding signal consisted
% of 5 Gold codes with the code length 127, 255, 511, 1023, and 2047
load swqGoldZ

% azimuth compression
sWqM127 = azimuth_processing(swqM127Z,127,1,5);
sWqM255 = azimuth_processing(swqM255Z,255,1,5);
sWqM511 = azimuth_processing(swqM511Z,511,1,5);
sWqM1023 = azimuth_processing(swqM1023Z,1023,1,5);
sWqM2047 = azimuth_processing(swqM2047Z, 2047,1,5);

%%  variation of the image entropy with the code length
% entropy2cl(sWqM127,sWqM255,sWqM511,sWqM1023,sWqM2047);














