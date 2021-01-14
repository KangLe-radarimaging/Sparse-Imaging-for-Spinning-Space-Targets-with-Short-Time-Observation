function [sn]=PskLiang(an,fc,Tb,Tspre,Tfpos,Td,Tv,tkM)
% Using Gold code to generate the orthogonal signal
%sn: echo signal
%sc: reference signal
%an: Gold code
%fc: carrier frequency
%Tb: code width
%Tspre: time before transmitting signal
%Tfpos: time after transmitting signal
%Td: time delay
%Tv: target speed
%tkM: oversampling rate

an1=an*2-1;

cj=sqrt(-1);
Bb=1/Tb;
Tl=Tb*length(an);
tks=0+Tspre;
tkf=Tl+Tfpos;
dtk=1/(tkM*Bb);
tk=tks:dtk:tkf;

tkd=tk-Td-Tv*tk;
tkdM=ones(length(an1),1)*tkd;
TbF=Tb*(0:(length(an1)-1))'*ones(1,length(tk));
TbB=TbF+Tb;
an1M=an1'*ones(1,length(tk));
sn=ones(1,length(an1))*(an1M.*exp(cj*2*pi*fc*tkdM).*(tkdM>TbF&tkdM<=TbB));
% figure
% plot(abs(fft(sn)));
% for i=1:length(an)-0%
%     %snb=(In(i)*cos(pi*tk/2/Tb).*cos(2*pi*fc*tk)-Pn(i)*sin(pi*tk/2/Tb).*sin(2*pi*fc*tk)).*(tk>(i-1)*Tb & tk<=(i)*Tb)+snb;
%     %snby=(Qn(i+0)*cos(pi*tk/2/Tb).*cos(2*pi*fc*tk)+In(i+0)*sin(pi*tk/2/Tb).*sin(2*pi*fc*tk)).*(tk>(i-1)*Tb & tk<=(i)*Tb)+snby;%经过分析，是正确的
%     snby=an1(i).*exp(cj*2*pi*fc*tkd).*(tkd>(i-1)*Tb & tkd<=(i)*Tb)+snby;%经过分析，是正确的
%     %sc  =((Qn(i+0)*cos(pi*tk /2/Tb)-cj*In(i+0)*sin(pi*tk /2/Tb)).*exp(cj*2*pi*fc*tk )).*(tk >(i-1)*Tb & tk <=(i)*Tb)+sc ;
% end






