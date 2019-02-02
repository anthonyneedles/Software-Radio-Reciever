function [decoded_msg, y]=Rx(r, rolloff, desired_user)
%Initialization Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load easy.mat
fs = 850000;        %Sampling frequency
fi = 2000000;       %Intermediate frequency
dataBW = 204000;    %FDM user slot bandwidth
frameLength = 2870; %Frame length in symbols
userSymNum = 875;   %Symbols per user
symPer = 6.4e-6;    %Symbol period
srrcWidth = 8;      %Width of tx SRRC pulse shape

%Initialization End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Downconversion Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Carrier Recovery, Dual Coastas Loop
t = (0:1/fs:(length(r)-1)/fs); recievedData = r';
lent = length(r); mu1=.003; mu2=.0003;  %mu1/mu2 found qualitatively
th1=zeros(1,lent); th2=zeros(1,lent);
for k=1:lent-1                      
    zs1=2*recievedData(k)*sin(2*pi*fi*t(k)+th1(k));
    zc1=2*recievedData(k)*cos(2*pi*fi*t(k)+th1(k));
    th1(k + 1) = th1(k) - mu1*zs1*zc1;
    %Uses gradient descent algorithm find th1/th2
    zs2=2*recievedData(k)*sin(2*pi*fi*t(k)+th1(k)+th2(k));
    zc2=2*recievedData(k)*cos(2*pi*fi*t(k)+th1(k)+th2(k));
    th2(k + 1) = th2(k) - mu2*zs2*zc2;
end

% figure(1) 
% subplot(3,1,1), plot(t,th1), title('Output of Costas Loop #1'), ylabel('th1')
% subplot(3,1,2), plot(t,th2), title('output of Costas Loop #2'), ylabel('th2')
% subplot(3,1,3), plot(t,th1+th2), title('Output of Both Costas Loops')
% xlabel('time (s)'), ylabel('th1+th2')

%Demodulation to Baseband
transBW = 3000;                 %Bandwidth of transition region of filter
transHigh = ((dataBW+transBW/2)/2)/(fs/2);  %Normalized upper/lower values
transLow = ((dataBW-transBW/2)/2)/(fs/2);   %+/-1500 Hz from cutoff freq.
fl = 1024;                      %LPF order
ff = [0 transLow transHigh 1];  %Frequency index for LPF 
fa = [1 1 0 0];                 %Magnitude index for LPF
lpf = firpm(fl,ff,fa);          %Generates user slot baseband isolation LPF
demodCos = 2*cos(2*pi*fi*t+th1+th2);%Demodulation cosine w/ tracked th1/th2
demodData = recievedData.*demodCos; %Demod to baseband
demodData = conv(demodData,lpf);    %LPF for baseband isolation
demodData = demodData(fl/2:end);    %Removal of zeros due to conv()

%Downconversion End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matched Filter Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create SRRC based on tx SRRC parameters (included with data set)
M = fs*symPer;  %Oversampling factor is # of clock cyles per symbol
matchedSRRC = srrc(srrcWidth/2, rolloff, M, 0); %SRRC matched to given tx
matchedSRRC = fliplr(matchedSRRC);              %SRRC, then flipped
matchData = conv(demodData,matchedSRRC);        %Data set matched

% figure(2), plotspec2(recievedData, demodData, matchData, 1/fs)

%Matched Filter End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolator Downsampler Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Output power gradient ascent algorithm
l = srrcWidth/2; mu=.15; delta=.07; n = round(length(matchData)/M);
tnow=srrcWidth/2*M+1; tau=0; symbolData=zeros(1,n); 
tausave=zeros(1,n); tausave(1)=tau; tauIndex=0;                         
while tnow<length(matchData)-l*M      
    tauIndex=tauIndex+1;
    symbolData(tauIndex)=interpsinc(matchData,tnow+tau,l,.3);%Value at
    x_plus=interpsinc(matchData,tnow+tau+delta,l,.3);        %chosen tau
    x_minus=interpsinc(matchData,tnow+tau-delta,l,.3); 
    dx=x_plus-x_minus;  %Derivative using value right and left of value            
    tau=tau+mu*dx*symbolData(tauIndex); %Gradient ascent algorithm to           
    tnow=tnow+M; tausave(tauIndex)=tau; %maximize output power      
end

% figure(3)
% subplot(2,1,1), plot(symbolData(1:tauIndex-2),'.') 
% title('constellation diagram'), ylabel('estimated symbol values')
% subplot(2,1,2), plot(tausave(1:tauIndex-2)) 
% ylabel('offset estimates'), xlabel('iterations')

%Interpolator Downsampler End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Decision Device Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Correlate preamble PAM to symbol list to find starting point
preambleText = 'A0Oh well whatever Nevermind';
preamblePAM = letters2pam2(preambleText);
correlated = xcorr(symbolData, preamblePAM); %Correlation to find locations
if(mean(correlated)<0)                       %of preambles in symbol list
    symbolData = -1.*symbolData;
    correlated = xcorr(symbolData, preamblePAM); 
end  %Symbol list flipped pos to neg detection, flips data back
posMax = max(correlated);   %Finds expected size of preamble peaks
corrThresh = posMax*.5;     %Threshold value to find preamble peaks
maxIdx = find(correlated > corrThresh);
preambleStart = maxIdx(1)-length(symbolData); 
%First peak below thresh is first preamble

%Compare expected 3's in preamble with actual values, find average "3"
preamblePAM3idx = find(preamblePAM == 3);
symbolData3Sum = 0;
for maxIndex = 1:length(preamblePAM3idx)-1
    symbolData3Sum = symbolData3Sum + symbolData(preamblePAM3idx(maxIndex)+preambleStart);
end
symbolData3Avg = symbolData3Sum/length(preamblePAM3idx);

%Use found average for expected "3's", normalize to 3 for gain control
gain = 3/symbolData3Avg;
symbolDataNorm = (gain).*symbolData;   

%Quantize normalized symbol set to 3-PAM
symbolQuant = quantalph(symbolDataNorm,[-3,-1,1,3]); 

% figure(4)
% subplot(3,1,1), plot(symbolData, '.'), title('Prenormalized Symbol List')
% subplot(3,1,2), plot(symbolDataNorm, '.'), title('Normalized Symbol List')
% subplot(3,1,3), plot(symbolQuant, '.'), title('Quantized Symbol List')

%Decision Device End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Decoder Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Decode quantized symbols for User #(desired_user)
frameNum = floor(length(symbolQuant)/frameLength);
userPAM = zeros(1, userSymNum*frameNum);
startPoint = preambleStart+length(preamblePAM)+(desired_user-1)*userSymNum+1;
for curIndex = 0:frameNum-1
    userPAM(curIndex*userSymNum+1:(curIndex+1)*userSymNum) = symbolQuant(startPoint+curIndex*frameLength:startPoint+curIndex*frameLength+userSymNum-1);
end
userMsg = pam2letters2(userPAM)

%Decoder End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output Declaration Start%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = symbolDataNorm(preambleStart:end);
decoded_msg = userMsg;
