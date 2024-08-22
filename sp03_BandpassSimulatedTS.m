function [ts, Kuramoto] = sp03_BandpassSimulatedTS(TS, low_f, high_f, TR)

N = size(TS,1);
T = size(TS,2);
Ts = T*TR;
freq = (0:T/TR-1)/Ts;

delt = TR; % sampling interval
k=2; % 2nd order butterworth filter
fnq=1/(2*delt); % Nyquist frequency
flp = low_f; % lowpass
fhi = high_f; % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn); % construct the filter

for seed = 1:N
    signaldata = TS(seed,:);
    x = sp04_demean(detrend(signaldata));
    ts_dummy = filtfilt(bfilt,afilt,x);
    ts(seed,:) = ts_dummy;
    H(seed,:) = hilbert(sp04_demean(ts_dummy));
end
Kuramoto = abs(sum(complex(cos(angle(H)),sin(angle(H))))/N);