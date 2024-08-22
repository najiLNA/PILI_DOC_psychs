function [ts, Power_Areas, freq, Kuramoto, Phases] = sp01_GetSpectrum(TS, low_f, high_f, TR)

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
    ts_dummy_z = zscore(ts_dummy);
    pw = abs(fft(ts_dummy_z));
    PowSpect = pw(1:floor(T/TR)).^2/(T/TR);
    Power_Areas(:,seed)=sp05_gaussfilt(freq,PowSpect,0.01);
    H(seed,:) = hilbert(sp04_demean(ts_dummy));
    Phases(seed,:) = angle(H(seed,:));
end
Kuramoto = abs(sum(complex(cos(angle(H)),sin(angle(H))))/N);