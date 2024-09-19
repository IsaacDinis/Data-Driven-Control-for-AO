% To test linear prediction like presented in OEF paper by Guyon & Males 2017, 
% or the extension in Males & Guyon 2018
%
%%
% ccc
%% Generate a spectrum and random temporal sequence
Nsample = 2e4;
fmax    = 4000;
fos     = 20; % Frequency oversampling to generate DIT fitlered data
seeing  = 0.75; % ["]
lbd     = 500e-9;
D       = 8;
dt      = 2; % Delay in frame

% inputType = 'compass';
% inputType = 'vibrations';
inputType = 'atm';

% dfolder = fullfile(oomaoFolders.res, 'AR_filter_tests');
% genUtilities.createFolder(dfolder)
compass_tilt = fitsread("P.fits");

%%
% TODO -- Add sensor integration  moving average box to realistically
% filter signal
r0 = (500e-9/(seeing*4.85e-6));
TT_Noll_var = 0.448*(D/r0)^(5./3.)*(500e-9/(2*pi))^2; % TT variance in radian, per axis


wn  = ifft(randn(1,Nsample)); % White noise
ttt = (1:Nsample)/fmax;

f   = linspace(fmax/Nsample, fmax, Nsample);
df = f(2)-f(1);

% Ddefine 2 vibration frequencies
f1 = 250;
avib = 0.025; % [lbd/D RMS]
stvib500 = 4*avib*lbd*cos(2*pi*f1*ttt);

f2 = 55;
avib = 0.025; % [lbd/D RMS]
stvib50 = 4*avib*lbd*cos(2*pi*f2*ttt);


switch inputType
    case 'atm'
        L0 = 25;
        
        f_   = linspace(fmax/Nsample, fmax*fos, Nsample*fos);
        sp_ = f_.^-(11./3);
        dit_filter = zeros(1,fos*Nsample);
        dit_filter(1:fos)=1./fos;
        dit_filter = ifft(dit_filter);
        sp_= sp_.*dit_filter;
        sp = sp_(1:Nsample);

        fcut = 400; % Cut spectrum at high frequencies for display purpose
        [~,icut] = min(abs(f-fcut));
        % spcut(icut:end) = 0*1e-9;
        
        st_atm = real(ifft(sqrt(sp).*wn,Nsample))*sqrt(Nsample);
        psd_amp_corr = TT_Noll_var/var(st_atm);
        st_atm = st_atm*sqrt(psd_amp_corr);
        % st_atm = st_atm*sqrt(psd_amp_corr);
        sp     = sp*psd_amp_corr;
        % st_atm = real(ifft(sqrt(sp).*wn,Nsample))*sqrt(Nsample);

        % Add the vibrations
        st0 = st_atm;
        st1 = st_atm+stvib500;
        st2 = st_atm+stvib50;


    case 'vibrations'
        % Generate a series of vibrations like in Guyon 2018
        sp = f*1e-30;

        Nvib = 20;
        fvib = linspace(1,300,Nvib);
        st_atm = sum(avib*lbd*cos(2*pi*fvib'*ttt+pi*randn(Nvib,1)),1);
        st0  = st_atm;
        st1 = st_atm+stvib500;
        st2 = st_atm+stvib50;
     case 'compass'
        % Generate a series of vibrations like in Guyon 2018
        sp = f*1e-30;

        Nvib = 20;
        fvib = linspace(1,300,Nvib);
        st_atm = compass_tilt;
        st0  = st_atm;
        st1 = st_atm;%+stvib500;
        st2 = st_atm;%+stvib50;
    otherwise
        error('Nope');
end

% std(st1)
% std(stvib500)
% std(stvib50)
% 

figure()
subplot(211)
plot(ttt, st1, '-')
hold on
plot(ttt, st2, '-')
hold off
f_ = linspace(min(f(:)), max(f(:))/2,Nsample/2);
sp_ = abs(fft(st1)).^2;
sp_ = sp_(1:Nsample/2);

subplot(212) % Check scaling is OK!
loglog(f_, abs(sp_))
hold on
loglog(f, abs(sp))
ylim([min(abs(sp)),max(abs(sp))])

%%-------------------------------------------------------------------------
%% Levinson method from Males -- No noise filtering?
% n = 20000; % Number of point to build the filter
% 
% stn = st1+randn(1,Nsample)*std(st1)*0.0;
% 
% [r,lg] = xcorr((stn(1:n)),'biased');
% % r = abs(fftshift(fft(f)))*Nsample;
% 
% % R=zeros(n);
% % for i=1:n
% %     R(i,:) = circshift(r(1:n),i-1);
% % end
% 
% figure(), loglog(abs(fft(r))), xlim([1, n/2])
% % figure, imagesc(R)
% 
% r(lg<0) = [];
% 
% 
% 
% filter_order = 300;
% % filter_order = n-1;
% % filter_order = 800;
% [a,e,k] = levinson(r,filter_order);
% rp = a(1:end-1).*(1:filter_order);
% 
% % pxx = pyulear(st(1:n),filter_order);
% % figure, loglog(pxx), hold on, loglog(sp_)
% 
% figure, plot(a)
% figure, hold on, plot(r)
% figure, plot(filter_order)
% figure, plot(e)
% 
% %%
% af = zeros(size(a)) + a;
% conf = sqrt(2)*erfinv(0.99)/sqrt(n);
% validDelay = abs(k)>=conf;
% sum(validDelay)
% st_=[];
% 
% % validDelay = 1:filter_order<3;
% 
% af(~validDelay)=0;
% 
% k_=k;
% k_(~validDelay)=nan;
% 
% figure
% stem(k,'o')
% hold on
% stem(k_,'o','filled')
% [X,Y] = ndgrid(xlim,conf*[-1 1]);
% plot(X,Y,'--r')
% % xlim([1,100])
% %%
% for i=1:(Nsample-n)
%     st_(i) = (af)*(stn(i:i+filter_order))';
% end
% 
% % st_=exp(st_);
% 
% 
% %% 
% [xc,lc]=xcorr(st1,st_);
% % figure(), plot(lc, xc)
% [~,il] = max(xc);
% lag=lc(il);
% lag=filter_order;
% 
% % st_=exp(st_);
% 
% 
% 
% [xc,lc]=xcorr(st1,st_);
% figure(), plot(lc, xc)
% [~,il] = max(xc);
% lag=lc(il);
% lag=filter_order;
% 
% figure()
% subplot(211)
% plot((lag:+lag+Nsample-n-1), st_/std(st_), '-')
% hold on
% plot((st1/std(st1)))
% plot((stn/std(stn)))
% % plot(st_/std(st_))
% hold off
% % xlim([1,100*lag])
% 
% psd_ = abs(fft(st_));
% psd_ = psd_(1:numel(psd_));
% subplot(212)
% loglog(psd_)



%% SVD method with noise fitlering -- Guyon's description
m = 1;     % Number of mdoes (1 for know)
n = 2;   % Filter order; order 1 = standard integrator; order 3 == PID
k = 10000; % Size of training set

noise = 1e-2*1; % White noise amplitude in temporal space, in unit of 
             % input disturbance sandard deviation

% Add measurement noise to training set
% st_training = st1;
st_training = st1+randn(1,Nsample)*std(st_atm)*noise;

% Create noisy data sets
st1n = st1 + randn(1,Nsample)*std(st_atm)*noise;
st2n = st2 + randn(1,Nsample)*std(st_atm)*noise;


D = zeros(n*m,k);
tic
for i=1:k
    h = (st_training(i:i+n-1)); % Data history to predict next step
    D(:,i) = h;
end
toc



% Define the objective values to reach with the filter Fi
P = (st_training(n+dt:n-1+k+dt));


% Tikohnov regularization
% lambda = sqrt(sum(P.^2))*0.001;
% if lambda>0
%     Dn = 0*D;
%     Dn(1:n,1:n)=lambda;
%     D = [D, Dn];
%     P = [P, zeros(m,k)];
% end

%% % Build AR filter
tic
[u,s,v]=svd(D, 'econ'); % U = u*s*v';
toc

snorm = s/max(s(:));
noise = 0
if noise==0
    threshold=1e-4;
else
    threshold = 1e-3;
end
% threshold = min(noise^2, 1e-2);
[~,nmax] = min(abs(diag(snorm)-threshold));

sc = zeros(size(s))+s;
% sc(nmax:end)=inf;

Dp = v(:,1:nmax)*pinv(sc(1:nmax,1:nmax))*u(:,1:nmax)';

% Tikohonov regularization
% Dp = (D'*D+lambda*eye(k))*D';

Fi = P*Dp;


%%
filename = [sprintf('AR%03.0f_delay_%.0f_noise_%.3f_f1_%.0f_f2_%.0f_', n, dt, noise, f1, f2), 'input_', inputType,'.png'];


%%
% Predictive reconstruction using past data
% % Test for 1 sample
% stnp1 = Fi*st1(1:n)';
% stnp1-st1(n-3:n+3)

% Building shifting data sets
% n == history size to predict step at iteration n+dt
hn  = zeros(n,Nsample);  % Prediction with training set statistics
h1  = zeros(n,Nsample);  % Prediction of noiseless original data set
h1n = zeros(n,Nsample); % Prediction of noisy original data set
h2  = zeros(n,Nsample);  % Prediction with another data set (nosieless)
h2n = zeros(n,Nsample); % Prediction with another data set (noisy)

for i=n:Nsample-dt
    % Prediction at step i+dt
    hn(:,i+dt)  = st_training(i-n+1:i);
    h1(:,i+dt)  = st1(i-n+1:i);
    h2(:,i+dt)  = st2(i-n+1:i);
    h1n(:,i+dt) = st1n(i-n+1:i);
    h2n(:,i+dt) = st2n(i-n+1:i);
end

% Reconstructed signal
st_   = Fi*hn;
st1_  = Fi*h1;
st2_  = Fi*h2;
st1n_ = Fi*h1n;
st2n_ = Fi*h2n;
stSI_ = circshift(st1n,2);


% figure, plot([st_; circshift(st1,-dt)]', '.-')

colLP   = "#D95319";
colSI   = "#EDB120";
colData = "#77AC30";
% colData = "green";

fig=figure();
fig.Position = [100,100,1500,1500];

subplot(421)
loglog(diag(snorm), 'k', 'LineWidth', 2)
hold on
loglog([nmax,nmax],[min(diag(snorm)), 1.2], 'k--')
hold off
grid()
xlim([1, 2*n])
xlabel('EV mode')
ylabel('EV amplitude')
title(sprintf('Filter order: %.0f, Nmodes: %.0f', n, nmax))

subplot(422)
stem(flip(Fi),'filled')
xlabel('Delay [frame]')
title('AR filter')

subplot(4,2,[3,4])
hold on
plot(st1, 'color', "#4DBEEE", 'LineWidth', 3, 'DisplayName', 'Disturbance')
plot(st_training,'.-', 'color', colData, 'DisplayName', 'Noisy data', 'MarkerSize', 10)
if strcmp(inputType, 'vibrations')
    plot(stSI_, '.-', 'color', colSI, 'LineWidth',1, 'DisplayName', 'SI')
end
plot(st1n_, '.-', 'color',  colLP, 'DisplayName', 'LP reconstruction')
grid()
% xscale('log')
xlim([10*n,Nsample])
ylabel('Signal [a.u.]')
legend('Location','best')
xlabel('Iteration')

resData = (st1 - st_training)/std(st1)*100;
errSI  = (st1 - stSI_)/std(st1)*100;
resLP1 = (st1 - st1n_)/std(st1)*100;

resData = (st1 - st_training)*1e9;
errSI  = (st1 - stSI_)*1e9;
resLP1 = (st1 - st1n_)*1e9;

subplot(4,2,5)
hold on
plot(resData, 'color', colData, 'DisplayName',sprintf('Data: %.2f nm RMS', std(resData)))
plot(errSI, 'color', colSI, 'DisplayName',sprintf('SI: %.2f nm RMS', std(errSI)))
plot(resLP1, 'color', colLP, 'DisplayName',sprintf('LP: %.2f nm RMS', std(resLP1((100+2*n):(end-2*n-100)))), 'LineWidth',2)
grid()
legend('Location', 'best')
ylabel('Error [nm RMS]')
xlim([10*n,10*n+max(n,1000)])
% ylim([-10,10])
% xscale('log')
xlabel('Iteration')
title(sprintf('Residual - f = %.0f Hz (training set)', f1))

resSI2 = (st2 - circshift(st2n,dt))/std(st2)*100;
resLP2 = (st2 - st2n_)/std(st2)*100; % Another signal reconstructed with original filter

subplot(4,2,6)
hold on
% plot(resData2, 'color', colData, 'DisplayName',sprintf('Data: %.2f RMS', std(resData)))
plot(resSI2, 'color', colSI, 'DisplayName',sprintf('SI: %.2f RMS', std(resSI2)))
plot(resLP2, 'color', colLP, 'DisplayName',sprintf('LP: %.2f RMS', std(resLP2(2*n:end))), 'LineWidth',2)
grid()
legend('Location', 'best')
ylabel('Error [%]')
xlim([2*n,2*n+1e3])
% ylim([-10,10])
% xscale('log')
xlabel('Iteration')
title(sprintf('Residual - f = %.0f Hz (new set)', f2))


subplot(427)
[xc,lc]=xcorr(st1, circshift(st1,dt),'biased');
plot(lc, xc/max(xc),'.-')
hold on
[xc,lc]=xcorr(st1, st1n_,'biased');
plot(lc, xc/max(xc),'.-')
[xc,lc]=xcorr(st2, st2n_,'biased');
plot(lc, xc/max(xc),'.-')
[xc,lc]=xcorr(st1, st1,'biased');
plot(lc, xc/max(xc),'-k')
hold off
legend('SI', 'LP1','LP2','Input', 'Location', 'best')
xlim([-20,20])
% ylim([0,1])
xlabel('Delay [frame]')
ylabel('Auto-correlation')
grid()


Nfft = 4000;
ttt = (1:Nsample)/fmax;

subplot(428)
pxx=periodogram(errSI,[],Nfft,fmax);
[psdSI, freq, ~]=genUtilities.PSDavg(ttt, errSI, Nfft);
loglog(freq(2:end), psdSI(2:end))
hold on
[psdLP1, freq, ~]=genUtilities.PSDavg(ttt, resLP1, Nfft);
plot(freq(2:end), psdLP1(2:end))
[psdLP2, freq, ~]=genUtilities.PSDavg(ttt, resLP2, Nfft);
loglog(freq(2:end), psdLP2(2:end))
grid()
hold off
% xscale('log')
% yscale('log')
xlabel('Frequency [Hz]')
ylabel('Residual spectrum [a.u.]')
ylim([min([prctile(psdSI,0.05),prctile(psdLP1,0.05)]), max(max(psdLP1), max(psdSI))])
xlim([freq(2), freq(end)])
legend('SI', 'LP (training set)','LP (new set)', 'Location', 'best')

% saveas(fig, fullfile(dfolder, filename))


%% Close loop test
gSI = 0.4;
gLP = 1;

errSI = zeros(1,Nsample); % Error signal
corSI = zeros(1,Nsample); % Integrator correction device
bufSI = zeros(1,dt+1);

errLP = zeros(1,Nsample); % Error signal
corLP = zeros(1,Nsample); % Integrator correction device
histLP = zeros(1,n);      % Buffer of previous measurements

distNoise = (st1-st1n)*1;

for i=1:Nsample-1
    % Simple integrator
    errSI(i) = st1(i) - corSI(i) + distNoise(i); % Sensor measurement
    bufSI(1) = errSI(i);           % Fill in buffer to apply the measurment with some delay
    corSI(i+1) = corSI(i) + gSI*bufSI(end); % Apply delayed correction
    bufSI = circshift(bufSI,-1);    % Circulate the buffer 1 step

    % Linear predictor
    if i<=n    % Start with SI solution
        histLP(1) = st1n(i); % Sensor measurement
        corLP(i+1) = st1n(i); % Sensor measurement
        histLP = circshift(histLP,1);    % Circulate the buffer 1 step
    else       % Run with LP
        errLP(i)   = st1(i) - corLP(i) + distNoise(i); % Sensor measurement
        histLP(1)  = corLP(i) + gLP*errLP(i); % Sensor measurement
        corLP(i+1) = Fi*histLP'; % Apply delayed correction
        histLP = circshift(histLP,-1);    % Circulate the buffer 1 step
    end

end

figure, plot([st1; st1; corSI; corLP; corSI-st1; corLP-st1]')
legend('Dist+noise', 'Dist.', 'SI', 'LP', ...
    sprintf('Res SI: %.2f nm RMS', std(corSI(n:end)-st1(n:end))*1e9), ...
    sprintf('Res LP: %.2f nm RMS', std(corLP(n:end)-st1(n:end))*1e9))
