%%
% path_to_fusion = "/home/isaac/mosek/9.3/tools/platform/linux64x86/bin/mosek.jar";
path_to_fusion = "/home/isaac/mosek/10.1/tools/platform/linux64x86/bin/mosek.jar";
fusion_chk = contains(javaclasspath('-dynamic'), "mosek", 'IgnoreCase', true);
if strlength(path_to_fusion) && ~sum(fusion_chk)
    javaaddpath(path_to_fusion);
end 
%%
fs = 1200;
TT_vibr = load('data/TT_pol_SPHERE.mat');
tilt_s1 = TT_vibr.TT_psol_s1(:,1);
% tilt_s1 = wdenoise(double(tilt_s1));
tilt_s1 = detrend(tilt_s1);
% tilt_s1 = tilt_s1(1:10000);
tilt_s1 = tilt_s1-mean(tilt_s1);
mas2um = 4.84e-9*8e6;
tilt_s1 = tilt_s1*mas2um;

size_fft = 500;

[tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
f = f(2:end);
tilt_s1_psd = tilt_s1_psd(2:end);

x = log10(f);
y = log10(tilt_s1_psd);
[fitresult, gof] = create_fit_pol4(x, y');
y_fit = feval(fitresult,x);

figure()
plot(x,y)
hold on
plot(x,y_fit)
legend('raw SPHERE data','amtosphere fit')
xlabel('freq (Hz)')
ylabel('tilt mag (dB)')
title('PSD')
make_it_nicer()
%%
atm_fit = 10.^y_fit;
atm_fit(1:8) = tilt_s1_psd(1:8);
psd = atm_fit;


figure()
semilogx(f,20*log10(tilt_s1_psd))
hold on
semilogx(f, 20*log10(atm_fit))
legend('raw SPHERE data','amtosphere fit')
xlabel('freq (Hz)')
ylabel('tilt mag (dB)')
title('PSD')
make_it_nicer()

figure()
semilogx(f,20*log10(tilt_s1_psd-atm_fit));

%%
G = tf([1],[1,0],1/fs); % 1 samples delay
bandwidth = 100;
order = 50;
max_control_gain = 2;

%% Load data
max_order =  max(order);
n_modes = 1;
Kdd_matrix = zeros(max_order+1,2,n_modes);

g = 0.2;
K0 = tf([g,0],[1,-1],1/fs);
S_int = feedback(1,G*K0);

tic
w = f*2*pi;
W1 = frd(psd,w,1/fs);
% W1 = frd(psd,w,1/fs)/abs(frd(S_int,w,'rad/s'));
val = freqresp(W1, bandwidth*2*pi);
% val = interp1(f,psd,bandwidth);
% val = 1;
W1 = W1/val;




% W13 = W1*W12;
W3 = tf(max_control_gain,1,1/fs);

% W32 = makeweight(0.001,200,10,1/fs);
% figure()
% bodemag(W1,W12,W13,W32)
% legend()

Kdd =  ao_dd_controller(fs,w,order,W1,W3,[],'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation


K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,G*Kdd);
sys_int = feedback(1,G*K0);

t = 0:1/fs:size(tilt_s1,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,tilt_s1,t);
[yint,~] = lsim(sys_int,tilt_s1,t);

rms_dd = rms(ydd(500:end));
rms_int = rms(yint(500:end));

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms_dd,rms_int);
gain = mean(100-rms_dd*100./rms_int);

[psd_dd, f] = compute_psd_fft(ydd,size_fft,fs);
[psd_int, f] = compute_psd_fft(yint,size_fft,fs);

figure()
semilogx(f,10*log10(psd_dd))

title('Tilt vibration PSD after filtering')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

grid()
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/residual_npnoise.pdf -transparent

%% Save controller coefficients
Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
% save('../Kdd_ProxCen','Kdd_matrix');
%%
S_dd = feedback(1,G*Kdd);
dummy = S_dd*W1*val;
figure()
bodemag(dummy,S_dd,W1)

S_int = feedback(1,G*K0);
dummy = S_int*W1*val;
figure()
bodemag(dummy,S_int,W1)

fitswrite(ydd(8:end),'tilt_vibration.fits')
f = f(2:end);
psd_dd = psd_dd(2:end);
%%


figure()
semilogx(f,20*log10(tilt_s1_psd))
hold on
semilogx(f,20*log10(psd_dd))
title('Tilt vibration identified PSD')
xlabel('freq (Hz)')
ylabel('magnitude (dB)')
legend('with atmosphere','atmosphere filtered out')
grid()
make_it_nicer()

set(gcf,'Position',[100 100 800 500])
set(gcf,'PaperType','A4')

