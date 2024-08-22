fs = 1000;
t = 0:1/fs:1-1/fs;
x = sin(2*pi*50*t);
% x = fitsread('data2/disturbance.fits');
tilt_res = fitsread('tilt_res.fits');
tilt_applied = fitsread('tilt_applied.fits');
x = tilt_applied(1:end-2)+tilt_res(1+2:end);

len = length(x);
order = 10;
arcoef = aryule(x,order);
sys = ar(x,order,'ls');
arcoef = -arcoef(2:end);
arcoef = -sys.A(2:end);
x_pred = zeros(1,len-order);
x_pred2 = zeros(1,len-order);
x_pred3 = zeros(1,len-order);

for i = order:len
    x_pred(i) = dot(arcoef,x(i:-1:i-order+1));
    x_pred2(i) = dot(arcoef(2:end),x(i:-1:i-order+2))+x_pred(i)*arcoef(1);
end

for i = order:len
    x_pred3(i) = dot(arcoef,x_pred(i:-1:i-order+1));
end


% figure()
% plot(x(53:70))
% hold on;
% plot(x_pred2(51:70))
% plot(x(51:70))
% legend('true value','prediction','measurement')

figure()
stem(arcoef)
% 
figure()
plot(x(53:70))
hold on;
plot(x_pred2(51:70))
plot(x(51:70))
legend('true value','prediction','measurement')


% figure()
% plot(x(53:end))
% hold on;
% plot(x_pred2(51:end))
% plot(x(51:end))
% legend('true value','prediction','measurement')