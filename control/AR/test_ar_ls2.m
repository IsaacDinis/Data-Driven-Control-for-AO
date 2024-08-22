fs = 4000;


tilt_res = fitsread('tilt_res.fits');
tilt_applied = fitsread('tilt_applied.fits');
x = tilt_applied(1:end-2)+tilt_res(1+2:end);
% x = x(1:10000)';
x = x';
y = x(6:10005);
m = 5;
N = length(x);

T = [0:1:N-1]/fs;

PHI2 = [x(1:10000),x(2:10001),x(3:10002),x(4:10003),x(5:10004)];

arcoef = (PHI2'*PHI2)\ PHI2' * y;
arcoef = arcoef(end:-1:1);
% y_hat = PHI*theta;
% loss = sum((x(2:end)-y_hat).^2);
% fprintf(" Loss function of our prediction : %f \n",loss);
% 

x_pred = zeros(1,N-m);
x_pred2 = zeros(1,N-m);


for i = m:N
    x_pred(i) = dot(arcoef,x(i:-1:i-m+1));
    x_pred2(i) = dot(arcoef(2:end),x(i:-1:i-m+2))+x_pred(i)*arcoef(1);
end

% x_pred2 = PHI*arcoef;
figure()
plot(x(53:70))
hold on;
plot(x_pred2(51:70))
plot(x(51:70))
legend('true value','prediction','measurement')
