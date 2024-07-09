
gam = sdpvar;
W1 = 19;
Y_n = sdpvar(2,1);
X_n = sdpvar(2,1);
Y_n(1) = 1;
XY_n = [X_n;Y_n];

ZFy = [1.0000 + 0.0050i,1.0000 + 0.0000i];
ZFx = [1.0000 + 0.0050i,1.0000 + 0.0000i];
F_a = W1*ZFy*Y_n;
Pc = 1.2000 + 0.0050i;
% Pc = 1000;
G = 10000.0000 - 100.0050i;

PHI = 2*real([G.*ZFx, ZFy].*conj(Pc))*XY_n-abs(Pc).^2;
% gam = 0.5*gamma_2;
% rcone = utils.rcone_serialized(abs(Pc)^2,gam,F_a);
rcone = utils.rcone_serialized(PHI,gam, F_a);
CON = [rcone,gam>=0, PHI>=0];

% CON = [CON, integ.'*gamma_2 <= obj_2 ];
OBJ = gam;
diagonstics = optimize(CON,OBJ);

