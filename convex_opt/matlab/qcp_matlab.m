x = sdpvar();
y = sdpvar();
z = sdpvar();


% gam = 0.5*gamma_2;
% rcone = utils.rcone_serialized(abs(Pc)^2,gam,F_a);
c1 = x + y + z == 1;
c2 = x + y - z <= 0;
% c3 = rcone(sqrt(2)*x,y,z);
c3 = utils.rcone_serialized(y,z,sqrt(2)*x);

CON = [c1,c2,c3];

% CON = [CON, integ.'*gamma_2 <= obj_2 ];
OBJ = 1/x;
diagonstics = optimize(CON,OBJ);