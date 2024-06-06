x = sdpvar(4,1);
y = sdpvar(1,1);
c = [-1.42,-0.15,-0.26,2.23]';
a = [-2.43,0.11,0.37,1.35]';


CON = cone([y;a.*x+c]);
CON = [CON,y>=0];
OBJ = y;
diagonstics = optimize(CON,OBJ);

norm(a.*double(x)+c)