x = [-10:.1:10];
y = normpdf(x,0,1);
y = y/max((y));
y1 = normpdf(x,-2,1);
y1 = y1/max((y1));
y2 = normpdf(x,2,1);
y2 = y2/max((y2));
y3 = normpdf(x,0,2.355);
y3 = y3/max((y3));
bump = normpdf(x,1,1);
bump = -bump/max((bump));

figure()
plot(x,y);
hold on
plot(x,y2)
plot(x,y1)
plot(x,y3)

figure()
plot(x,0.91*y+y1+y2);

%%
figure()
phase = y3+0.2*y1+0.2*y2;
plot(x,phase)

%%
y4 = normpdf(x,1,2.355);
y4 = y4/max((y4));
figure()
phase = y4-y-y2;
plot(x,phase)

%%
figure()
plot(x,0.91*y+y1+y2+bump);
%%
figure()
plot(x,y4+bump);

%%
figure()
plot(x,bump+y+y2);