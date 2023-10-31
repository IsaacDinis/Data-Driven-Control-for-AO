act_dead = fitsread('target_phase_croped_35.fits');
pos_act = fitsread('pos_act.fits');

figure()
surf(act_dead)
hold on 

stem3(pos_act(:,1)',pos_act(:,2),3.5*ones(4,4),'filled','red')