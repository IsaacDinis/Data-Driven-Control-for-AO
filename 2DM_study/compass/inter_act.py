pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
pos_LODM = pos_LODM.copy()
pos_LODM -= np.min(pos_LODM,axis = 0)
pos_LODM /= 32
iop = np.arange(0,99)
pos_LODM = pos_LODM.astype(int)
zbra[pos_LODM[:,0],pos_LODM[:,1]]= iop

pupil = np.zeros((11,11))
pupil[pos_LODM[:,0],pos_LODM[:,1]] = 1
pupil = pupil.astype(int)

pupil_roll = np.roll(pupil,1,axis = 0)
pupil_roll[0,:] = 0
zbra_roll = np.roll(zbra,1,axis = 0)

inter = np.zeros((11,11))
inter_flat = zbra[pupil*pupil_roll == 1]-zbra_roll[pupil*pupil_roll == 1]

inter[pos_LODM[:,0],pos_LODM[:,1]] = inter_flat

plt.imshow(inter)

pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T






pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
pos_LODM -= np.min(pos_LODM,axis = 0)
step = pos_LODM[1,0]-pos_LODM[0,0]
pos_LODM /= step
pos_LODM = pos_LODM.astype(int)

iop = np.arange(0,100)

zbra[pos_LODM[:,0],pos_LODM[:,1]]= iop

pupil = np.zeros((11,11))
pupil[pos_LODM[:,0],pos_LODM[:,1]] = 1
pupil = pupil.astype(int)

pupil_roll = np.roll(pupil,1,axis = 0)
pupil_roll[0,:] = 0
zbra_roll = np.roll(zbra,1,axis = 0)
inter_x = zbra[pupil*pupil_roll == 1]-zbra_roll[pupil*pupil_roll == 1]
print(inter_x)

pupil_roll = np.roll(pupil,1,axis = 1)
pupil_roll[:,0] = 0
zbra_roll = np.roll(zbra,1,axis = 1)
inter_y = zbra[pupil*pupil_roll == 1]-zbra_roll[pupil*pupil_roll == 1]
print(inter_y)




inter = np.array([inter_x,inter_y])