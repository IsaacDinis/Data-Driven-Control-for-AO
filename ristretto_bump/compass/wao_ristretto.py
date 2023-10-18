#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/ristretto_bump/compass/ristretto_param.py

from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
import astropy.io.fits as pfits

# print(wao.supervisor.config.p_dms[0].get_ntotact())
# print(wao.supervisor.config.p_dms[1].get_ntotact())
# print(wao.supervisor.rtc.get_slopes(0).shape)

M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V_DM0.npy')
S2M_DM0 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/S2M_DM0.npy')
S2M_DM1 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/S2M_DM1.npy')
command = wao.supervisor.rtc.get_command(0)
M2V = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V.npy')

n_act0 = wao.supervisor.config.p_dms[0].get_ntotact()
n_act1 = wao.supervisor.config.p_dms[1].get_ntotact()
n_act2 = wao.supervisor.config.p_dms[2].get_ntotact()

n_modes = 500

wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,n_modes]), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

wao.supervisor.rtc.set_command(0,M2V[:,n_act0+800]*0.01)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)

modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[1])
print(modes_DM1[1])

wao.supervisor.rtc.set_command(0,M2V[:,0]*0)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,n_modes]*0), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

n_modes = 500
slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,:n_modes]@modes_DM1[:n_modes]), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()



print(modes_DM0[1])
print(modes_DM1[1])

norm = np.zeros(800)
for i in range(800):
    wao.supervisor.rtc.set_command(0,M2V[:,n_act0+i])
    wao.supervisor.next()
    wao.supervisor.next()
    wao.supervisor.next()
    target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
    norm[i] = np.std(target_phase)
norm /=np.max(norm)   
plt.plot(norm)
plt.show()



pos_LODM = np.array([wao.supervisor.config.p_dms[0].get_xpos(),wao.supervisor.config.p_dms[0].get_ypos()]).T
pos_HODM = np.array([wao.supervisor.config.p_dms[1].get_xpos(),wao.supervisor.config.p_dms[1].get_ypos()]).T
pos_HODM -= np.min(pos_HODM,axis = 0)

command *= 0
command[-1] = 1
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
plt.imshow(target_phase)


# 685 686 606 765
p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

middle = (np.max(xpos)-np.min(xpos))/2
# xpos0 = 240.5 # 960 961
# ypos0 = 184.5
xpos0 = 240.0 # 960 961
ypos0 = 185.8

plt.scatter(xpos, ypos, marker='.', color="red")
plt.scatter(xpos0, ypos0, marker='.', color="green")

influ = p_dm._influ
plt.imshow(influ[:,:,800])
influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)
i10 = xpos0-20.5
j10 = ypos0-20.5 
xcenter = p_geom.cent
ycenter = p_geom.cent

xpos0 += p_dm._n1
ypos0 += p_dm._n1
i10 += p_dm._n1
j10 += p_dm._n1
file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i10,j10,influ0,xpos0,ypos0,xcenter,ycenter,pixsize,diam)





p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
diam = p_tel.diam

xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

# place bump
xpos0 = 240.0 
ypos0 = 185.8

# check bump position
plt.scatter(xpos, ypos, marker='.', color="red")
plt.scatter(xpos0, ypos0, marker='.', color="green")


influ = p_dm._influ
influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)

i10 = xpos0-p_dm._xpos[0]-p_dm._n1 
j10 = ypos0-p_dm._xpos[0]-p_dm._n1

xcenter = p_geom.cent
ycenter = p_geom.cent

xpos0 += p_dm._n1
ypos0 += p_dm._n1
i10 += p_dm._n1
j10 += p_dm._n1

file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i10,j10,influ0,xpos0,ypos0,xcenter,ycenter,pixsize,diam)




p_dm._i1
wao.supervisor.config.p_dms[2].set_n1(331)