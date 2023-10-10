#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/ristretto/compass/ristretto_param.py

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