from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 

pupil_diam = wao.supervisor.config.p_geom.get_pupdiam()
pupil_grid = make_pupil_grid(pupil_diam)
zernike_basis = make_zernike_basis(3, 1, pupil_grid)
tilt = zernike_basis[2].shaped
pupil_valid = zernike_basis[0].shaped

command = wao.supervisor.rtc.get_command(0)
slopes = wao.supervisor.rtc.get_slopes(0)

command = np.array([0.001,0])
wao.supervisor.rtc.set_command(0,command)

a = wao.supervisor.target.get_tar_phase(0,pupil=True)
# print(np.std(a))
print(np.max(a))
print(np.min(a))

slopes = wao.supervisor.rtc.get_slopes(0)

M2S = np.zeros((slopes.shape[0], 2))
ampli = 0.001
command = np.array([ampli,0])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)/ampli
M2S[:,0] = slopes.copy()

command = np.array([0,ampli])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)/ampli
M2S[:,1] = slopes.copy()

S2M = np.linalg.pinv(M2S)



command = np.array([1,0])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)
modes= np.dot(S2M,slopes)
print(modes[0])
print(modes[1])


lambd = 1.65
D = 8.




command = np.array([1.65/5/D,0])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)*2*np.pi/lambd
# print(np.std(a))
print(np.max(target_phase))
print(np.min(target_phase))
print(np.std(target_phase,where = pupil_valid.astype(bool)))

pupil_valid2 = np.pad(pupil_valid,8)
a = wao.supervisor.wfs.get_wfs_phase(0)*pupil_valid2
# print(np.std(a))
print(np.max(a))
print(np.min(a))