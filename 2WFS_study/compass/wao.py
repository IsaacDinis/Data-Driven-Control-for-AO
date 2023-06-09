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





command = np.array([0.165,0])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
a = wao.supervisor.target.get_tar_phase(0,pupil=True)
# print(np.std(a))
print(np.max(a))
print(np.min(a))
