pupil = supervisor.get_s_pupil()

phase2nm = 1e3/np.sqrt(np.sum(pupil))
ampli=1
mode = 1
command = np.concatenate((np.zeros(n_actus_DM0), M2V_DM1[:,mode]*ampli,np.zeros(n_actus_bump)), axis=0)
# command = np.concatenate((M2V_DM0[:,mode]*ampli,np.zeros(n_actus_DM1+n_actus_bump)), axis=0)


supervisor.rtc.set_command(0, command) 
supervisor.next()
supervisor.next()
supervisor.next()
supervisor.next()
phase = supervisor.target.get_tar_phase(0)
phase_rms = np.std(phase,where = pupil.astype(bool))
print(phase_rms)