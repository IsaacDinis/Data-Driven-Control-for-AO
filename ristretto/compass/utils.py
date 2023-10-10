
def save_perf(path,exposure_time,strehl,phase_rms):
	with open(path+'perf.txt', 'w') as f:
	    f.write('exp time = {:.1f}s, strehl = {:.3f}, phase rms = {:.3f}um'.format(exposure_time,strehl,phase_rms))

