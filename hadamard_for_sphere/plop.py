n_act = wao.config.p_dms[0].get_xpos().shape[0] # get number of visible actuators

n_slopes = wao.supervisor.rtc.get_slopes(0).shape[0]

ampli = 0.1

hadamard_size = 2**n_act.bit_length()
hadamard_matrix = hadamard(hadamard_size)

C = np.zeros((n_slopes,hadamard_size))

# # ampli = 0.01
# slopes = wao.supervisor.rtc.get_slopes(0)
# imat = np.zeros((slopes.shape[0], nmodes))
# # imat = np.zeros((slopes.shape[0], nmodes+2))



#-----------------------------------------------
# compute the command matrix [nmodes , nslopes]
#-----------------------------------------------
for i in range(hadamard_size):
    wao.supervisor.rtc.set_perturbation_voltage(0, "", hadamard_matrix[:n_act,i]*ampli)
    wao.supervisor.next()
    wao.supervisor.next()
    slopes = wao.supervisor.rtc.get_slopes(0)/ampli
    C[:,i] = slopes

D = C @ np.linalg.inv(hadamard_matrix)
# D = D[:,:n_act_square]
D = D[:,:n_act]
S2A = np.linalg.pinv(D)

v *= 0
x = 67
v[x] = 0.1
wao.supervisor.rtc.set_perturbation_voltage(0, "", v)
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)
amplitude = S2A@slopes
amplitude[x-1:x+2]