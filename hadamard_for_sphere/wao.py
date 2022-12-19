ipython -i ~/compass/shesha/shesha/widgets/widget_ao.py compass_param.py 


n_act_square = wao.config.p_dms[0].get_nact()
v = np.zeros(n_act_square**2)

wao.supervisor.rtc.set_perturbation_voltage(0, "", v)
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)
dummy = S2A @ slopes

supervisor.rtc.set_perturbation_voltage(0, "", hadamard_matrix[:n_act_square,1]*ampli)

inf_mat, _ = wao.supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]

wao.supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,0])

x_pos = supervisor.config.p_dms[0].get_xpos()