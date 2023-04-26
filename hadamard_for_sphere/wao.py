ipython -i ~/compass/shesha/shesha/widgets/widget_ao.py saxo.py 

ipython -i ~/shesha/shesha/widgets/widget_ao.py saxo.py 


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

inf_mat, _ = wao.supervisor.basis.compute_modes_to_volts_basis("KL2V") #
   ...:  or "KL2V" [nvolts ,nmodes]
Computing KL2V basis...
Computing IF sparse...
[/root/miniconda3/conda-bld/compass_1654803282844/work/libcarma/src.cpp/carma_sparse_obj.cpp@87]: Warning : empty CarmaObj cannot be sparsed


