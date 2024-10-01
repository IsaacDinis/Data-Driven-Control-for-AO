import numpy as np
from dd4compass import K_dd

class K_dd_ris(K_dd):
    def __init__(self,order,S2M,M2V,boostrap_n_iter,fs,offload_mat = np.empty(0),stroke = np.inf,offload_ratio = 0):
        delay = 2
        super().__init__(order, delay, S2M, M2V, boostrap_n_iter, fs) 
        self.S2M = S2M
        self.M2V = M2V
        self.res = np.zeros(self.S2M.shape[0])
        self.V2M = np.linalg.pinv(M2V)
        self.stroke = stroke
        self.offload_mat = offload_mat
        self.offload_count = 0
        if (offload_mat.any()):   
            self.offload_mat_T = np.linalg.pinv(offload_mat)
        self.u_offload_HODM = 0
        self.u_offload_LODM = 0
        self.offload_ratio = offload_ratio 

    def update_command(self,slopes):
        if self.status == "initialized":
            u = self.bootstrap(slopes)
        elif self.status == "trained":
            u = self.compute_voltage(slopes)

        if (self.offload_mat.any()):
            if self.offload_count % self.offload_ratio == 0:
                self.u_offload_LODM = self.offload_mat@u
                self.u_offload_HODM = self.offload_mat_T@self.u_offload_LODM
            u -= self.u_offload_HODM
            self.offload_count += 1

        act_sat = np.argwhere(np.abs(u)>self.stroke)
        u[act_sat] = np.sign(u[act_sat])*self.stroke

        if (self.offload_mat.any()):
            u_tot = u+self.u_offload_HODM
        else:
            u_tot = u
        if (self.offload_mat.any()):
            return u,self.u_offload_LODM
        else:
            return u
