import numpy as np

class K:
    def __init__(self,order,a,b,S2M,M2V, offload_mat = np.empty(0), stroke = np.inf):
        self.order = order
        
        self.a = a
        self.b = b
        self.S2M = S2M
        self.M2V = M2V
        self.res = np.zeros(self.S2M.shape[0])
        
        self.state_mat = np.zeros((self.order+1,self.order+1,M2V.shape[1]))
        self.V2M = np.linalg.pinv(M2V)
        self.stroke = stroke
        self.offload_mat = offload_mat
        self.count = 0
        if (offload_mat.any()):   
            self.offload_mat_T = np.linalg.pinv(offload_mat)

    def compute_modal_res(self,slopes):
        self.res = self.S2M @ slopes

    def update_command(self,slopes):
        self.compute_modal_res(slopes)
        self.state_mat[1:,:,:] = self.state_mat[0:-1,:,:]
        self.state_mat[0,0,:] = self.res
        self.state_mat[0,1,:] = 0
        modal_u = np.dot(self.b,self.state_mat[:,0,:]) - np.dot(self.a,self.state_mat[:,1,:])
        
        u = -self.M2V @ modal_u
        if (self.offload_mat.any()):
            u_offload_LODM = self.offload_mat@u
            u_offload_HODM = self.offload_mat_T@u_offload_LODM
            u -= u_offload_HODM

        act_sat = np.argwhere(np.abs(u)>self.stroke)
        u[act_sat] = np.sign(u[act_sat])*self.stroke

        if (self.offload_mat.any()):
            u_tot = u+u_offload_HODM
        else:
            u_tot = u

        self.state_mat[0,1,:] = -self.V2M@u_tot
        # self.state_mat[0,1,:] = modal_u
        if (self.offload_mat.any()):
            return u,u_offload_LODM
        else:
            return u
