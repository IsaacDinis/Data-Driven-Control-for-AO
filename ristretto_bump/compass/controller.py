import numpy as np

class K:
    def __init__(self,order,a,b,S2M,M2V):
        self.order = order
        
        self.a = a
        self.b = b
        self.S2M = S2M
        self.M2V = M2V
        self.res = np.zeros(self.S2M.shape[0])
        self.u = np.zeros(self.M2V.shape[0])
        self.state_mat = np.zeros((self.order+1,self.order+1,M2V.shape[1]))

    def compute_modal_res(self,slopes):
        self.res = self.S2M @ slopes

    def update_command(self,slopes):
        self.compute_modal_res(slopes)
        self.state_mat[1:,:,:] = self.state_mat[0:-1,:,:]
        self.state_mat[0,0,:] = self.res
        self.state_mat[0,1,:] = 0
        modal_u = np.dot(self.b,self.state_mat[:,0,:]) - np.dot(self.a,self.state_mat[:,1,:])
        self.state_mat[0,1,:] = modal_u
        self.u = -self.M2V @ modal_u
        return self.u
