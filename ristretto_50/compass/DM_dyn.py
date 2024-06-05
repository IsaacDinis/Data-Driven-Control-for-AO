import numpy as np

class DM_dyn:
    def __init__(self, size):
        self.order = 2
        
        self.a = np.array([1, -0.00373488546341600,    3.48734235620900e-06])
        self.b = np.array([0.684151982031302, 0.310839005201411, 0.00127761464622797])
        self.state_mat = np.zeros((self.order+1,2,size))

 
  


    def update_command(self,voltages):

        self.state_mat[1:,:,:] = self.state_mat[0:-1,:,:]
        self.state_mat[0,0,:] = voltages
        self.state_mat[0,1,:] = 0
        new_voltages = np.dot(self.b,self.state_mat[:,0,:]) - np.dot(self.a,self.state_mat[:,1,:])
        


        self.state_mat[0,1,:] = new_voltages

        return new_voltages

