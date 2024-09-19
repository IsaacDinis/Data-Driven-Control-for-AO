import numpy as np
import astropy.io.fits as pfits
class eof:
    def __init__(self,order,S2M, M2V,l,n_modes = 1, sys_delay = 2):
        self.n_modes = n_modes
        self.order = order
        self.history_size = int(self.n_modes*order)
        self.l = l
        self.S2M = S2M[:self.n_modes,:]
        self.F = np.zeros((self.n_modes,self.history_size))    
        self.D = np.zeros((self.history_size,l))
        self.P = np.zeros((self.n_modes,l))
        self.h = np.zeros(self.history_size)
        self.train_count = 0
        self.is_trained = 0
        self.M2V = M2V[:,:self.n_modes]
        self.V2M = np.linalg.pinv(M2V)
        self.sys_delay = sys_delay
        self.command_buf = np.zeros((self.n_modes,sys_delay))
        self.command = 0

    # def train(self,slopes,command):

    #     # modes = -self.S2M@slopes + self.command_buf[:,0]
    #     modes = self.S2M@slopes - self.V2M@command

    #     self.command_buf = np.roll(self.command_buf,1,axis=1)
    #     self.command_buf[:,0] =  -self.V2M@command

    #     if self.train_count < self.sys_delay:
    #         self.train_count += 1

    #     elif self.train_count < self.l+self.sys_delay+self.order:
    #         self.h = np.roll(self.h,self.n_modes)
    #         self.D = np.roll(self.D,1,axis=1)
    #         self.h[:self.n_modes] = modes
    #         self.D[:,0] = self.h
    #         self.P = np.roll(self.P,1,axis = 1)
    #         self.P[:,0] = modes
    #         self.train_count += 1

    #     elif self.train_count < self.l + 2*self.sys_delay +self.order:
    #         self.h = np.roll(self.D[:,0],self.n_modes)
    #         self.h[:self.n_modes] = modes
    #         self.P = np.roll(self.P,1,axis = 1)
    #         self.P[:,0] = modes
    #         self.train_count += 1

    #     else:
    #         self.D /= np.linalg.norm(self.D,axis = 1).reshape(-1,1)
    #         self.P /= np.linalg.norm(self.P,axis = 1).reshape(-1,1)
    #         D_inv = np.linalg.pinv(self.D.T, 1e-3)
    #         self.F = (D_inv @ self.P.T).T
    #         self.is_trained = 1
    #         print("EOF trained")
    #         self.save()
    #         # self.F /= np.max(self.F)
    #         self.h *= 0
    #         self.command_buf *= 0

    def train(self,slopes,command):

        # modes = -self.S2M@slopes + self.command_buf[:,0]

        self.h[:self.n_modes] = self.S2M@slopes + self.command
        self.command = -self.V2M@command

        if self.train_count < self.sys_delay:
            self.train_count += 1

        elif self.train_count < self.l+self.sys_delay+self.order:
            
            self.D = np.roll(self.D,1,axis=1)
            self.D[:,0] = self.h
            self.P = np.roll(self.P,1,axis = 1)
            self.P[:,0] = self.h[:self.n_modes]
            self.train_count += 1

        elif self.train_count < self.l + 2*self.sys_delay +self.order:
            self.P = np.roll(self.P,1,axis = 1)
            self.P[:,0] = self.h[:self.n_modes]
            self.train_count += 1
        else:
            self.D /= np.linalg.norm(self.D,axis = 1).reshape(-1,1)
            self.P /= np.linalg.norm(self.P,axis = 1).reshape(-1,1)
            D_inv = np.linalg.pinv(self.D.T, 1e-3)
            self.F = (D_inv @ self.P.T).T
            self.is_trained = 1
            print("EOF trained")
            self.save()

        self.h = np.roll(self.h,self.n_modes) 

    def update_command(self, slopes):

            # self.h = np.roll(self.h,self.n_modes)
            # modes = self.S2M@slopes + self.command_buf[:,0]
            # self.h[:self.n_modes] = modes
            # # self.h[:self.n_modes] = self.command_buf[:,0]
            # # self.h[self.n_modes:2*self.n_modes] = modes
            # self.command_buf = np.roll(self.command_buf,1,axis=1)
            # command = self.F@self.h
            # self.command_buf[:,0] = command

            self.h[:self.n_modes] = self.S2M@slopes + self.command
            self.command = self.F@self.h
            self.h = np.roll(self.h,self.n_modes)

            return -self.M2V@self.command

    def save(self):
        pfits.writeto("F.fits", self.F, overwrite = True)

        pfits.writeto("D.fits", self.D, overwrite = True)

        pfits.writeto("P.fits", self.P, overwrite = True)
