import numpy as np
from matplotlib import pyplot as plt
import astropy.io.fits as pfits

class DM_stroke_plot:
    def __init__(self,title, refresh_rate, n_act, n_iter, act_pos, cross_act):

        # self.order = order

        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.stroke = np.zeros((n_iter,n_act))
        self.inter_stroke = np.zeros(n_iter)
        self.count = 1
        self.line_max_stroke, = self.ax.plot(self.stroke[:,0])
        self.line_std_stroke, = self.ax.plot(self.stroke[:,0])
        self.line_inter_stroke, = self.ax.plot(self.inter_stroke)
        self.ax.set_ylabel("[um]")
        self.ax.set_xlabel("iter")
        self.line_max_stroke.set_label('max stroke')
        self.line_inter_stroke.set_label('max inter actuator stroke')
        self.line_std_stroke.set_label('std stroke')
        self.n_iter = n_iter

        act_pos -= np.min(act_pos,axis = 0)
        step = np.max([act_pos[1,0]-act_pos[0,0],act_pos[0,1]-act_pos[0,0]])
        act_pos /= step
        act_pos = act_pos.astype(int)
        self.act_pos = act_pos

        self.command_2D  = np.zeros((cross_act,cross_act))

        pupil = np.zeros((cross_act,cross_act))
        pupil[self.act_pos[:,0],self.act_pos[:,1]] = 1
        pupil = pupil.astype(int)

        pupil_roll = np.roll(pupil,1,axis = 0)
        pupil_roll[0,:] = 0
        self.pupil_axis_0 = pupil*pupil_roll == 1

        pupil_roll = np.roll(pupil,1,axis = 1)
        pupil_roll[:,0] = 0
        self.pupil_axis_1 = pupil*pupil_roll == 1

    def plot(self, DM_command,iter_n):
        
        self.stroke[iter_n,:] = np.abs(DM_command)

        self.command_2D[self.act_pos[:,0],self.act_pos[:,1]] = DM_command
        command_2D_roll_axis_0 = np.roll(self.command_2D,1,axis = 0)
        command_2D_roll_axis_1 = np.roll(self.command_2D,1,axis = 1)


        inter_stroke_x = self.command_2D[self.pupil_axis_0]-command_2D_roll_axis_0[self.pupil_axis_0]
        inter_stroke_y = self.command_2D[self.pupil_axis_1]-command_2D_roll_axis_1[self.pupil_axis_1]
        inter_stroke = np.array([inter_stroke_x,inter_stroke_y])
        self.inter_stroke[iter_n] = np.max(np.abs(inter_stroke))

        if self.count == 0:
            max_stroke = np.max(self.stroke[:iter_n,:],axis = 1)
            std_stroke = np.std(self.stroke[:iter_n,:],axis = 1)


            self.line_max_stroke.set_ydata(max_stroke)
            self.line_inter_stroke.set_ydata(self.inter_stroke[:iter_n])
            self.line_std_stroke.set_ydata(std_stroke)

            self.ax.set_ylim(0,np.max(np.array([max_stroke,self.inter_stroke[:iter_n],std_stroke])))
            self.ax.set_xlim(0,iter_n)

            self.line_max_stroke.set_xdata(np.arange(iter_n))
            self.line_std_stroke.set_xdata(np.arange(iter_n))
            self.line_inter_stroke.set_xdata(np.arange(iter_n))

            self.ax.set_ylabel("[um]")
            self.ax.set_xlabel("iter")
            self.ax.legend()
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.stroke *= 0
        self.inter_stroke *= 0

    def save(self, path_name):
        pfits.writeto(path_name, self.stroke, overwrite = True)

    def load(self,path_name):
        self.stroke = pfits.getdata(path_name)
        
    def save_plot(self, path_name):
        self.fig.savefig(path_name) 