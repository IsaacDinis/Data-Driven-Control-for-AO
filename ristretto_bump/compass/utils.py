import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
import compute_psd
import astropy.io.fits as pfits

def save_perf(path,exposure_time,strehl,phase_rms):
	with open(path+'perf.txt', 'w') as f:
	    f.write('exp time = {:.1f}s, strehl = {:.3f}, phase rms = {:.3f}um'.format(exposure_time,strehl,phase_rms))

class zernike_plot:
    def __init__(self, title, refresh_rate, order, pup_diam, pupil, n_iter):
        pupil_grid = make_pupil_grid(pup_diam)
        zernike_basis = make_zernike_basis(order, 1, pupil_grid)
        transform_matrix = zernike_basis.transformation_matrix
        transform_matrix = np.linalg.pinv(transform_matrix)
        transform_matrix /= np.std(transform_matrix[1,:],where = transform_matrix[1,:]!=0)
        transform_matrix /= np.sum(pupil)
        # self.order = order
        self.transform_matrix = transform_matrix[1:,:]
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        
        # self.line, = self.ax.plot(self.zernike_res)
        self.ax.stem(np.zeros(order-1))
        self.ax.set_ylabel("[nm]")
        self.ax.set_xlabel("order")
        self.zernike_res = np.zeros((n_iter,order-1))
        self.n_iter = n_iter


    def plot(self, phase,iter_n):
        
        if np.sum([phase!=0]) !=0:
            zernike_res = self.transform_matrix@phase.flatten()*1e3
        else:
            zernike_res = np.zeros_like(self.zernike_res[0,:])
        self.zernike_res[iter_n,:] = zernike_res

        subtitle = 'tot = {:.3f}'.format((np.sqrt(np.sum(np.square(zernike_res)))))

        if self.count == 0:
            # self.line.set_ydata(zernike_rms)
            self.ax.cla()
            self.ax.stem(np.abs(zernike_res))
            self.ax.set_ylim(0,np.max(np.abs(zernike_res)))
            self.ax.set_title(subtitle)
            self.ax.set_ylabel("[nm]")
            self.ax.set_xlabel("order")
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.zernike_res *= 0


    def plot_psd(self,mode,fs,title,path_name = ''):
        n_average = 10
        window_size = int(np.floor(self.n_iter/n_average))
        psd, freq = compute_psd.compute_psd_fft(self.zernike_res, n_average, window_size,fs)
        fig, ax = plt.subplots(constrained_layout=True)
        # psd /= np.max(psd)
        psd = 10*np.log10(psd)
        fig.suptitle(title)
        ax.set_title('mode {:d}'.format(mode))
        ax.set_xscale('log')
        ax.set_ylabel("mag. [dB]")
        ax.set_xlabel("freq. [Hz]")
        ax.plot(freq,psd[:,mode])
        if (path_name != ''):
            fig.savefig(path_name) 

    def save(self, path_name):
        pfits.writeto(path_name, self.zernike_res, overwrite = True)

    def load(self,path_name):
        self.zernike_res = pfits.getdata(path_name)
        
    def save_std_plot(self, path_name):
        zernike_res_rms = np.std(self.zernike_res,axis = 0)
        subtitle = 'tot = {:.3f}'.format((np.sqrt(np.sum(np.square(zernike_res_rms)))))
        self.ax.cla()
        self.ax.stem(zernike_res_rms)
        self.ax.set_ylim(0,np.max(zernike_res_rms))
        self.ax.set_title(subtitle)
        self.ax.set_ylabel("[nm]")
        self.ax.set_xlabel("order")
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        self.fig.savefig(path_name) 


class phase_plot:
    def __init__(self, title, refresh_rate):
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.im = self.ax.imshow(np.zeros((8,8)))
        self.cbar = self.fig.colorbar(self.im, ax=self.ax)
        self.cbar.set_label(label="[um]", size=12)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        self.opd = 0
        self.opd_rms = 0

    def plot(self, phase, subtitle = "", iter_n = 0):
        # self.opd += np.sum(np.square(phase)/np.sum([phase!=0]))
        if np.sum([phase!=0]) !=0:
            self.opd += np.sum(np.abs(phase)/np.sum([phase!=0]))
        iter_n += 1
        # self.opd_rms = np.sqrt(self.opd/iter_n)*1e3
        self.opd_rms = (self.opd/iter_n)*1e3

        if self.count == 0:
            # im = self.ax.imshow(phase)
            self.im.set_data(phase)
            self.im.autoscale()

            self.ax.set_title(subtitle)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.opd = 0
        self.opd_rms = 0

class modal_plot:
    def __init__(self, title, refresh_rate, order, n_iter):

        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        # self.line, = self.ax.plot(self.modal_res)
        self.ax.stem(np.zeros(order))
        self.ax.set_ylabel("[au]")
        self.ax.set_xlabel("order")
        self.modal_res = np.zeros((n_iter,order))
        self.n_iter = n_iter

    def plot(self, modal_res,iter_n, bool_update = True):
        if bool_update:
            self.modal_res[iter_n,:] = modal_res

        subtitle = 'tot = {:.3f}'.format(np.sqrt(np.sum(np.square(modal_res))))

        if self.count == 0:
            # self.line.set_ydata(modal_rms)
            self.ax.cla()
            self.ax.stem(np.abs(modal_res),'-')
            self.ax.set_ylim(0,np.max(np.abs(modal_res)))
            self.ax.set_title(subtitle)
            self.ax.set_ylabel("[au]")
            self.ax.set_xlabel("order")
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.modal_res *= 0

    def plot_psd(self,mode,fs,title,path_name = ''):
        n_average = 10
        window_size = int(np.floor(self.n_iter/n_average))
        psd, freq = compute_psd.compute_psd_fft(self.modal_res, n_average, window_size,fs)
        # psd /= np.max(psd)
        psd = 10*np.log10(psd)
        fig, ax = plt.subplots(constrained_layout=True)
        fig.suptitle(title)
        ax.set_title('mode {:d}'.format(mode))
        ax.set_xscale('log')
        ax.set_ylabel("mag. [dB]")
        ax.set_xlabel("freq. [Hz]")
        ax.plot(freq,psd[:,mode])
        if (path_name != ''):
            fig.savefig(path_name) 

    def save(self, path_name):
        pfits.writeto(path_name, self.modal_res, overwrite = True)

    def load(self,path_name):
        self.modal_res = pfits.getdata(path_name)

    def save_std_plot(self, path_name):
        modal_res_rms = np.std(self.modal_res,axis = 0)
        self.plot(modal_res_rms,0,False)
        self.fig.savefig(path_name) 


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
        # self.line_std_stroke, = self.ax.plot(self.stroke[:,0])
        # self.line_inter_stroke, = self.ax.plot(self.inter_stroke)
        self.ax.set_ylabel("[um]")
        self.ax.set_xlabel("iter")
        self.line_max_stroke.set_label('max stroke')
        # self.line_inter_stroke.set_label('max inter actuator stroke')
        # self.line_std_stroke.set_label('std stroke')
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
        
        # self.stroke[iter_n,:] = np.abs(DM_command)

        self.stroke[iter_n,:] = np.max(np.abs([DM_command[941],DM_command[942],DM_command[980],DM_command[981]]))
        
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
            # self.line_inter_stroke.set_ydata(self.inter_stroke[:iter_n])
            # self.line_std_stroke.set_ydata(std_stroke)

            # self.ax.set_ylim(0,np.max(np.array([max_stroke,self.inter_stroke[:iter_n],std_stroke])))
            self.ax.set_ylim(0,np.max(np.array([max_stroke])))
            self.ax.set_xlim(0,iter_n)

            self.line_max_stroke.set_xdata(np.arange(iter_n))
            # self.line_std_stroke.set_xdata(np.arange(iter_n))
            # self.line_inter_stroke.set_xdata(np.arange(iter_n))

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

def write_dm_custom_fits(file_name, i1, j1, influ_cube, xpos, ypos, xcenter, ycenter,pixsize,pupm, *,pitchm=None):
    """Write a custom_dm fits file based on user provided data (see args)

    Args:
        file_name : (string) : name of the custom dm fits file

        i1 : (np.ndarray) : x coordinate of the first pixel for each of the 2D maps

        j1 : (np.ndarray) : y coordinate of the first pixel for each of the 2D maps

        influ_cube : (np.ndarray) : 2D maps of the influence functions

        xpos : (np.ndarray) : x coordinate of the physical location of the actuator

        ypos : (np.ndarray) : y coordinate of the physical location of the actuator

        xcenter : (float) : x coordinate of the centre of the pupil, expressed in pixels

        ycenter : (float) : y coordinate of the centre of the pupil, expressed in pixels

        pixsize : (float) : size of the pixels on the maps in meters

        pupm : (float) : diameter of pupil stop (meters)

    Kwargs:
        pitchm : (float) : size of the DM pitch in meters. Defaults to None.

    Returns:
        (HDUList) : custom_dm data
    """
    dm_fits_content="""The DM FITS file is compatible with COMPASS DM database.
    The primary header contains the keywords:
    * PIXSIZE : the size of the pixels on the maps in meters.

    * XCENTER, YCENTER are the coordinates of the centre of the pupil, expressed in pixels, in a reference frame conformable to (i,j) coords. The translation from pixels to meters can be done using:
        meters = (pixels - XCENTER) * PIXSIZE

    * PUPM is the diameter of pupil stop (meters).

    * Additionally the header provides the user with:
        PITCHM is the size of the DM pitch in meters (may be handy is some cases and useful when using this file with COMPASS software)

    This FITS file contains 3 extensions:
    * Extension 'I1_J1' are the coordinate (i,j) of the first pixel for each of the 2D maps (see Extension 2), so that they can be inserted in a larger map.

    * Extension 'INFLU' are the 2D maps of the influence functions.

    * Extension 'XPOS_YPOS' are the coordinates (xpos, ypos) of the  physical location of the actuator, in pixels. This data is provided for information only and does not directly participate to build the DM. The present coordinates are positions in M1 space, i.e. include the distorsion due to telescope optics.
    """

    fits_version=1.2
    primary_hdu = pfits.PrimaryHDU()
    primary_hdu.header['VERSION'] = (fits_version,'file format version')
    primary_hdu.header['XCENTER'] = (xcenter     ,'DM centre along X in pixels')
    primary_hdu.header['YCENTER'] = (ycenter     ,'DM centre along Y in pixels')
    primary_hdu.header['PIXSIZE'] = (pixsize     ,'pixel size (meters)')
    primary_hdu.header['PUPM']    = (pupm        ,'nominal pupil diameter (meters)')
    if(pitchm is not None):
        primary_hdu.header['PITCHM'] = (pitchm,'DM pitch (meters)')

    for line in dm_fits_content.splitlines():
        primary_hdu.header.add_comment(line)

    image_hdu = pfits.ImageHDU(np.c_[i1 , j1 ].T, name="I1_J1")
    image_hdu2 = pfits.ImageHDU(influ_cube, name="INFLU")
    image_hdu3 = pfits.ImageHDU(np.c_[xpos, ypos].T, name="XPOS_YPOS")

    dm_custom = pfits.HDUList([primary_hdu, image_hdu, image_hdu2, image_hdu3])

    dm_custom.writeto(file_name,overwrite=1)
    return dm_custom