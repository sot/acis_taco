#!/usr/bin/env python
import Tkinter as Tk
import sys
import cPickle as pickle

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

import Ska.quatutil
from Chandra.Time import DateTime


def destroy(e): sys.exit()

#def update_alpha(*args):
#    image.set_alpha(alpha_slider.value.get())
#    canvas.draw()

def get_index_lims():
    center = int(date_slider.value.get())
    wid = int(width_slider.value.get() * 6)
    i0 = center - wid
    i1 = center + wid
    if i0 < 0:
        i0 = 0
    if i1 >= len(imgs):
        i1 = len(imgs)-1
    return i0, i1

def update_image(*args):
    i0, i1 =get_index_lims()
    img = imgs[i0:i1+1].sum(0)
    img_rgba = matplotlib.cm.jet(img / 7.2)
    img_rgba[:, :, 3] = alpha_slider.value.get()
    image.set_data(img_rgba)
    canvas.draw()

def motion(event):
    if event.inaxes is None:
        return
    x = event.xdata
    y = event.ydata
    # print antisun.img2sky(x, y, sun_eci)
    
class SolarSystemObject(object):
    def __init__(self, name, times, xyzs, color, npoly=40):
        self.name = name
        self.times = times
        self.xyzs = xyzs
        self.color = color
        self.patches = []

        for idx in range(len(times)):
            xyz = ephem_xyzs[name][:, idx]
            dist = np.sqrt(np.sum(xyz**2))
            open_angle = np.arcsin(radius[name] / dist) + np.radians(limb_margin[name])
            phis = np.linspace(0, 2*np.pi, npoly)
            theta = open_angle + phis * 0.0
            sin_theta = np.sin(theta)
            ecis = np.array([np.cos(theta) * np.ones(len(phis)),
                            np.sin(phis) * sin_theta,
                            np.cos(phis) * sin_theta])  # OFLS uses -cos(phi)
            quat_x_to_obj = Ska.quatutil.quat_x_to_vec(xyz)
            obj_ecis = np.dot(quat_x_to_obj.transform, ecis)
            x, y = antisun.eci2img(obj_ecis, sun_eci)
            line = ax.plot(x, y, color=color, visible=True, alpha=0.1)[0]
            #patch = matplotlib.patches.Polygon(
            #    np.array([x,y]).transpose(), edgecolor=color, facecolor='none',
            #    fill=False, visible=True, alpha=0.1)
            # ax.add_patch(patch)
            self.patches.append(line)

    def update(self, *args):
        i0, i1 = get_index_lims()
        for i, patch in enumerate(self.patches):
            if i0 <= i <= i1:
                patch.set_visible(True)
                alpha = 1.0 - abs(i - (i0+i1)/2) / ((i1-i0+1.0)/2.0)
                print alpha
                patch.set_alpha(alpha)
            else:
                patch.set_visible(False)
                                               
class Slider(object):
    def __init__(self, minval, maxval, label_command=None, side=Tk.TOP, anchor='w', **kwargs):
        self.label_command = label_command
        
        self.frame = Tk.Frame()
        self.frame.pack(side=side, anchor=anchor)

        self.value = Tk.DoubleVar()
        self.value.set((minval + maxval) / 2.0)
        self.scale = Tk.Scale(self.frame,
                              from_=minval, to=maxval,
                              variable=self.value,
                              orient=Tk.HORIZONTAL,
                              command=self.value_changed,
                              showvalue=False,
                              **kwargs)
        self.scale.pack(side=Tk.LEFT, anchor='w')

        if label_command is not None:
            self.label_var = Tk.StringVar()
            self.label = Tk.Label(self.frame, textvariable=self.label_var)
            self.label.pack(side=Tk.LEFT)

    def value_changed(self, scaleval):
        if self.label_command is not None:
            self.label_var.set(self.label_command(self.value.get()))
        
def get_date(idx_img):
    idx_img = int(idx_img)
    if idx_img < 0:
        idx_img = 0
    if idx_img >= len(times):
        idx_img = len(times) - 1
    return DateTime(times[idx_img]).date[:-4]

# Load data
dat = pickle.load(open('cube3.pkl'))
imgs = dat['illums']
times = dat['times']
antisun = dat['antisun']
sun_eci = dat['sun_eci']
ephem_xyzs = dat['ephem_xyzs']
# Ephemerides in frame aligned with ECI but centered on Chandra
ephem_xyzs['sun'] = ephem_xyzs['solar'] - ephem_xyzs['orbit']
ephem_xyzs['earth'] = -ephem_xyzs['orbit']
ephem_xyzs['moon'] = ephem_xyzs['lunar'] - ephem_xyzs['orbit']
radius = {'sun': 695500e3,
          'moon': 1747e3,
          'earth': 6371e3}
limb_margin = dict(sun=45,
                   moon=5,
                   earth=10)

root = Tk.Tk()
root.wm_title("ESA viewer")
root.bind("<Destroy>", destroy)

fig = Figure(figsize=(8,7), dpi=100)
ax = fig.add_subplot(1, 1, 1, axis_bgcolor='k')
#ax.set_xticklabels([])
#ax.set_yticklabels([])

menu_frame = Tk.Frame()
menu_frame.pack(side=Tk.TOP, anchor='w')
quit_button = Tk.Button(master=menu_frame, text='Quit', command=sys.exit)
quit_button.pack(side=Tk.LEFT)

# Main image drawing frame (connected to Matplotlib canvas)
image_frame = Tk.Frame()
image_frame.pack(side=Tk.TOP, expand=1, fill='both')
canvas = FigureCanvasTkAgg(fig, master=image_frame)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

canvas.mpl_connect('motion_notify_event', motion)

toolbar = NavigationToolbar2TkAgg(canvas, image_frame)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

# Draw image for the first time
maxscale = 0.4 * 3 * 6   # illum = 3 hours at 0.4
image = ax.imshow(imgs[0], interpolation='bilinear', animated=True, vmin=0,
                  vmax=maxscale, alpha=1.0, origin='lower')
image.set_cmap('jet')
ax.set_autoscale_on(False)

# Make / draw Earth and Moon constraints
earth = SolarSystemObject('earth', times, ephem_xyzs, color='r')
moon = SolarSystemObject('moon', times, ephem_xyzs, color='y')

# Make colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(image, cax=cax)

# Date slider control
date_slider = Slider(minval=0, maxval=len(imgs), label_command=get_date, length=300)
date_slider.value.trace('w', earth.update)
date_slider.value.trace('w', moon.update)
date_slider.value.trace('w', update_image)

# Width slider control
width_slider = Slider(minval=0, maxval=14.0, resolution=0.1,
                      label_command=lambda x: '{0:.1f}'.format(x))
width_slider.value.trace('w', earth.update)
width_slider.value.trace('w', moon.update)
width_slider.value.trace('w', update_image)

# Alpha slider control
alpha_slider = Slider(minval=0.0, maxval=1.0, resolution=0.01,
                      label_command=lambda x: '{0:.2f}'.format(x))
alpha_slider.value.set(1.0)
alpha_slider.value.trace('w', update_image)

update_image(None)
Tk.mainloop()


#sfreq.on_changed(update)
#samp.on_changed(update)

#rax = axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
#radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)

