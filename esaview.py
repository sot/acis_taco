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
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate

def destroy(e): sys.exit()

#def update_alpha(*args):
#    image.set_alpha(alpha_slider.value.get())
#    canvas.draw()

def get_index_lims():
    center = int(date_slider.value.get())
    wid = int(width_slider.value.get() * 3)
    i0 = center - wid
    i1 = center + wid
    if i0 < 0:
        i0 = 0
    if i1 >= len(imgs):
        i1 = len(imgs) - 1
    return i0, i1, center

def update_image(*args):
    i0, i1, i_center = get_index_lims()
    img = imgs[i0:i1+1].sum(0)
    img_rgba = matplotlib.cm.spectral(img / 7.2)
    img_rgba[:, :, 3] = alpha_slider.value.get()
    image.set_data(img_rgba)
    ax.set_title(get_date(i_center))
    canvas.draw()


class ImageCoords(object):
    def __init__(self):
        self.frame = Tk.Frame()
        self.textvar = dict()
        self.textvar['ra'] = Tk.StringVar()
        self.textvar['dec'] = Tk.StringVar()
        self.textvar['pitch'] = Tk.StringVar()
        self.textvar['phi'] = Tk.StringVar()
        Tk.Label(self.frame, text='RA, Dec').grid(row=0)
        Tk.Label(self.frame, text='Pitch, Phi').grid(row=1)
        self.ra = Tk.Label(self.frame, textvariable=self.textvar['ra'], width=12)
        self.dec = Tk.Label(self.frame, textvariable=self.textvar['dec'], width=12)
        self.pitch = Tk.Label(self.frame, textvariable=self.textvar['pitch'])
        self.phi = Tk.Label(self.frame, textvariable=self.textvar['phi'])
        self.ra.grid(row=0, column=1)
        self.dec.grid(row=0, column=2)
        self.pitch.grid(row=1, column=1)
        self.phi.grid(row=1, column=2)
    
    def update(self, event):
        if event.inaxes != ax:
            return
        x = event.xdata
        y = event.ydata
        ra, dec = antisun.img2sky(x, y, sun_eci)
        r, phi = antisun.img2polar(x, y)
        pitch = 180 - r
        self.textvar['ra'].set('{0:.4f}'.format(ra))
        self.textvar['dec'].set('{0:.4f}'.format(dec))
        self.textvar['pitch'].set('{0:.1f}'.format(pitch))
        self.textvar['phi'].set('{0:.1f}'.format(np.degrees(phi)))

class TimePlot(object):
    def __init__(self, fig, rect, times, ephem_xyzs):
        self.ax = fig.add_axes(rect)
        orbit_rs = np.sqrt(np.sum(ephem_xyzs['earth']**2, 0))
        plot_cxctime(times, orbit_rs, fig=fig, ax=self.ax)
        self.ax.grid()
        self.ax.set_autoscale_on(False)
        self.update()

    def update(self, *args):
        pd0, pd1, pd_center = cxctime2plotdate(times[np.array(get_index_lims())])
        if not hasattr(self, 'patch'):
            self.patch = matplotlib.patches.Rectangle(
                (pd0, 0), width=(pd1-pd0), height=1.6e8, zorder=-100,
                facecolor='y', alpha=0.5)
            self.ax.add_patch(self.patch)
        else:
            self.patch.set_xy((pd0, 0))
            self.patch.set_width(pd1-pd0)

class SolarSystemObject(object):
    def __init__(self, name, times, xyzs, color, ax, npoly=40):
        self.name = name
        self.times = times
        self.xyzs = xyzs
        self.color = color
        self.ax = ax
        self.regions = {}
        self.idxs_visible = set()

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
            self.regions[idx] = dict(x=x, y=y, linewidth=1)

    def update(self, *args):
        i0, i1, i_center = get_index_lims()
        stride = (i1 - i0) / 20
        if stride == 0:
            stride = 1
        idxs = np.arange(0, len(self.times), stride)
        idx_center = idxs[np.argmin(np.abs(idxs - i_center))]
        idxs = set(idxs[(idxs >= i0) & (idxs <= i1)])
        # Disable regions that are currently visible but not in next view
        for idx in self.idxs_visible - idxs:
            self.regions[idx]['line'].set_visible(False)
            
        for idx in idxs - self.idxs_visible:
            region = self.regions[idx]
            try:
                region['line'].set_visible(True)
            except KeyError:
                region['line'] = self.ax.plot(region['x'], region['y'],
                                              color=self.color, visible=True)[0]

        self.idxs_visible = idxs

#         for i, line in enumerate(self.lines):
#             if i0 <= i <= i1:
#                 line.set_visible(True)
#                 if i == i_center:
#                     line.set_linewidth(2.0)
#                 else:
#                     line.set_linewidth(1.0)
#             else:
#                 line.set_visible(False)
                                               
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

# Load data and set some globals (hopefully minimize this later)
filename = sys.argv[1]
dat = pickle.load(open(filename))
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

matplotlib.rc("axes", labelsize=10)
matplotlib.rc("xtick", labelsize=10)
matplotlib.rc("ytick", labelsize=10)

fig = Figure(figsize=(8,9), dpi=100)
ax = fig.add_axes([0.1, 0.25, 0.7, 0.7], axisbg='b')
ax.format_coord = lambda x,y: ""

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

toolbar = NavigationToolbar2TkAgg(canvas, image_frame)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

# Draw image for the first time
maxscale = 0.4 * 3 * 6   # illum = 3 hours at 0.4
image = ax.imshow(imgs[0], interpolation='bilinear', animated=True, vmin=0,
                  vmax=maxscale, alpha=1.0, origin='lower')
image.set_cmap('spectral')
ax.set_autoscale_on(False)

# Make / draw Earth and Moon constraints
earth = SolarSystemObject('earth', times, ephem_xyzs, color='r', ax=ax)
moon = SolarSystemObject('moon', times, ephem_xyzs, color='y', ax=ax)

# Make colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(image, cax=cax)

# Date slider control
date_slider = Slider(minval=0, maxval=len(imgs)-1, label_command=get_date, length=500)
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

time_plot = TimePlot(fig, [0.1, 0.05, 0.7, 0.15], times, ephem_xyzs)
date_slider.value.trace('w', time_plot.update)
width_slider.value.trace('w', time_plot.update)

image_coords = ImageCoords()
image_coords.frame.pack()
canvas.mpl_connect('motion_notify_event', image_coords.update)

update_image(None)
earth.update()
moon.update()
Tk.mainloop()
