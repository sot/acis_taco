#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import Tkinter as Tk
import sys
import cPickle as pickle

import argparse
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import Ska.quatutil
from Chandra.Time import DateTime
import Ska.Matplotlib
from Ska.Matplotlib import cxctime2plotdate
import Ska.Sun
from Quaternion import Quat

def plot_cxctime(times, y, fmt='-b', fig=None, ax=None, yerr=None, xerr=None, tz=None,
                 state_codes=None, interactive=True, **kwargs):
    """Make a date plot where the X-axis values are in CXC time.  If no ``fig``
    value is supplied then the current figure will be used (and created
    automatically if needed).  If yerr or xerr is supplied, ``errorbar()`` will be
    called and any additional keyword arguments will be passed to it.  Otherwise
    any additional keyword arguments (e.g. ``fmt='b-'``) are passed through to
    the ``plot()`` function.  Also see ``errorbar()`` for an explanation of the possible
    forms of *yerr*/*xerr*.

    If the ``state_codes`` keyword argument is provided then the y-axis ticks and
    tick labels will be set accordingly.  The ``state_codes`` value must be a list
    of (raw_count, state_code) tuples, and is normally set to ``msid.state_codes``
    for an MSID object from fetch().

    If the ``interactive`` keyword is True (default) then the plot will be redrawn
    at the end and a GUI callback will be created which allows for on-the-fly
    update of the date tick labels when panning and zooming interactively.  Set
    this to False to improve the speed when making several plots.  This will likely
    require issuing a plt.draw() or fig.canvas.draw() command at the end.

    :param times: CXC time values for x-axis (date)
    :param y: y values
    :param fmt: plot format (default = '-b')
    :param fig: pyplot figure object (optional)
    :param yerr: error on y values, may be [ scalar | N, Nx1, or 2xN array-like ] 
    :param xerr: error on x values in units of DAYS (may be [ scalar | N, Nx1, or 2xN array-like ] )
    :param tz: timezone string
    :param state_codes: list of (raw_count, state_code) tuples
    :param interactive: use plot interactively (default=True, faster if False)
    :param **kwargs: keyword args passed through to ``plot_date()`` or ``errorbar()``

    :rtype: ticklocs, fig, ax = tick locations, figure, and axes object.
    """

    if fig is None:
        fig = matplotlib.pyplot.gcf()

    if ax is None:
        ax = fig.gca()

    if yerr is not None or xerr is not None:
        ax.errorbar(cxctime2plotdate(times), y, yerr=yerr, xerr=xerr, fmt=fmt, **kwargs)
        ax.xaxis_date(tz)
    else:
        ax.plot_date(cxctime2plotdate(times), y, fmt=fmt, **kwargs)
    ticklocs = Ska.Matplotlib.set_time_ticks(ax)
    fig.autofmt_xdate()

    if state_codes is not None:
        counts, codes = zip(*state_codes)
        ax.yaxis.set_major_locator(Ska.Matplotlib.FixedLocator(counts))
        ax.yaxis.set_major_formatter(Ska.Matplotlib.FixedFormatter(codes))

    # If plotting interactively then show the figure and enable interactive resizing
    if interactive and hasattr(fig, 'show'):
        # fig.canvas.draw()
        ax.callbacks.connect('xlim_changed', Ska.Matplotlib.remake_ticks)

    return ticklocs, fig, ax

def get_args():
    parser = argparse.ArgumentParser(description='Run Earth Solid Angle viewer')
    parser.add_argument('infile', type=str,
                       help='Input data file root (e.g. FEB1411)')
    args = parser.parse_args()
    return args

def get_input_data():
    args = get_args()
    infiles = [args.infile,
               args.infile + '.pkl',
               os.path.join(os.path.dirname(__file__), args.infile),
               os.path.join(os.path.dirname(__file__), args.infile + '.pkl')]
    for infile in infiles:
        if os.path.exists(infile):
            try:
                print 'Reading {0}'.format(infile)
                dat = pickle.load(open(infile))
                return dat
            except:
                print "ERROR: failed to load data file {0}".format(infile)
                raise
    else:
        print "ERROR: could not find input data file {0} or {0}.pkl".format(args.infile)
        sys.exit(1)

def get_date(idx_img):
    idx_img = int(idx_img)
    if idx_img < 0:
        idx_img = 0
    if idx_img >= n_times:
        idx_img = n_times - 1
    return DateTime(times[idx_img]).date[:-4]

def get_index_lims():
    center = int(date_slider.value.get())
    wid = int(width_slider.value.get() * 6) # hardcoded for 5-min intervals
    i0 = center - wid
    i1 = center + wid
    if i0 < 0:
        i0 = 0
    if i1 >= n_times:
        i1 = n_times - 1
    return i0, i1, center

class IllumImage(object):
    def __init__(self, fig):
        self.ax = fig.add_axes([0.1, 0.2, 0.8, 0.75], axisbg='k')
        self.ax.format_coord = lambda x,y: ""
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])

        # Draw image for the first time
        maxscale = 0.4 * 3 * 6   # illum = 3 hours at 0.4
        self.image = self.ax.imshow(imgs[0], interpolation='bilinear', animated=True, vmin=0,
                          vmax=maxscale, alpha=1.0, origin='lower')
        self.image.set_cmap('spectral')
        self.draw_pitch_contours()
        self.ax.set_autoscale_on(False)

        # Make colorbar
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(self.image, cax=cax)

    def update(self, *args):
        i0, i1, i_center = get_index_lims()
        i0 = imgs_idx_map[i0]
        i1 = imgs_idx_map[i1]
        img = imgs[i0:i1+1].sum(0)

        img_rgba = matplotlib.cm.spectral(img / 7.2)
        img_rgba[:, :, 3] = alpha_slider.value.get()
        self.image.set_data(img_rgba)
        self.ax.set_title(get_date(i_center))
        # image_canvas.draw()

    def draw_pitch_contours(self):
        phi = 0
        xc, yc = antisun.phys2img(0, 0)
        for pitch in (170, 156, 152, 135, 120, 105, 90, 75, 60, 45):
            x, y = antisun.phys2img(180-pitch, 0)
            rad = np.sqrt((x - xc)**2 + (y - yc)**2)
            patch = matplotlib.patches.Circle((xc, yc), rad, edgecolor='w', fill=False, linewidth=0.5)
            self.ax.add_patch(patch)

class Taco3dView(object):
    def __init__(self):
        self.window = None

    def open_window(self):
        if self.window:
            return

        # Taco3d window and mpl figure
        self.window = Tk.Toplevel()
        self.window.protocol('WM_DELETE_WINDOW', self.destroy)
        fig = Figure(figsize=(5,5), dpi=100)
        self.canvas = FigureCanvasTkAgg(fig, master=self.window)
        self.canvas.get_tk_widget().pack(side=Tk.LEFT)
        self.ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection='3d', azim=0.0, elev=0.0)
        self.draw_taco3d()
        self.canvas.show()

    def destroy(self):
        self.window.destroy()
        self.window = None
        
    def draw_taco3d(self):
        TACO_X_OFF = 250
        TACO_Y_OFF = 689
        RAD_Z_OFF = 36

        y_pnts = np.array([-689.0,  -333,  -63,  553, 689])
        z_pnts = np.array([0.0, 777, 777, 421, 0]) - RAD_Z_OFF

        y = np.array([y_pnts, y_pnts])
        z = np.array([np.zeros_like(z_pnts) - RAD_Z_OFF, z_pnts])
        x = np.zeros_like(z)

        self.ax.plot_surface(x - 250, y, z, shade=True, color='y')
        self.ax.plot_surface(x + 250, y, z, shade=True, color='r')

        y = np.linspace(-241, 241, 2)
        x = np.linspace(-216, 216, 2)
        xx, yy = np.meshgrid(x, y)
        zz = np.zeros_like(xx)

        self.ax.plot_surface(xx, yy, zz, shade=True, color='b')

        x = np.linspace(-250, 250, 2)
        y = np.linspace(-689, 689, 2)
        xx, yy = np.meshgrid(x, y)
        zz = np.zeros_like(xx) - RAD_Z_OFF
        self.ax.plot_surface(xx, yy, zz, shade=True, color='y')
        
        self.ax.set_xlim3d(-700, 700)
        self.ax.set_ylim3d(-700, 700)
        self.ax.set_zlim3d(-700, 700)
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        
    def update(self, ra, dec):
        if self.window:
            self.ax.view_init(dec, ra)
            self.canvas.draw()

class ImageCoords(object):
    def __init__(self, ax, master):
        self.ax = ax
        self.frame = Tk.Frame(master=master, relief=Tk.RAISED, borderwidth=2)
        self.textvar = dict()
        self.textvar['ra'] = Tk.StringVar()
        self.textvar['dec'] = Tk.StringVar()
        self.textvar['pitch'] = Tk.StringVar()
        self.textvar['phi'] = Tk.StringVar()
        self.textvar['earth_ra_cb'] = Tk.StringVar()
        self.textvar['earth_dec_cb'] = Tk.StringVar()
        Tk.Label(self.frame, text='RA, Dec').grid(row=0)
        Tk.Label(self.frame, text='Pitch, Phi').grid(row=1)
        Tk.Label(self.frame, text='Earth Alt, Az').grid(row=2)
        self.ra = Tk.Label(self.frame, textvariable=self.textvar['ra'], width=12)
        self.dec = Tk.Label(self.frame, textvariable=self.textvar['dec'], width=12)
        self.pitch = Tk.Label(self.frame, textvariable=self.textvar['pitch'])
        self.phi = Tk.Label(self.frame, textvariable=self.textvar['phi'])
        self.earth_ra_cb = Tk.Label(self.frame, textvariable=self.textvar['earth_ra_cb'])
        self.earth_dec_cb = Tk.Label(self.frame, textvariable=self.textvar['earth_dec_cb'])
        self.ra.grid(row=0, column=1)
        self.dec.grid(row=0, column=2)
        self.pitch.grid(row=1, column=1)
        self.phi.grid(row=1, column=2)
        self.earth_ra_cb.grid(row=2, column=2)
        self.earth_dec_cb.grid(row=2, column=1)
    
    def update(self, event):
        if event.inaxes != self.ax:
            return
        x = event.xdata
        y = event.ydata
        i_ephem = int(date_slider.value.get())
        sun_eci = ephem_xyzs['sun'][:, i_ephem]
        earth_eci = ephem_xyzs['earth'][:, i_ephem]
        ra, dec = antisun.img2sky(x, y, sun_eci)
        r, phi = antisun.img2polar(x, y)
        roll = Ska.Sun.nominal_roll(ra, dec, times[i_ephem])
        q_att = Quat([ra, dec, roll])
        earth_cb = np.dot(q_att.transform.transpose(), earth_eci)
        earth_ra_cb, earth_dec_cb = Ska.quatutil.eci2radec(earth_cb)
        pitch = 180 - r
        self.textvar['ra'].set('{0:.4f}'.format(ra))
        self.textvar['dec'].set('{0:.4f}'.format(dec))
        self.textvar['pitch'].set('{0:.1f}'.format(pitch))
        self.textvar['phi'].set('{0:.1f}'.format(np.degrees(phi)))
        self.textvar['earth_ra_cb'].set('{0:.1f}'.format(earth_ra_cb))
        self.textvar['earth_dec_cb'].set('{0:.1f}'.format(earth_dec_cb))
        taco3d_view.update(earth_ra_cb, earth_dec_cb)

class TimePlot(object):
    def __init__(self, fig, rect, times, ephem_xyzs):
        self.ax = fig.add_axes(rect)
        self.orbit_rs = np.sqrt(np.sum(ephem_xyzs['earth']**2, 0))
        plot_cxctime(times, self.orbit_rs, fig=fig, ax=self.ax)
        self.ax.grid()
        self.ax.set_autoscale_on(False)

    def get_perigee_idxs(self):
        alt = self.orbit_rs
        local_min = (alt[2:] > alt[1:-1]) & (alt[1:-1] < alt[:-2])
        return np.where(local_min)[0]

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
    RADIUS = {'sun': 695500e3,
              'moon': 1747e3,
              'earth': 6371e3}
    LIMB_MARGIN = dict(sun=45,
                       moon=6,
                       earth=10)
    def __init__(self, name, times, xyzs, color, ax, npoly=40):
        self.name = name
        self.times = times
        self.xyzs = xyzs
        self.color = color
        self.ax = ax
        self.regions = {}
        self.idxs_visible = set()

        for idx in range(n_times):
            xyz = ephem_xyzs[name][:, idx]
            dist = np.sqrt(np.sum(xyz**2))
            open_angle = (np.arcsin(self.RADIUS[name] / dist)
                          + np.radians(self.LIMB_MARGIN[name]))
            phis = np.linspace(0, 2*np.pi, npoly)
            theta = open_angle + phis * 0.0
            sin_theta = np.sin(theta)
            ecis = np.array([np.cos(theta) * np.ones(len(phis)),
                            np.sin(phis) * sin_theta,
                            np.cos(phis) * sin_theta])  # OFLS uses -cos(phi)
            quat_x_to_obj = Ska.quatutil.quat_x_to_vec(xyz)
            obj_ecis = np.dot(quat_x_to_obj.transform, ecis)
            sun_eci = ephem_xyzs['sun'][:, idx]
            x, y = antisun.eci2img(obj_ecis, sun_eci)
            self.regions[idx] = dict(x=x, y=y, linewidth=1)

    def update(self, *args):
        i0, i1, i_center = get_index_lims()
        stride = (i1 - i0) / 30 + 1
        idxs = np.arange(0, n_times, stride)
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
                region['line'] = self.ax.plot(region['x'], region['y'], linewidth=1,
                                              color=self.color, visible=True)[0]
        try:
            self.regions[self.last_idx_center]['line'].set_linewidth(1)
        except AttributeError:
            pass
        self.regions[idx_center]['line'].set_linewidth(2.5)
        self.last_idx_center = idx_center
        self.idxs_visible = idxs

class Slider(object):
    def __init__(self, minval, maxval, label_command=None,
                 side=Tk.TOP, anchor='w', master=None, **kwargs):
        self.label_command = label_command
        
        self.frame = Tk.Frame(master=master)
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
        
# Load data and set some globals (hopefully minimize this later)
dat = get_input_data()
times = dat['times']
n_times = len(times)
antisun = dat['antisun']
ephem_xyzs = dat['ephem_xyzs']
idxs = np.arange(n_times)
imgs = dat['illums']
imgs_idxs = dat['illum_idxs']
imgs_idx_map = np.searchsorted(imgs_idxs, idxs)
np.clip(imgs_idx_map, 0, n_times - 1)

# Ephemerides in frame aligned with ECI but centered on Chandra
ephem_xyzs['sun'] = ephem_xyzs['solar'] - ephem_xyzs['orbit']
ephem_xyzs['earth'] = -ephem_xyzs['orbit']
ephem_xyzs['moon'] = ephem_xyzs['lunar'] - ephem_xyzs['orbit']

taco3d_view = Taco3dView()

matplotlib.rc("axes", labelsize=10)
matplotlib.rc("xtick", labelsize=10)
matplotlib.rc("ytick", labelsize=10)

root = Tk.Tk()
root.wm_title("ESA viewer")
root.bind("<Destroy>", lambda x: sys.exit)
root.geometry('800x1000')

# Top menu bar
menu_frame = Tk.Frame(master=root)
menu_frame.pack(side=Tk.TOP, anchor='w')
quit_button = Tk.Button(master=menu_frame, text='Quit', command=sys.exit)
quit_button.pack(side=Tk.LEFT)
taco3d_button = Tk.Button(master=menu_frame, text='Taco3d', command=taco3d_view.open_window)
taco3d_button.pack(side=Tk.LEFT)

# Frame containing matplotlib figures
mpl_figs_frame = Tk.Frame(master=root)
mpl_figs_frame.pack(side=Tk.TOP, anchor='w', fill=Tk.BOTH, expand=1)

# Main image drawing frame (connected to Matplotlib canvas)
fig = Figure(figsize=(3, 4), dpi=100)
illum_image = IllumImage(fig)
time_plot = TimePlot(fig, [0.1, 0.05, 0.8, 0.15], times, ephem_xyzs)
image_frame = Tk.Frame(master=mpl_figs_frame)
image_frame.pack(side=Tk.LEFT, expand=1, fill='both')
image_canvas = FigureCanvasTkAgg(fig, master=image_frame)
image_canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

image_toolbar = NavigationToolbar2TkAgg(image_canvas, image_frame)
image_toolbar.update()
image_canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

# Make / draw Earth and Moon constraints
earth = SolarSystemObject('earth', times, ephem_xyzs, color='r', ax=illum_image.ax)
moon = SolarSystemObject('moon', times, ephem_xyzs, color='y', ax=illum_image.ax)

controls_frame = Tk.Frame(master=root)
controls_frame.pack(side=Tk.LEFT, fill=Tk.X, expand=1)
sliders_frame = Tk.Frame(master=controls_frame)
sliders_frame.pack(side=Tk.LEFT)

# Sliders
date_slider = Slider(minval=0, maxval=n_times-1, length=350, master=sliders_frame,
                     label_command=lambda x: 'Time: {0}'.format(get_date(x)))
width_slider = Slider(minval=0, maxval=14.0, resolution=0.1, length=350, master=sliders_frame,
                      label_command=lambda x: 'Width: {0:.1f} hours'.format(x))
alpha_slider = Slider(minval=0.0, maxval=1.0, resolution=0.01, master=sliders_frame,
                      label_command=lambda x: 'Alpha: {0:.2f}'.format(x))
date_slider.value.set(time_plot.get_perigee_idxs()[0])
width_slider.value.set(3.0)
alpha_slider.value.set(1.0)

image_coords = ImageCoords(illum_image.ax, master=controls_frame)
image_coords.frame.pack(side=Tk.RIGHT)

# Set up callbacks
image_canvas.mpl_connect('motion_notify_event', image_coords.update)

# trace seems to call callbacks in reverse order from initialization order
date_slider.value.trace('w', lambda *args: image_canvas.draw())
date_slider.value.trace('w', illum_image.update)
date_slider.value.trace('w', earth.update)
date_slider.value.trace('w', moon.update)
date_slider.value.trace('w', time_plot.update)

width_slider.value.trace('w', lambda *args: image_canvas.draw())
width_slider.value.trace('w', illum_image.update)
width_slider.value.trace('w', earth.update)
width_slider.value.trace('w', moon.update)
width_slider.value.trace('w', time_plot.update)

alpha_slider.value.trace('w', lambda *args: image_canvas.draw())
alpha_slider.value.trace('w', illum_image.update)

def change_date(incr):
    date_slider.value.set(date_slider.value.get()+incr)

root.bind("<Right>", lambda x: change_date(1))
root.bind("<Left>", lambda x: change_date(-1))
root.bind("<Shift-Right>", lambda x: change_date(64*12))
root.bind("<Shift-Left>", lambda x: change_date(-64*12))

# Initial updates
illum_image.update(None)
earth.update()
moon.update()
time_plot.update()
image_canvas.draw()

# Do it
Tk.mainloop()
