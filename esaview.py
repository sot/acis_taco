#!/usr/bin/env python
import cPickle as pickle
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from Chandra.Time import DateTime

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

import Tkinter as Tk
import sys

def destroy(e): sys.exit()

def update_image(val):
    center = int(date_scale.get())
    wid = int(width_scale.get() * 6)
    i0 = center - wid
    i1 = center + wid
    if i0 < 0:
        i0 = 0
    if i1 >= len(imgs):
        i1 = len(imgs)-1
    i_center = (i0 + i1) // 2
    date_label_var.set(DateTime(times[i_center]).date[:-4])
    width_label_var.set("%.1f" % width_scale.get())
    image.set_data(imgs[i0:i1+1].sum(0))
    canvas.draw()

def motion(event):
    if event.inaxes is None:
        return
    #print event.name, event.xdata, event.ydata, event.inaxes
    
class Slider(object):
    def __init__(self, minval, maxval, command=None, side=Tk.TOP, anchor='w', **kwargs):
        self.frame = Tk.Frame()
        self.frame.pack(side=side, anchor=anchor)

        self.scale = Tk.Scale(self.frame, from_=minval, to=maxval,
                              variable=(minval + maxval) / 2.0,
                              orient=Tk.HORIZONTAL,
                              command=self.update,
                              showvalue=False,
                              **kwargs)
        self.scale.pack(side=Tk.LEFT, anchor='w')
        self.scale.set(50)

        self.label_var = Tk.StringVar()
        self.label = Tk.Label(self.frame, textvariable=self.label_var)
        self.label.pack(side=Tk.LEFT)

    def update(self, scaleval):
        pass
        
# Load data
dat = pickle.load(open('cube.pkl'))
imgs = dat['illums']
ephem_rs = dat['ephem_r']
times = dat['times']

root = Tk.Tk()
root.wm_title("ESA viewer")
root.bind("<Destroy>", destroy)

fig = Figure(figsize=(8,7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.set_xticklabels([])
ax.set_yticklabels([])

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

# Date slider control
date_slider_frame = Tk.Frame()
date_slider_frame.pack(side=Tk.TOP, anchor='w')

date_scale = Tk.Scale(date_slider_frame, from_=0, to=len(imgs), variable=50.,
                      orient=Tk.HORIZONTAL, command=update_image, showvalue=False,
                      length=300)
date_scale.pack(side=Tk.LEFT, anchor='w')
date_scale.set(50)

date_label_var = Tk.StringVar()
date_label = Tk.Label(date_slider_frame, textvariable=date_label_var)
date_label.pack(side=Tk.LEFT)

# Width slider control
width_slider_frame = Tk.Frame()
width_slider_frame.pack(side=Tk.TOP, anchor='w')

width_scale = Tk.Scale(width_slider_frame, from_=0, to=14.0, resolution=0.1, variable=5.,
                      orient=Tk.HORIZONTAL, command=update_image, showvalue=False)
width_scale.pack(side=Tk.LEFT, anchor='w')
width_scale.set(3.0)

width_label_var = Tk.StringVar()
width_label_var.set('hello world')
width_label = Tk.Label(width_slider_frame, textvariable=width_label_var)
width_label.pack(side=Tk.LEFT)

maxscale = 0.4 * 3 * 6   # illum = 3 hours at 0.4
image = ax.imshow(imgs[50], interpolation='bilinear', animated=True, vmin=0, vmax=maxscale, alpha=1.0)
image.set_cmap('jet')
ax.get_xaxis().set_ticklabels([])
ax.get_yaxis().set_ticklabels([])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(image, cax=cax)

update_image(None)
Tk.mainloop()


#sfreq.on_changed(update)
#samp.on_changed(update)

#rax = axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
#radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)

