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

#def update_alpha(*args):
#    image.set_alpha(alpha_slider.value.get())
#    canvas.draw()

def update_image(*args):
    center = int(date_slider.value.get())
    wid = int(width_slider.value.get() * 6)
    i0 = center - wid
    i1 = center + wid
    if i0 < 0:
        i0 = 0
    if i1 >= len(imgs):
        i1 = len(imgs)-1
    i_center = (i0 + i1) // 2
    img = imgs[i0:i1+1].sum(0)
    if 1:
        img_rgba = matplotlib.cm.jet(img / 7.2)
        img_rgba[:, :, 3] = alpha_slider.value.get()
        image.set_data(img_rgba)
    else:
        image.set_data(img)
    canvas.draw()

def motion(event):
    if event.inaxes is None:
        return
    #print event.name, event.xdata, event.ydata, event.inaxes
    
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
date_slider = Slider(minval=0, maxval=len(imgs), label_command=get_date, length=300)
date_slider.value.trace('w', update_image)

# Width slider control
width_slider = Slider(minval=0, maxval=14.0, label_command=lambda x: '{0:.1f}'.format(x))
width_slider.value.trace('w', update_image)

# Alpha slider control
alpha_slider = Slider(minval=0.0, maxval=1.0, resolution=0.01,
                      label_command=lambda x: '{0:.2f}'.format(x))
alpha_slider.value.trace('w', update_image)

maxscale = 0.4 * 3 * 6   # illum = 3 hours at 0.4
image = ax.imshow(imgs[50], interpolation='nearest', animated=True, vmin=0, vmax=maxscale, alpha=1.0)
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

