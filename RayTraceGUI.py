#!/usr/bin/env python
# coding: utf-8

#Version 22.10.2021

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as anim
import matplotlib
import tkinter as tk
from tkinter import ttk
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
plt.ioff()


# ### Constants


# ***************
Rmin = 2.15
Rmax = 19.
shape = 400
max_reflections = 500
decay_rate = .02
alpha_channel = .2

alpha = 4.8
beta = 0.
xoffset = 0.
yoffset = 0.
Rx = 246.0
Ry = 269.4
distance = 199.91
theta = 62.55
phi = 45.

# ***************


# ### Matrix Formalism


get_D_matrix = lambda d: np.kron(np.eye(2,dtype=int),np.eye(2) + d * np.eye(2,k=1))


#get_R_matrix = lambda Rx, Ry: lag.block_diag(np.eye(2) - 2 / Rx *  np.eye(2,k=-1), np.eye(2) - 2 / Ry * np.eye(2,k=-1))
get_R_matrix = lambda Rx, Ry: np.array([[1., 0., 0., 0.,],[-2/Rx, 1., 0., 0.], [0.,0.,1.,0.], [0., 0., -2/Ry, 1.]])


get_T_matrix = lambda theta: np.cos(theta) * np.eye(4) + np.sin(theta) * (np.eye(4,k=2) - np.eye(4,k=-2))


def hz_fast(v0,be,x,y,d,Rx,Ry,theta,phi,max_reflects = [500]):
    X = np.repeat(np.array([[x,np.sin(v0),y,np.sin(be),]]).T[None,:,:],shape,0)
    theoX = X[0,:,:].copy()
    X[:,0,0] += np.random.normal(0,0.1,shape)
    X[:,2,0] += np.random.normal(0,0.1,shape)
    T_IM = get_T_matrix(phi)
    T_IM_inv = get_T_matrix(-phi)
    X = T_IM @ X
    theoX = T_IM @ theoX
    T = get_T_matrix(theta)
    Tm = get_T_matrix(-theta)
    R = get_R_matrix(Rx,Ry)
    D = get_D_matrix(d)
    PI = R @ D
    PB = Tm @ R @ T @ D
    XIM = []
    YIM = []
    XBM = []
    YBM = []
    npass = 2
    n = 0
    for max_refl in max_reflects:
        xim = []
        yim = []
        xbm = []
        ybm = []
        while (n < max_refl):
            X = PB @ X 
            XS = T_IM_inv @ X
            xbm.extend(XS[:,0,0])
            ybm.extend(XS[:,2,0])
            alive = np.less((X[:,:1,:1]**2 + X[:,2:3,:1]**2),Rmax**2)
            X *= alive / alive
            X = PI @ X
            XS = T_IM_inv @ X
            xim.extend(XS[:,0,0])
            yim.extend(XS[:,2,0])
            alive *= np.greater((X[:,:1,:1]**2 + X[:,2:3,:1]**2),Rmin**2)
            alive *= np.less((X[:,:1,:1]**2 + X[:,2:3,:1]**2),Rmax**2)
            alive *= np.greater(np.random.rand(shape),decay_rate)[:,None,None]
            X *= alive / alive
            theoX = PI @ PB @ theoX
            if np.mean(alive) > (0.5 * decay_rate**n):
                npass += 2
                n += 1
                #theoX = PI @ PB @ theoX
            elif np.sum(alive) == 0:
                break
        XIM.append(np.array(xim))
        YIM.append(np.array(yim))
        XBM.append(np.array(xbm))
        YBM.append(np.array(ybm))
        if np.sum(alive) == 0:
            break
    return XBM, YBM, XIM, YIM,npass, (T_IM_inv @ theoX)[:,0]



def update_fig(xbs,ybs,xis,yis,):
    for implot, bmplot in zip(implots, bmplots):
        implot.set_data([],[])
        bmplot.set_data([],[])
    for implot,bmplot,xb,yb,xi,yi in zip(implots,bmplots,xbs,ybs,xis,yis):
        implot.set_data(xi,yi)
        bmplot.set_data(-xb,yb)
    fig.canvas.restore_region(background)
    for t in objs:
        fig.draw_artist(t)
    fig.canvas.blit()

def get_list_from_csv(csv):
    return [int(s) for s in csv.split(',')]

def redraw_ray():
    try:
        #v0_deg = float(alpha_var.get())
        rx = float(rx_var.get()) + rx_fin.get()
        rx_res_var.set('{:7.2f}'.format(rx))
        ry = float(ry_var.get()) + ry_fin.get()
        ry_res_var.set('{:7.2f}'.format(ry))
        d = float(dist_var.get())+ dist_fin.get()+ dist_fin2.get()
        dist_res_var.set('{:6.2f}'.format(d))
        v0_deg = float(alpha_var.get()) + alpha_fin.get()
        alpha_res_var.set('{:7.2f}'.format(v0_deg))
        th_deg = float(theta_var.get()) + theta_fin.get()
        theta_res_var.set('{:7.2f}'.format(th_deg))
        ph_deg = float(phi_var.get())+ phi_fin.get()
        phi_res_var.set('{:7.2f}'.format(ph_deg))
        x = float(x_off_var.get()) + x_off_fin.get()
        x_off_res_var.set('{:6.2f}'.format(x))
        y = float(y_off_var.get()) + y_off_fin.get()
        y_off_res_var.set('{:6.2f}'.format(y))
        be_deg = float(beta_var.get()) + beta_fin.get()
        beta_res_var.set('{:7.2f}'.format(be_deg))
    except Exception as err:
        err_var.set(err)
        return 1
    v0 = v0_deg * np.pi / 180
    be = be_deg * np.pi / 180
    th = th_deg * np.pi / 180
    ph = ph_deg * np.pi / 180
    xbs,ybs,xis,yis, npass, X = hz_fast(v0,be,x,y,d,rx,ry,th-ph,ph,max_refl_list)
    nps_var.set('{:d}'.format(npass))
    repos = X[::2]
    reang = np.arcsin(X[1::2])/np.pi*180
    repos_var.set('({:.2f} {:.2f})'.format(*repos))
    reang_var.set('({:.2f}° {:.2f}°)'.format(*reang))
    performance = np.array([np.sum(np.square(repos)), np.abs((reang[0] - v0_deg)/v0_deg), np.abs(reang[1])])
    badbeam = np.sum(np.greater(performance, [2.,1.,1.])) > 0
    goodbeam = np.sum(np.greater(performance, [.25,.1,.1])) == 0
    fg = '#0c0' if goodbeam else '#f00' if badbeam else '#000'
    nps_res.config(foreground=fg)
    update_fig(xbs,ybs,xis,yis,)
    return 0

def redraw_fig():
    global objs, implots, bmplots, background, max_refl_list
    try:
        max_refl_list = get_list_from_csv(max_var.get())
        max_refl_list.append(max_reflections)
    except Exception as err:
        err_var.set(err)
        return 1
    objs, implots, bmplots, background = reset_figure(Rmax,Rmin,max_refl_list)
    return 0

def val_changed(dummy=None):
    if redraw_ray():
        return
    err_var.set('')

def force_update(dummy=None):
    if redraw_fig():
        return 
    if redraw_ray():
        return
    err_var.set('')

# GUI

def exit_gui():
    root.destroy()
    root.quit()

root = tk.Tk()
root.title("Herriott Cell Calculator")
root.protocol("WM_DELETE_WINDOW", exit_gui)

fig,axs = plt.subplots(1,2,figsize=(10,5),dpi=150,sharey=True)
fig.canvas.draw()
empty = fig.canvas.copy_from_bbox(fig.bbox)

canvas = FigureCanvasTkAgg(fig, master=root)
plot_widget = canvas.get_tk_widget()

def reset_figure(Rmax, Rmin, max_reflects):
    fig.canvas.restore_region(empty)
    fig.set_facecolor('white')
    ic = plt.Circle((0.,0.),Rmin,fc='white',ec='k')
    for ax in axs:
        ax.clear()
        ax.set_aspect('equal')
        ax.set_facecolor('white')
        ax.set(xlabel='x',xlim=(-Rmax*1.1,Rmax*1.1),ylim=(-Rmax*1.1,Rmax*1.1))
        oc = plt.Circle((0.,0.),Rmax,fc='white',ec='k')
        ax.add_artist(oc)
    axs[0].add_artist(ic)
    axs[0].set(title='IM',ylabel='y')
    axs[1].set(title='BM')
    implots = []
    bmplots = []
    objs = []
    for j,max_refl in enumerate(max_reflects): 
        i = j+1
        implot, = axs[0].plot([],[],',',color=(i % 2, (i // 2) % 2, (i // 4) % 2, alpha_channel))
        bmplot, = axs[1].plot([],[],',',color=(i % 2, (i // 2) % 2, (i // 4) % 2, alpha_channel))
        implots.append(implot)
        bmplots.append(bmplot)
        objs.append(implot)
        objs.append(bmplot)
    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(fig.bbox)
    return objs, implots, bmplots, background

max_refl_list = [max_reflections]
objs, implots, bmplots, background = reset_figure(Rmax, Rmin, max_refl_list)

frame = ttk.Frame(master=root,)

# Geometry

rx_lbl = ttk.Label(master=frame,text='Rx [mm]')
rx_var = tk.StringVar(value='{:5.1f}'.format(Rx))
rx_ent = ttk.Entry(master=frame,textvariable=rx_var,)
rx_res_var = tk.StringVar(value='{:6.2f} mm'.format(Rx))
rx_res = ttk.Label(master=frame,textvariable=rx_res_var)
rx_fin = tk.Variable(value=0.)
rx_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=rx_fin)

ry_lbl = ttk.Label(master=frame,text='Ry [mm]')
ry_var = tk.StringVar(value='{:5.1f}'.format(Ry))
ry_ent = ttk.Entry(master=frame,textvariable=ry_var,)
ry_res_var = tk.StringVar(value='{:6.2f} mm'.format(Ry))
ry_res = ttk.Label(master=frame,textvariable=ry_res_var)
ry_fin = tk.Variable(value=0.)
ry_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=ry_fin)

dist_lbl = ttk.Label(master=frame,text='D [mm]')
dist_var = tk.StringVar(value='{:6.2f}'.format(distance))
dist_ent = ttk.Entry(master=frame,textvariable=dist_var,)
dist_res_var = tk.StringVar(value='{:6.2f} mm'.format(distance))
dist_res = ttk.Label(master=frame,textvariable=dist_res_var)
dist_fin = tk.Variable(value=0.)
dist_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=dist_fin)
dist_fin2 = tk.Variable(value=0.)
dist_scl2 = ttk.Scale(master=frame,command=val_changed,from_=-.3125,to=.3125,length=900,orient=tk.HORIZONTAL,variable=dist_fin2)

# Angles

alpha_lbl = ttk.Label(master=frame,text='alpha [°]')
alpha_var = tk.StringVar(value='{:7.2f}'.format(alpha))
alpha_ent = ttk.Entry(master=frame,textvariable=alpha_var,)
alpha_res_var = tk.StringVar(value='{:7.2f} °'.format(alpha))
alpha_res = ttk.Label(master=frame,textvariable=alpha_res_var)
alpha_fin = tk.Variable(value=0.)
alpha_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=alpha_fin)

theta_lbl = ttk.Label(master=frame,text='theta [°]')
theta_var = tk.StringVar(value='{:7.2f}'.format(theta))
theta_ent = ttk.Entry(master=frame,textvariable=theta_var,)
theta_res_var = tk.StringVar(value='{:7.2f}'.format(theta))
theta_res = ttk.Label(master=frame,textvariable=theta_res_var)
theta_fin = tk.Variable(value=0.)
theta_scl = ttk.Scale(master=frame,command=val_changed,from_=-3.,to=3.,length=200,orient=tk.HORIZONTAL,variable=theta_fin)

phi_lbl = ttk.Label(master=frame,text='phi [°]')
phi_var = tk.StringVar(value='{:7.2f}'.format(phi))
phi_ent = ttk.Entry(master=frame,textvariable=phi_var,)
phi_res_var = tk.StringVar(value='{:7.2f}'.format(phi))
phi_res = ttk.Label(master=frame,textvariable=phi_res_var)
phi_fin = tk.Variable(value=0)
phi_scl = ttk.Scale(master=frame,command=val_changed,from_=-3.,to=3.,length=200,orient=tk.HORIZONTAL,variable=phi_fin)

# Errors

beta_lbl = ttk.Label(master=frame,text='beta [°]')
beta_var = tk.StringVar(value='{:7.2f}'.format(beta))
beta_ent = ttk.Entry(master=frame,textvariable=beta_var,)
beta_res_var = tk.StringVar(value='{:7.2f} °'.format(beta))
beta_res = ttk.Label(master=frame,textvariable=beta_res_var)
beta_fin = tk.Variable(value=0.)
beta_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=beta_fin)

x_off_lbl = ttk.Label(master=frame,text='x0 [mm]')
x_off_var = tk.StringVar(value='{:6.2f}'.format(xoffset))
x_off_ent = ttk.Entry(master=frame,textvariable=x_off_var,)
x_off_res_var = tk.StringVar(value='{:6.2f} mm'.format(xoffset))
x_off_res = ttk.Label(master=frame,textvariable=x_off_res_var)
x_off_fin = tk.Variable(value=0.)
x_off_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=x_off_fin)

y_off_lbl = ttk.Label(master=frame,text='y0 [mm]')
y_off_var = tk.StringVar(value='{:6.2f}'.format(yoffset))
y_off_ent = ttk.Entry(master=frame,textvariable=x_off_var,)
y_off_res_var = tk.StringVar(value='{:6.2f} mm'.format(yoffset))
y_off_res = ttk.Label(master=frame,textvariable=y_off_res_var)
y_off_fin = tk.Variable(value=0.)
y_off_scl = ttk.Scale(master=frame,command=val_changed,from_=-2.,to=2.,length=200,orient=tk.HORIZONTAL,variable=y_off_fin)


# Results

max_lbl = ttk.Label(master=frame,text='Color Splits')
max_var = tk.StringVar(value='')
max_ent = ttk.Entry(master=frame,textvariable=max_var,)

nps_lbl = ttk.Label(master=frame,text='Num Passes')
nps_var = tk.StringVar()
nps_res = ttk.Label(master=frame,textvariable=nps_var,font=('Helvetica', 10, 'bold',))

repos_lbl = ttk.Label(master=frame,text='Exit beam pos')
repos_var = tk.StringVar()
repos_res = ttk.Label(master=frame,textvariable=repos_var,)

reang_lbl = ttk.Label(master=frame,text='Exit beam angle')
reang_var = tk.StringVar()
reang_res = ttk.Label(master=frame,textvariable=reang_var,)



meta_btn = ttk.Button(master=frame,text='Update',command=force_update)

err_var = tk.StringVar()
err_lbl = ttk.Label(master=frame,textvariable=err_var)

frame.pack(padx=5,pady=5,side=tk.TOP)
plot_widget.pack(padx=5,pady=5,side=tk.TOP,)


dist_scl2.grid(column=0,row=6,columnspan=12,padx=5,pady=5)

col = 0
rx_lbl.grid(column=0+col,row=0,padx=5,pady=5)
rx_ent.grid(column=1+col,row=0,padx=5,pady=5)
rx_res.grid(column=2+col,row=0,padx=5,pady=5)
rx_scl.grid(column=0+col,row=1,columnspan=3,padx=5,pady=5)

ry_lbl.grid(column=0+col,row=2,padx=5,pady=5)
ry_ent.grid(column=1+col,row=2,padx=5,pady=5)
ry_res.grid(column=2+col,row=2,padx=5,pady=5)
ry_scl.grid(column=0+col,row=3,columnspan=3,padx=5,pady=5)

dist_lbl.grid(column=0+col,row=4,padx=5,pady=5)
dist_ent.grid(column=1+col,row=4,padx=5,pady=5)
dist_res.grid(column=2+col,row=4,padx=5,pady=5)
dist_scl.grid(column=0+col,row=5,columnspan=3,padx=5,pady=5)


col = 3

alpha_lbl.grid(column=0+col,row=0,padx=5,pady=5)
alpha_ent.grid(column=1+col,row=0,padx=5,pady=5)
alpha_res.grid(column=2+col,row=0,padx=5,pady=5)
alpha_scl.grid(column=0+col,row=1,columnspan=3,padx=5,pady=5)

theta_lbl.grid(column=0+col,row=2,padx=5,pady=5)
theta_ent.grid(column=1+col,row=2,padx=5,pady=5)
theta_res.grid(column=2+col,row=2,padx=5,pady=5)
theta_scl.grid(column=0+col,row=3,columnspan=3,padx=5,pady=5)

phi_lbl.grid(column=0+col,row=4,padx=5,pady=5)
phi_ent.grid(column=1+col,row=4,padx=5,pady=5)
phi_res.grid(column=2+col,row=4,padx=5,pady=5)
phi_scl.grid(column=0+col,row=5,columnspan=3,padx=5,pady=5)

col = 6

beta_lbl.grid(column=0+col,row=0,padx=5,pady=5)
beta_ent.grid(column=1+col,row=0,padx=5,pady=5)
beta_res.grid(column=2+col,row=0,padx=5,pady=5)
beta_scl.grid(column=0+col,row=1,columnspan=3,padx=5,pady=5)

x_off_lbl.grid(column=0+col,row=2,padx=5,pady=5)
x_off_ent.grid(column=1+col,row=2,padx=5,pady=5)
x_off_res.grid(column=2+col,row=2,padx=5,pady=5)
x_off_scl.grid(column=0+col,row=3,columnspan=3,padx=5,pady=5)

y_off_lbl.grid(column=0+col,row=4,padx=5,pady=5)
y_off_ent.grid(column=1+col,row=4,padx=5,pady=5)
y_off_res.grid(column=2+col,row=4,padx=5,pady=5)
y_off_scl.grid(column=0+col,row=5,columnspan=3,padx=5,pady=5)

col = 9

max_lbl.grid(column=0+col,row=0,padx=5,pady=5)
max_ent.grid(column=1+col,row=0,padx=5,pady=5,columnspan=2)

meta_btn.grid(column=1+col,row=1,padx=5,pady=5,)

err_lbl.grid(column=0+col,row=2,columnspan=3,padx=5,pady=5)

nps_lbl.grid(column=0+col,row=3,padx=5,pady=5)
nps_res.grid(column=1+col,row=3,columnspan=2,padx=5,pady=5)

repos_lbl.grid(column=0+col,row=4,padx=5,pady=5)
repos_res.grid(column=1+col,row=4,columnspan=2,padx=5,pady=5)

reang_lbl.grid(column=0+col,row=5,padx=5,pady=5)
reang_res.grid(column=1+col,row=5,columnspan=2,padx=5,pady=5)



val_changed()

root.mainloop()