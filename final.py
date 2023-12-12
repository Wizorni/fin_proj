from tkinter import *
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from astropy.io import fits
from math import sqrt

def rabota():

    Radius = 0
    rad_fon = []
    star_cent = []
    noise = 0
    eng = 0
    fin = 0
    sc = vod2.get()
    sc = sc.split(' ')
    way = ''

    way = vod1.get()
    fot = fits.open(way)
    scidata = fot[0].data
    exptime = float(fot[0].header['exptime'])
    for i in range(len(sc)):
        star_cent.append(int(sc[i]))

    Radius = int(vod3.get())

    rr = vod4.get()
    rr = rr.split(' ')
    for i in range(len(rr)):
        rad_fon.append(int(rr[i]))
    if rad_fon[0] < rad_fon[1]:
        rad_fon[0], rad_fon[1] = rad_fon[1], rad_fon[0]

    summa = 0
    s_pix = 0
    for i in range(-Radius, Radius):
        for j in range(-Radius, Radius):
            if sqrt((i**2) + (j**2)) <= Radius:
                summa += float(scidata[star_cent[1]+j, star_cent[0]+i])
                s_pix += 1
    n_pix = 0
    for i in range(-rad_fon[1], rad_fon[1]):
        for j in range(-rad_fon[0], rad_fon[0]):
                if sqrt((i**2) + (j**2)) <= rad_fon[0] and sqrt((i**2) + (j**2)) >= rad_fon[1]:
                    noise += float(scidata[star_cent[1]+j, star_cent[0]+i])
                    n_pix += 1
    fin = str((summa / exptime) - (noise * s_pix/(exptime * n_pix)))

    lb_fin['text'] = fin

def graf():
    global exptime, star_cent, Radius, rad_fon, way, noise, eng, fin
    global fot, scidata, exptime

    Radius = 0
    star_cent = [0, 0]
    way = ''
    rad_fon = []

    way = vod1.get()
    fot = fits.open(way)
    scidata = fot[0].data

    rr = vod4.get()
    rr = rr.split(' ')
    for i in range(len(rr)):
        rad_fon.append(int(rr[i]))
    if rad_fon[0] < rad_fon[1]:
        rad_fon[0], rad_fon[1] = rad_fon[1], rad_fon[0]

    sc = vod2.get()
    sc = sc.split(' ')
    for i in range(len(sc)):
        star_cent[i] = int(sc[i])

    Radius = int(vod3.get())

    x_star = []
    en_x = []
    y_star = []
    en_y = []

    if ch_V.get() == 1:
        for i in range(star_cent[0] - rad_fon[0], star_cent[0] + rad_fon[0]):
            en_x.append(scidata[star_cent[1], i])
            x_star.append(i)
        plt.figure()
        plt.plot(x_star, en_x)
        plt.xlim(star_cent[0] - rad_fon[0], star_cent[0] + rad_fon[0])
        plt.ylabel('Energy')
        plt.xlabel("X coord")
        plt.title('2D (by X axis)')
        plt.show()

    elif ch_V.get() == 2:
        for j in range(star_cent[1] - rad_fon[0], star_cent[1] + rad_fon[0]):
            en_y.append(scidata[j, star_cent[0]])
            y_star.append(j)
        plt.figure()
        plt.plot(y_star, en_y)
        plt.xlim(star_cent[1] - rad_fon[0], star_cent[1] + rad_fon[0])
        plt.ylabel('Energy')
        plt.xlabel("Y coord")
        plt.title('2D (by Y axis)')
        plt.show()

    elif ch_V.get() == 3:
        for i in range(star_cent[0] - rad_fon[0], star_cent[0] + rad_fon[0]):
            x_star.append(i)
        for j in range(star_cent[1] - rad_fon[0], star_cent[1] + rad_fon[0]):
            y_star.append(j)

        ful_en = np.zeros([len(y_star), len(x_star)], dtype = float)
        for i in range(star_cent[0]-rad_fon[0],star_cent[0] + rad_fon[0]):
            for j in range(star_cent[1]-rad_fon[0], star_cent[1] + rad_fon[0]):
                ful_en[j - (star_cent[1]-rad_fon[0])][i - (star_cent[0]-rad_fon[0])] = float(scidata[ j][i])
        print(ful_en)
        x_plot, y_plot = np.meshgrid(x_star, y_star)
        fig = (plt.figure())
        ax = fig.add_subplot(projection='3d')
        surf = ax.plot_surface(x_plot, y_plot, ful_en, cmap = cm.plasma)
        ax.set_xlabel('X coord')
        ax.set_ylabel('Y coord')
        ax.set_zlabel('Energy')
        plt.title('3d')
        plt.show()


    R = 0
    star_cent = [0, 0]

#############################################################################################################################################
#way = ''
#Radius = 0
#rad_fon = []
#star_cent = []
#noise = 0
#fin = 0
#############################################################################################################################################
root = Tk()
for c in range(3):
    root.columnconfigure(index=c)
for r in range(20):
    root.rowconfigure(index=r)
root.title("Calc of brightness")
root.geometry("1000x400")
root.resizable(False, True)

lb1 = ttk.Label(text = 'ВВедите расположение фитки', font = ("Smeshariki 2007 fixed", 12))
lb1.grid(row = 0, column =0)
vod1 = ttk.Entry()
vod1.grid(row = 0, column =2)

lb2 = ttk.Label(text = 'ВВедите координаты центра звезды (X,Y через пробел)', font = ("Smeshariki 2007 fixed", 12))
lb2.grid(row = 2, column =0)
vod2 = ttk.Entry()
vod2.grid(row = 2, column =2)

lb3 = ttk.Label(text = 'ВВедите радиус звезды в пикселях', font = ("Smeshariki 2007 fixed", 12))
lb3.grid(row = 4, column =0)
vod3 = ttk.Entry()
vod3.grid(row = 4, column =2)

lb4 = ttk.Label(text = 'ВВедите внешний и внутренний радиус фона звезды через пробел', font = ("Smeshariki 2007 fixed", 12))
lb4.grid(row = 6, column =0)
vod4 = ttk.Entry()
vod4.grid(row = 6, column =2)

bt_fin = ttk.Button(text = 'Enter', command=rabota)
bt_fin.grid(row = 8, column = 2)

lb_fin = ttk.Label(text = '                      ', font=("Smeshariki 2007 fixed", 12))
lb_fin.grid(row = 10, column = 1)

lb_quest = ttk.Label(text = 'Какие графики нужны?', font = ("Smeshariki 2007 fixed", 12))
lb_quest.grid(row = 12, column = 0, columnspan = 3)

ch_V = IntVar()
ch1= ttk.Radiobutton(text ="2D, X", variable=ch_V, value= 1)
ch1.grid(row= 14, column=0)
ch2= ttk.Radiobutton(text ="2D, Y", variable=ch_V, value= 2)
ch2.grid(row= 14, column=1)
ch3= ttk.Radiobutton(text ="3D", variable=ch_V, value= 3)
ch3.grid(row= 14, column=2)
bt_gr = ttk.Button(text='Enter', command=graf)
bt_gr.grid(row = 16, column = 1)

root.mainloop()
###############################################################################################################################################