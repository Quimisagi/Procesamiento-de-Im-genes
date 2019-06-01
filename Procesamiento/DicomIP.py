
import warnings
import pydicom
import math

import os
import numpy as np
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import random

from tkinter import *
from tkinter import filedialog
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import messagebox
from PIL import ImageTk, Image

intensity = [0]*2500
image = []
rows = 0
columns = 0


def loadFiles():
    #resetSelector()
    fileName.set(filedialog.askopenfilename(filetypes=[("Dicom files", "*.dcm")]))
       
    lstFilenames = []
    processImage()
    

def resetSelector():
    global lstFilesDCM
    lstFilesDCM.clear()
    info_text.delete('1.0', END)
    file_cb.set("Seleccionar Archivo")
    functions_cb.set("Funcion a aplicar")
    apply_bt.configure(state=DISABLED)    
    
def showHeaderInfo(header):
    info_text.delete('1.0', END)
    info_text.config(state=NORMAL)  
    info = "\nPatient ID: "+header.PatientID
    info += "\nManufacturer: "+header.Manufacturer
    info += "\nStudy Description: "+header.StudyDescription
    info += "\nMR Acquisition Type: "+header.MRAcquisitionType
    info += "\nSpacing Between Slices: "+str(header.SpacingBetweenSlices)
    info += "\nPixel Bandwidth: "+str(header.PixelBandwidth)
    info += "\nRows: "+str(header.Rows)
    info += "\nColumns: "+str(header.Columns)
    info += "\nPixel Spacing : "+str(header.PixelSpacing)
    info_text.delete('1.0', END)
    info_text.insert('1.0', info)
    info_text.config(state=DISABLED)
    
def hist():
    global canvas, intensity, image, rows, columns
    ds = pydicom.read_file(fileName.get())
    rows = int(ds.Rows)
    columns = int(ds.Columns)

    for i in range(rows):
	    for j in range(columns):
		    intensity[image[i][j]] += 1
            
	
    intensity = np.asarray(intensity)

    
    plt.plot(intensity)
    plt.title("Histograma")
    plt.show()


def convu(pollo):
    global canvas, image   

    kernel = pollo

    ds = pydicom.read_file(fileName.get())
    rows = int(ds.Rows)
    columns = int(ds.Columns)
    pixArray = ds.pixel_array

    count = 0

    convu = np.zeros((rows-2,columns-2))


    for y in range(1,rows-1):
        for x in range(1, columns-1):
             count = 0
             for i in range(3):
                for j in range(3):
                    
                    count += kernel[i][j] * image[y+i-1][x+j-1]

             #print(count)
             convu[y-1][x-1] = count
             image[y-1][x-1] = count
             
    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(convu))     
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=RIGHT) 


def gauss(N = 3, sigma = 5):

    
    gaussian = np.zeros((N,N))
    total = 0

    for y in range (-1,2):
        for x in range(-1,2):
            #We write down the left side (scalar value) of the equation: (1/ 2 * PI * SIGMA)
            left_side = (1 / (np.pi * 2 * np.power(sigma, 2)))
            #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
            aux = -(np.power(x, 2) + np.power(y, 2)) / (2 * np.power(sigma, 2))
            #Finally, we use the Exp function to calculate the e^(aux)
            right_side = np.exp(aux)
            total += left_side * right_side
            gaussian[y+1][x+1] = left_side * right_side
    pan = 1/total
    
    for y in range (3):
        for x in range(3):
            gaussian[x][y] =  gaussian[x][y] * pan

    print(gaussian)
    convu(gaussian)

def rayleigh(N = 3, sigma = 8):
    
    rayMat = np.zeros((N,N))
    total = 0

    for y in range (3):
        for x in range(3):
            #We write down the left side (scalar value) of the equation: (1/ 2 * PI * SIGMA)
            left_side = np.divide(x * y, np.power(sigma, 4))
            #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
            aux = -(np.power(x, 2) + np.power(y, 2)) / (2 * np.power(sigma, 2))
            #Finally, we use the Exp function to calculate the e^(aux)
            right_side = np.exp(aux)
            total += left_side * right_side
            rayMat[y][x] = left_side * right_side


    print(rayMat)
    pan = 1/total
    
    for y in range (3):
        for x in range(3):
            rayMat[x][y] =  rayMat[x][y] * pan

    convu(rayMat)


def median():
    global canvas, image

    ds = pydicom.read_file(fileName.get())
    rows = int(ds.Rows)
    columns = int(ds.Columns)
    image = ds.pixel_array
    members = [(0,0)] * 9

    for i in range(rows-1):
        for j in range(columns-1):
            if i==0 or i==rows-1 or j==0 or j==columns-1:
                image[i,j]=0
            else:
                members[0] = image[i-1][j-1]
                members[1] = image[i-1][j]
                members[2] = image[i-1][j+1]
                members[3] = image[i][j-1]
                members[4] = image[i][j]
                members[5] = image[i][j+1]
                members[6] = image[i+1][j-1]
                members[7] = image[i+1][j]
                members[8] = image[i+1][j+1]
                members.sort()
                image[i][j] = members[4]

    
    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(image))     
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=RIGHT) 


def sobel():
    global canvas, image

    ds = pydicom.read_file(fileName.get())
    rows = int(ds.Rows)
    columns = int(ds.Columns)
    image = ds.pixel_array

    sobelMat = np.zeros((rows-2,columns-2))

    gradX = [[-1, 0, 1],
             [-2, 0, 2],
             [-1, 0, 1]]

    gradY = [[-1, -2, -1],
             [-0, 0, 0],
             [1, 2, 1]]

    for y in range(1, rows-1):
        for x in range(1, columns-1):
            Gx= 0
            Gy= 0
            for i in range(3):
                for j in range(3):
                    Gx += gradX[i][j] * image[y+i-1][x+j-1]
                    Gy += gradY[i][j] * image[y+i-1][x+j-1]

            image[y-1][x-1] = math.sqrt((Gx * Gx) + (Gy * Gy))



    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(image))     
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=RIGHT) 

def otsu():
    global canvas, intensity, image

    val_max = -99999
    thr = -1

    for t in range(1,2500):
        # Non-efficient implementation
        q1 = np.sum(intensity[:t])
        q2 = np.sum(intensity[t:])
        m1 = np.sum(np.array([i for i in range(t)])*intensity[:t])/q1
        m2 = np.sum(np.array([i for i in range(t,2500)])*intensity[t:])/q2
        val = q1*(1-q1)*np.power(m1-m2,2)
        if val_max < val:
            val_max = val
            thr = t

    print("Threshold: {}".format(thr))

    for i in range (rows):
        for j in range (columns):
            if image[i][j] > thr:
                image[i][j] = 0
            else:
                image[i][j] = 2499

    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(image))     
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=RIGHT) 


def kmeans(k = 5):
    global canvas, image

    ds = pydicom.read_file(fileName.get())
    rows = int(ds.Rows)
    columns = int(ds.Columns)

    centroids = [0] * k
    data = [0] * rows*columns
    clusters = [0] * k
    clustTotal = [0] * k


    for i in range (k):
        pollo = random.randint(40, 900)
        centroids[i] = pollo
        temp = 0

        clusters[i] = []

    finalImage = [0] * rows
    centroids[0] = 20

    for i in range (rows):
        finalImage[i] = [0] * columns
        for j in range (columns):
            data[temp] = image[i][j]
            temp += 1
            finalImage[i][j] = [0, 0, 0]

    data.sort()

    isSalami = False

    while isSalami == False:
        for i in range (k):
            clusters[i] = []
            clustTotal[i] = 0
        for i in range(rows*columns):
            nearest = 100000
            queso = 0
            for j in range(k):
            
                if abs(data[i] - centroids[j]) <= nearest:
                

                    nearest = abs(data[i] - centroids[j])
                    queso = j
            
            clusters[queso].append(data[i])
            clustTotal[queso] += data[i]

        for a in range(k):
            pollo = round(clustTotal[a]/len(clusters[a]))
            if pollo != centroids[a]:
                isSalami = False
                centroids[a] = pollo
            else:
                isSalami = True

    print(len(clusters[0]))
    print(len(clusters[1]))

    for i in range (rows):
        for j in range (columns):
            for a in range(k):
                if image[i][j] >= clusters[a][0] and image[i][j] <= clusters[a][len(clusters[a])-1]:
                    #image[i][j] = centroids[a]

                    if a == 0:
                        finalImage[i][j] = [0,0,0]
                    if a == 1:
                        finalImage[i][j] = [255,0,0]
                    if a == 2:
                        finalImage[i][j] = [0,255,0]
                    if a == 3:
                        finalImage[i][j] = [0,0,255]
                    if a == 4:
                        finalImage[i][j] = [255,255,30]



   


    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(finalImage))     
    #.imshow(np.flipud(finalImage))  
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=RIGHT)





def processImage():
    global canvas, intensity, image, rows, columns       
    text_fr.pack(padx=10, pady=10, side=RIGHT)
    histo_bt.configure(state=NORMAL)
    ds = pydicom.read_file(fileName.get())
    showHeaderInfo(ds)

    image = ds.pixel_array
    intensity = [0] * 2500

    rows = int(ds.Rows)
    columns = int(ds.Columns)

    for i in range(rows):
	    for j in range(columns):
		    intensity[image[i][j]] += 1
            
	
    intensity = np.asarray(intensity)

    plt.set_cmap(plt.gray())
    f = Figure()
    a = f.add_subplot(111)    
    a.imshow(np.flipud(ds.pixel_array))     
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, images_fr)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP)  
    
 


#TK components
root = tk.Tk()
root.title("Image Proccesing")
root.configure(background='#f8f8f8')
root.geometry('%dx%d+%d+%d' % (1000, 720, 5, 5))

bannerPath = "banner.png"
img = ImageTk.PhotoImage(Image.open(bannerPath))
panel = tk.Label(root, image = img)
panel.configure(background='#f8f8f8')
panel.pack(side="top", pady=3, fill = "both", expand = "no")


applyEfct_fr = Frame(root, width=100, height=450)
listEffcts_fr = Frame(applyEfct_fr, relief=GROOVE, borderwidth=2)
listEffcts_fr.configure(background='#f8f8f8')
histo_bt = tk.Button(listEffcts_fr, text="Histogram", command=hist, state=DISABLED, bg='white')
gauss_bt = tk.Button(listEffcts_fr, text="Gauss", command=gauss,  bg='white')
ray_bt = tk.Button(listEffcts_fr, text="Rayleigh", command=rayleigh,  bg='white')
median_bt = tk.Button(listEffcts_fr, text="Median", command=median,  bg='white')
sobel_bt = tk.Button(listEffcts_fr, text="Sobel", command=sobel,  bg='white')
otsu_bt = tk.Button(listEffcts_fr, text="Otsu", command=otsu,  bg='white')
kmeans_bt = tk.Button(listEffcts_fr, text="Kmeans", command=kmeans,  bg='white')
histo_bt.pack(side=TOP, padx=10, pady=4)
gauss_bt.pack(side=TOP, padx=10, pady=4)
ray_bt.pack(side=TOP, padx=10, pady=4)
median_bt.pack(side=TOP, padx=10, pady=4)
sobel_bt.pack(side=TOP, padx=10, pady=4)
otsu_bt.pack(side=TOP, padx=10, pady=4)
kmeans_bt.pack(side=TOP, padx=10, pady=4)
listEffcts_fr.place(relx=0.01, rely=0.125, anchor=NW)
Label(applyEfct_fr, text='Apply effects', bg='#f8f8f8').place(relx=.06, rely=0.125,anchor=W)
applyEfct_fr.pack(side=LEFT,padx=30, pady=20)
applyEfct_fr.configure(background='#f8f8f8')

topFrame_fr = Frame(root)
topFrame_fr.configure(background='#f8f8f8')
topFrame_fr.pack(padx=10, pady=3)

selectFile = Frame(topFrame_fr)
selectFile.configure(background='#f8f8f8')
selectFile.pack(padx=10, pady=5, side=TOP)

images_fr = Frame(topFrame_fr)
images_fr.configure(background='#f8f8f8')
images_fr.pack(padx=10, pady=10, side=TOP)


fileName = StringVar()
fileName.set("No se ha seleccionado ningún archivo")

selectFolder_bt = tk.Button(selectFile, text="Abrir archivo", command=loadFiles, bg='white')
selectFolder_bt.pack(pady=5)



files_fr =Frame(topFrame_fr)
files_fr.configure(background='#f8f8f8')
files_fr.pack(padx=10, pady=10, side=LEFT)





text_fr = Frame(topFrame_fr)
text_fr.configure(background='gray')


info_text = tk.Text(images_fr, width = 40, height = 11)
info_text.pack(padx=10, pady=5, side=BOTTOM)

f = Figure()
canvas = FigureCanvasTkAgg(f, master=images_fr)


root.mainloop()