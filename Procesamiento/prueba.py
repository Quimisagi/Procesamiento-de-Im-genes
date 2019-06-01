
import warnings
import pydicom

import os
import numpy as np
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def convu():

    kernel = [[1,1,1],
              [1,1,1],
              [1,1,1]]

    prueba = [[1,2,1,2, 3],
              [2,1,1,2, 3],
              [1,2,1,2,4],
              [4,3,2,1,2],
              [1,2,3,4, 5]]

    count = 0

    convu = np.zeros((3,3))

    for y in range(1,4):
        for x in range(1,4):
             count = 0
             for i in range(3):
                for j in range(3):
                    
                    count += kernel[i][j] * prueba[y+i-1][x+j-1]
             if count > 255:
                 count = 255
             if count < 0:
                 count = 0;
             print(y-1)
             convu[y-1][x-1] = count
    
    print(convu)

def gauss(N = 3, sigma = 1):

    
    #gaussMatY = []
    gaussian = np.zeros((N,N))
    total = 0

    for y in range (-1,2):
        #gaussMatX = []
        for x in range(-1,2):
            #We write down the left side (scalar value) of the equation: (1/ 2 * PI * SIGMA)
            left_side = (1 / (np.pi * 2 * np.power(sigma, 2)))
            #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
            aux = -(np.power(x, 2) + np.power(y, 2)) / (2 * np.power(sigma, 2))
            #Finally, we use the Exp function to calculate the e^(aux)
            right_side = np.exp(aux)
            total += left_side * right_side
            #gaussMatX.append(left_side * right_side)
            gaussian[y+1][x+1] = left_side * right_side
        #gaussMatY.append(gaussMatX)


    print(total)
    pan = 1/total
    
    for y in range (3):
        for x in range(3):
            gaussian[x][y] =  gaussian[x][y] * pan

    print(gaussian)

def rayleight(N = 3, sigma = 1):
    
    rayMat = np.zeros((N,N))
    total = 0

    for y in range (3):
        #gaussMatX = []
        for x in range(3):
            #We write down the left side (scalar value) of the equation: (1/ 2 * PI * SIGMA)
            left_side = np.divide(x * y, np.power(sigma, 4))

            #We start writting down the right side: -(X^2 + Y^2) / 2 * SIGMA^2
            aux = -(np.power(x, 2) + np.power(y, 2)) / (2 * np.power(sigma, 2))
            #Finally, we use the Exp function to calculate the e^(aux)
            right_side = np.exp(aux)
            total += left_side * right_side
            #gaussMatX.append(left_side * right_side)
            rayMat[y][x] = left_side * right_side
        #gaussMatY.append(gaussMatX)


    print(total)
    pan = 1/total
    
    for y in range (3):
        for x in range(3):
            rayMat[x][y] =  rayMat[x][y] * pan

    print(rayMat)


convu()
