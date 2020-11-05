import os.path
import sys
from os import path
fileName = input("Ingrese el nombre del archivo a cargar: ")
if(path.isfile(fileName)):
    test = open("lafalopita.txt", "r")
    iMax = 0
    iAvg = 0
    cont = 0
    for i in test:
        if ("iterations" in i):
            cont += 1
            i = i.split()
            tempVal = int(i[0])
            if(tempVal > iMax):
                iMax = tempVal
                iAvg += tempVal
                iAvg /= (cont + 1)
                print("iMax = {} iAvg = {}".format(str(iMax), str(iAvg)))
else:
    print("No se encontro un archivo con ese nombre")
