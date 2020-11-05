iMax = 0
iAvg = 0
test = open("lafalopita.txt", "r")
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
