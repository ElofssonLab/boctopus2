import math, sys, os, string


listwhat = ["pssm"]
neighbor = 20

def myneighbourMain(protein, readfilename, neighborbaseHome):#{{{
    for what in listwhat:
        tempTong = []
        if what == "pssm":
            f     = open(neighborbaseHome + "tmp/" + readfilename, "r")
            lines = f.readlines()
            f.close()

        for i in range(0, len(lines)):
            lines[i] = lines[i].strip()
            temp = string.split(lines[i])
            tempEntry = []
            for j in range(0,len(temp)):
                tempEntry.append(temp[j])
            tempTong.append(tempEntry)

        ### From "tempTong" To "culled"
        unitL  = len(tempTong[0])-1
        culled = []
        for i in range(0,len(tempTong)):
            temp = []
            if i < neighbor:
                for j in range(0,i+neighbor+1):
                    if len(tempTong) > j:
                        if abs(int(tempTong[j][0])-int(tempTong[i][0])) <= neighbor:
                            temp.append(tempTong[j])
            elif len(tempTong)-i-1 < neighbor:
                for j in range(i-neighbor,len(tempTong)):
                    if abs(int(tempTong[j][0])-int(tempTong[i][0])) <= neighbor:
                        temp.append(tempTong[j])
            else:
                for j in range(i-neighbor,i+neighbor+1):
                    if abs(int(tempTong[j][0])-int(tempTong[i][0])) <= neighbor:
                        temp.append(tempTong[j])

            tempEntry = []
            for j in range(int(tempTong[i][0])-neighbor, int(tempTong[i][0])+neighbor+1):
                presence = 0
                for k in range(0,len(temp)):
                    if int(temp[k][0]) == j:
                        tempEntry += temp[k][1:]
                        presence = 1
                        break
                if presence == 0:
                    for k in range(1,len(tempTong[0])):
                        tempEntry.append("NA")

            Len = len(tempEntry)
            oTempEntry  = [tempTong[i][0]]
            oTempEntry +=  tempEntry[(Len/unitL/2)*unitL:(Len/unitL/2)*unitL+unitL]          ### Len = 60, unitL = 20   tempEntry[20:40]
            for j in range(1,(Len/unitL/2)+1):                                               ### j = 1
                oTempEntry += tempEntry[(Len/unitL/2-j)*unitL:(Len/unitL/2-j)*unitL+unitL]   ### tempEntry[0:20]
                oTempEntry += tempEntry[(Len/unitL/2+j)*unitL:(Len/unitL/2+j)*unitL+unitL]  ### tempEntry[40:60]
            culled.append(oTempEntry)


        ### Recording the culled data
        f = open(neighborbaseHome + "pssm/" + protein + "_" + str(2*neighbor+1) + ".filtered", "w")
        for i in range(0,len(culled)):
            for j in range(0,len(culled[i])):
                f.write(culled[i][j]+" ")
            f.write("\n")
        f.close()

    return
#}}}
