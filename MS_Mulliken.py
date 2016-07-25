import numpy as np
import matplotlib.pyplot as plt

filename = "E:\\temp\\For python\\Mulliken\\SO2.outmol"

# Step1: generate an array
MullikenFlag = False
num = 0   # TEST
listLikeArray = []
for line in open(filename):
    num += 1   # TEST
    #print(str(num))   # TEST
    if MullikenFlag == False:           # Start flag
        if "Population analysis for representation" in line:
            MullikenFlag = True
            continue
    else:
        try:
            #print(str(line))  # TEST
            assert '.' in line
            if line[3] != ' ':                                   # Avoid line break
                atomNum = int(line[0:4])
                atom_n = int(line[5])
                atom_l = int(line[6])
                atom_ms = int(line[7:9])
                occupation = float(line[9:16])

                info = line[16:].strip()
                splitting = info.split('.')
            else:
                splitting += line.strip().split('.')
                del listLikeArray[-1]                             # Delete the incomplete array unit

            while '' in splitting:
                splitting.remove('')

            splitting = map(lambda x: int(x), splitting)
            listLikeArray.append(splitting)

            #print(listLikeArray)
        except AssertionError:                                     # End of the Mulliken array
            break

#print(listLikeArray)
assert len(listLikeArray) == len(listLikeArray[-1])
MullikenPopulationArray = np.zeros([len(listLikeArray), len(listLikeArray[-1])])
for i in xrange(len(listLikeArray)):
    for j in xrange(len(listLikeArray[i])):
        MullikenPopulationArray[i,j] = listLikeArray[i][j]
MullikenPopulationArray /= 2000.0


# Step 2: plotting
#MullikenPopulationArray = MullikenPopulationArray[1:5, 1:5]
plt.imshow(MullikenPopulationArray, interpolation = "nearest", cmap = "bone")
plt.colorbar()
plt.show()