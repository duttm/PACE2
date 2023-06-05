import math, sys
import numpy as np
import random
import WaterTrapping
import Vesicle
import GroConverter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
	#Convert gro file to .dat peptide + solvent files and read data from peptide file
    file = 'peptideFile.dat'
    GroConverter.main(sys.argv[1])
    peptideFile = open(file, 'r')
    peptideData = peptideFile.readlines()
    data = readPeptideData(peptideData)


    #Remove outliers
    newData = data
    print 'Number of data points after dropping outliers:',len(newData) 

    #Reassign identifierIDs
    for i in range(len(newData)):
        newData[i]['identifierID'] = i+1

    #Grab box dimensions
    box_dimensions = getBoxDimensions()

    #Plot system before shifting
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.scatter3D(newData[:]['x'], newData[:]['y'], newData[:]['z'], cmap='hsv')
    fig.savefig('beforeShift.png')

    #Temporarily rename data for shifting
    shiftingData = newData

    #Perform all of shifting 2 times
    for i in range(2):

        #Grab random peptide particle
        randomPeptideIndex = random.randint(0,len(shiftingData)-1)

        xShiftFactor = shiftingData[randomPeptideIndex]['x'] - box_dimensions[0]/2
        yShiftFactor = shiftingData[randomPeptideIndex]['y'] - box_dimensions[1]/2
        zShiftFactor = shiftingData[randomPeptideIndex]['z'] - box_dimensions[2]/2

        #Shift first based on random peptide particle
        data_firstShift = shift(shiftingData,xShiftFactor,yShiftFactor,zShiftFactor,box_dimensions)

        fig = plt.figure()
        ax = plt.axes(projection="3d")

        ax.scatter3D(data_firstShift[:]['x'], data_firstShift[:]['y'], data_firstShift[:]['z'], cmap='hsv')
        fig.savefig('after1stshift.png')

        #Temp reassign data to shift 40 more times
        shiftingData = data_firstShift

        #Find center of mass
        x_mean = np.mean(shiftingData[:]['x'])
        y_mean = np.mean(shiftingData[:]['y'])
        z_mean = np.mean(shiftingData[:]['z'])

        #Shift
        shiftingData = shift(shiftingData,x_mean - box_dimensions[0]/2,y_mean - box_dimensions[1]/2,z_mean - box_dimensions[2]/2,box_dimensions)

        fig = plt.figure()
        ax = plt.axes(projection="3d")

        ax.scatter3D(shiftingData[:]['x'], shiftingData[:]['y'], shiftingData[:]['z'], cmap='hsv')
        fig.savefig('after2ndShift.png')

        #Iteratively shift 39 more times based on center of mass
        for i in range(40):
    	    #Find center of mass
    	    x_mean = np.mean(shiftingData[:]['x'])
    	    y_mean = np.mean(shiftingData[:]['y'])
    	    z_mean = np.mean(shiftingData[:]['z'])

    	    #Shift
    	    shiftingData = shift(shiftingData,x_mean + box_dimensions[0]/2,y_mean + box_dimensions[1]/2,z_mean + box_dimensions[2]/2,box_dimensions)

    #Final post-shifting data
    finalPeptideData = shiftingData

    #Plot final post-shifting data
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.scatter3D(shiftingData[:]['x'], shiftingData[:]['y'], shiftingData[:]['z'], cmap='hsv')
    fig.savefig('after40thShift.png')

    #Read data from solvent file
    completeData = readSolventData(finalPeptideData)

    #Check if structure traps water
    traps = WaterTrapping.main(completeData,box_dimensions)

    #Check vesicle characterization
    if traps:
        if Vesicle.main(newData,box_dimensions):
            file = open("output.dat","w")
            file.write("v")
        else:
            file = open("output.dat","w")
            file.write("d")
    else:
            file = open("output.dat","w")
            file.write("d")



## **HELPER FUNCTIONS**

#Read peptide data (only backbone), returned as structured numpy array
def readPeptideData(masterData):
    #Find pattern of backbones
    patternOfBackbones = 0
    for i in range(1,len(masterData)):
        curLine = masterData[i].split()
        patternOfBackbones+=1

        if curLine[1] == "BB":
            break

    #Populate structured numpy array using masterData
    dt = np.dtype([('residueID', np.unicode_, 16), ('moleculeID', np.unicode_, 16), ('identifierID', int), ('x', np.float64), ('y', np.float64), ('z', np.float64)])
    data = np.empty(len(masterData)/patternOfBackbones,dtype = dt)

    #Only record data that is backbone
    index = 0
    for i in range(len(masterData)):
        if (i%patternOfBackbones == 0):
            curLine = masterData[i].split()
            data[index]['residueID'] = curLine[0]
            data[index]['moleculeID'] = curLine[1]
            data[index]['identifierID'] = curLine[2]
            data[index]['x'] = round(float(curLine[3]),3)
            data[index]['y'] = round(float(curLine[4]),3)
            data[index]['z'] = round(float(curLine[5]),3)
            index+=1

    return data

#Read box dimensions
def getBoxDimensions():
    #Read solvent txt data as second argument
    file2 = "solventFile.dat"
    solventFile = open(file2,'r')
    solventData = solventFile.readlines()

    #Extract box dimensions from end of file
    temp = solventData[len(solventData)-1].split()

    box_dimensions = []
    box_dimensions.append(round(float(temp[0]),3))
    box_dimensions.append(round(float(temp[1]),3))
    box_dimensions.append(round(float(temp[2]),3))

    return box_dimensions

#Read solvent data, building array with both peptide and solvent data and returning box dimensions 
def readSolventData(peptideData):
    #Read solvent txt data as second argument
    file2 = "solventFile.dat"
    solventFile = open(file2,'r')
    solventData = solventFile.readlines()

    solventData = np.delete(solventData,len(solventData)-1)

    dt = np.dtype([('residueID', np.unicode_, 16), ('moleculeID', np.unicode_, 16), ('identifierID', int), ('x', np.float64), ('y', np.float64), ('z', np.float64)])
    completeData = np.empty(len(peptideData) + len(solventData),dtype = dt)

    for i in range(len(peptideData)):
        completeData[i] = peptideData[i]

    for i in range(len(solventData)):
        curLine = solventData[i].split()
        completeData[i+len(peptideData)]['residueID'] = 'solvent'
        completeData[i+len(peptideData)]['moleculeID'] = 'S'
        completeData[i+len(peptideData)]['identifierID'] = i+len(peptideData)

        if (len(curLine) == 6):
            completeData[i+len(peptideData)]['x'] = curLine[3]
            completeData[i+len(peptideData)]['y'] = curLine[4]
            completeData[i+len(peptideData)]['z'] = curLine[5]
        else:
            completeData[i+len(peptideData)]['x'] = curLine[2]
            completeData[i+len(peptideData)]['y'] = curLine[3]
            completeData[i+len(peptideData)]['z'] = curLine[4]

    return completeData


#Screen out outliers
def dropOutliers(data):
    #Calculate Q1 and Q3 in each dimension
    x_q1 = np.percentile(data[:]['x'],25)
    x_q3 = np.percentile(data[:]['x'],75)
    y_q1 = np.percentile(data[:]['y'],25)
    y_q3 = np.percentile(data[:]['y'],75)
    z_q1 = np.percentile(data[:]['z'],25)
    z_q3 = np.percentile(data[:]['z'],75)

    #Calculate IQR in each dimension
    x_IQR = x_q3 - x_q1
    y_IQR = y_q3 - y_q1
    z_IQR = z_q3 - z_q1

    #Calculate allowed tolerance in each dimension
    x_max = x_q3 + 1.5*x_IQR
    x_min = x_q1 - 1.5*x_IQR
    y_max = y_q3 + 1.5*y_IQR
    y_min = y_q1 - 1.5*y_IQR
    z_max = z_q3 + 1.5*z_IQR
    z_min = z_q1 - 1.5*z_IQR

    #Populate list of indices that follow out of allowed tolerance in any dimension, and delete those indices
    outlierIndices = []
    for i in range(len(data)):
        if ((data[i]['x'] < x_min) | (data[i]['x'] > x_max) | (data[i]['y'] < y_min) | (data[i]['y'] > y_max) | (data[i]['z'] < z_min) | (data[i]['z'] > z_max)):
            outlierIndices.append(i)

    print('Indices where there are outliers:')
    for i in range(len(outlierIndices)):
        print(outlierIndices[i])

    newData = np.delete(data,outlierIndices,0)
    return newData

#Shifting
def shift(data,xShiftFactor,yShiftFactor,zShiftFactor,box_dimensions):
    for i in range(len(data)):
        data[i]['x'] = round(data[i]['x'] - xShiftFactor,3)
        data[i]['y'] = round(data[i]['y'] - yShiftFactor,3)
        data[i]['z'] = round(data[i]['z'] - zShiftFactor,3)

        #Account for periodic boundary conditions
        if (data[i]['x'] > box_dimensions[0]):
            data[i]['x'] = round(data[i]['x'] - box_dimensions[0],3)
        elif (data[i]['x'] < 0):
            data[i]['x'] = data[i]['x'] + box_dimensions[0]

        if (data[i]['y'] > box_dimensions[1]):
            data[i]['y'] = round(data[i]['y'] - box_dimensions[1],3)
        elif (data[i]['y'] < 0):
            data[i]['y'] = round(data[i]['y'] + box_dimensions[1],3)

        if (data[i]['z'] > box_dimensions[2]):
            data[i]['z'] = round(data[i]['z'] - box_dimensions[2],3)
        elif (data[i]['z'] < 0):
            data[i]['z'] = round(data[i]['z'] + box_dimensions[2],3)

    return data
    
if __name__ == "__main__":
    main()