import math
from math import pow
import numpy as np

def main(data,dimensions):

	#Calculate average distance in each direction
	avgDistYZ = 0
	avgDistXZ = 0
	avgDistXY = 0

	for i in range(len(data)):
		avgDistYZ = avgDistYZ + round(math.sqrt(pow(abs(data[i]['y']),2) + pow(abs(data[i]['z']),2)),3)
		avgDistXZ = avgDistXZ + round(math.sqrt(pow(abs(data[i]['x']),2) + pow(abs(data[i]['z']),2)),3)
		avgDistXY = avgDistXY + round(math.sqrt(pow(abs(data[i]['x']),2) + pow(abs(data[i]['y']),2)),3)

	avgDistYZ /= len(data)
	avgDistXZ /= len(data)
	avgDistXY /= len(data)

	vesicleSuccessCounter = 0

	#Define ellipse based on average distances in each direction 
	for i in range(len(data)):
		equationOfVesicle = pow(data[i]['x'],2)/pow(avgDistYZ,2) + pow(data[i]['y'],2)/pow(avgDistXZ,2) + pow(data[i]['z'],2)/pow(avgDistXY,2)
		#Region of acceptable ellipse is defined between 0.5 and 2 
		if 0.5 < equationOfVesicle and equationOfVesicle < 2:
			vesicleSuccessCounter+=1
	print "Percentage of points that fit nanovesicle structure: ",vesicleSuccessCounter/float(len(data)) * 100, "%"

	#If the percentage of points that fit the ellipse is greater than or equal to 70%, then we consider it a vesicle
	if (vesicleSuccessCounter/float(len(data)) * 100) > 70:
		return 1
	else:
		return 0
