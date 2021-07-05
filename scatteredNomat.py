import numpy as np
import scipy.io
from scipy import interpolate


def scatteredPython(tissue_ids,mySliceOnXY_model_spatial,xq,sliceLength,xqLength):
	
	print("Start")

	mySliceOnXY_model_spatial = np.reshape(np.transpose(np.array(mySliceOnXY_model_spatial)),[int(sliceLength),3])

	xq = np.reshape(np.transpose(np.array(xq)),[int(xqLength),3])

	tissue_ids = np.transpose(np.array(tissue_ids))
	
	print(mySliceOnXY_model_spatial)

	result = interpolate.griddata(mySliceOnXY_model_spatial,tissue_ids,xq,method='nearest')	
	return(result)	
	print('Success')



###### MATLAB CODE FOR FUNCTION ######


#cd ~/Downloads/
#load('mySliceOnXY_model_spatial.mat')
#load('tissue_ids.mat')
#load('xq.mat')
#mySliceOnXY_model_spatial = transpose(mySliceOnXY_model_spatial);
#xq = transpose(xq);
#tic
#result = py.scatteredNomat.scatteredTest(tissue_ids',mySliceOnXY_model_spatial(:).',xq(:).');
#result = double(py.array.array('d', py.numpy.nditer(result)));
#toc



