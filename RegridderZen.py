import numpy as np
from scipy.sparse import csr_matrix

def RegridderZen(inArray, outSize):
    inSize = inArray.shape  # get the size of the inarray
    xFact = inSize[0] / outSize[0]  # y resize factor
    yFact = inSize[1] / outSize[1]  # x resize factor
    
    
    # ==========Y MATRIX CREATION -YES IT USES ALL X STUFF (maybe I named it wrong at the start lol?)=======
    xRange = np.arange(0, inSize[0] + xFact, xFact) #get a range of pixel slices to make
    xRangeAdj = xRange.copy()  # copy the xRange array to xRangeAdj
    

    for i in range(len(xRangeAdj)):
        if xRangeAdj[i] == np.fix(xRangeAdj[i]):  # catch true integers
            xRangeAdj[i] = xRangeAdj[i] - 0.1  # force any int values to int64(floor to the previous integer for counting reasons
            
    if inSize[0] != outSize[0]:
        yMatrix = np.zeros((outSize[0], inSize[0]))
        for i in range(1, outSize[0]+1):
            xPixelsInvolved = np.arange(np.fix(xRangeAdj[i]), np.fix(xRangeAdj[i + 1]) + 1)  # get the pixels involved
            xPixelPercent = np.ones_like(xPixelsInvolved)  # get the percentages involved
            
            if len(xPixelsInvolved) == 1:
                xPixelPercent[0] = xRange[i+1] - xRange[i]  
            else:
                if i == 0:
                    xPixelPercent[0] = np.ceil(xRange[i]) + 1 - xRange[i]  
                else:
                    xPixelPercent[0] = np.ceil(xRange[i]) - xRange[i]  

                xPixelPercent[-1] = xRange[i+1] - np.fix(xRangeAdj[i+1])
            
            for j in range(len(xPixelsInvolved)):
                yMatrix[i, xPixelsInvolved[j] + 1] = xPixelPercent[j]  

    else:
        yMatrix = np.diag(np.ones(inSize[0]))
        
    #========== X MATRIX CREATION - YES IT USES ALL Y STUFF (maybe I named it wrong at the start lol?) ==========#
    yRange = np.arange(0, inSize[1] + 1, yFact) #?
    yRangeAdj = np.copy(yRange)

    for i in range(len(yRangeAdj)):
        if yRangeAdj[i] == np.fix(yRangeAdj[i]):
            yRangeAdj[i] = yRangeAdj[i] - 0.1
        
    if inSize[1] != outSize[1]:
        xMatrix = np.zeros((outSize[1], inSize[1]))
        for i in range(1, outSize[1]+1):
            yPixelsInvolved = np.arange(np.fix(yRangeAdj[i]), np.fix(yRangeAdj[i + 1]) + 1)  # get the pixels involved
            yPixelPercent = np.ones_like(yPixelsInvolved)  # get the percentages involved
            
            if len(yPixelsInvolved) == 1:
                yPixelPercent[0] = yRange[i+1] - yRange[i]
            else:
                if i == 0:
                    yPixelPercent[0] = np.ceil(yRange[i]) + 1 - yRange[i]  
                else:
                    yPixelPercent[0] = np.ceil(yRange[i]) - yRange[i]  

                yPixelPercent[-1] = yRange[i+1] - np.fix(yRangeAdj[i+1])
            
            for j in range(len(yPixelsInvolved)):
                xMatrix[i, yPixelsInvolved[j] + 1] = yPixelPercent[j]  

    else:
        xMatrix = np.diag(np.ones(inSize[0]))#??
    
    # Convert to sparse matricies for speed (mostly 0's)
    if inSize[0]>1000 or inSize[1]>1000 or outSize[0]>1000 or outSize[1]>1000:
        # if array is large, use sparse arrays
        yMatrix = csr_matrix(yMatrix)
        xMatrix = csr_matrix(xMatrix)
        # ========== THE RESIZE ==========
        outArray = yMatrix*inArray*xMatrix
    else:
        outArray = yMatrix*inArray*xMatrix
    
    
    return outArray
