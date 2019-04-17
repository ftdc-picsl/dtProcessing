# Load this first, so it doesn't mask ANTsR origin function
library(RNifti)
library(ANTsR)


## This works as long as the voxel space is consistent, ie voxel (i,j,k) is the same place in NIFTI and ITK


warpTracts <- function(seedFile, refFile, tractFile, outputTractFile, transforms, invert) {

  counter = 1

  ## Length of buffer - 8 bytes per value, so 1024*1024*3*10 == 240 Mb
  ## The amount of RAM used will be several times this, as we have to store various intermediate point arrays
  ##
  ## Buffer length is a multiple of 3 so it can be reshaped into 3D point matrix 
  bufferLength = 1024*1024*3*4

  ## Main buffer holds points to be transformed
  ##
  ## Each call to transformPoints will read transforms from disk; we want this to be minimized, hence the buffering
  ##
  buffer = vector("numeric", bufferLength)

  ## estimate lengths for numPoints and seedIndex buffers
  numPointsBuffer = vector("numeric", bufferLength / 100)
  seedIndexBuffer = vector("numeric", bufferLength / 100)

  ## Next available index to place data
  bufferIndex = 1;
  bufferNumTractIndex = 1;

  imNifti = readNifti(seedFile)

  imITK = antsImageRead(seedFile)

  refITK = antsImageRead(refFile)

  refNifti = readNifti(refFile)

  tractSource = file(tractFile, open = "rb", raw = T)

  outputStream = file(outputTractFile, open = "wb", raw = T)

  numPoints = readBin(tractSource, "numeric", 1, endian = "big", size = 4)

  while (length(numPoints) > 0) {

    if (bufferIndex + 3 * numPoints > bufferLength) {
      ## Time to write some output
      warpedPoints = transformTractPoints(buffer, (bufferIndex - 1) / 3, imNifti, imITK, refNifti, refITK, transforms, invert)
      writeTracts(bufferNumTractIndex - 1, numPointsBuffer, seedIndexBuffer, warpedPoints, outputStream)
      rm(warpedPoints)
      gc()
      bufferIndex = 1
      bufferNumTractIndex = 1
    }
    
    numPointsBuffer[bufferNumTractIndex] = numPoints
    
    seedIndexBuffer[bufferNumTractIndex] = readBin(tractSource, "numeric", 1, endian = "big", size = 4)
   
    buffer[bufferIndex:(bufferIndex - 1 + 3 * numPoints)] = readBin(tractSource, "numeric", numPoints * 3, endian = "big", size = 4)

    bufferIndex = bufferIndex + 3 * numPoints
    
    bufferNumTractIndex = bufferNumTractIndex + 1

    if (bufferNumTractIndex > length(numPointsBuffer)) {
      numPointsBuffer = c(numPointsBuffer, rep(0, length(numPointsBuffer)))
      seedIndexBuffer = c(seedIndexBuffer, rep(0, length(seedIndexBuffer)))
    }
    
    numPoints = readBin(tractSource, "numeric", 1, endian = "big", size = 4)
    
  }

  ## Write final output and close
  warpedPoints = transformTractPoints(buffer, (bufferIndex - 1) / 3, imNifti, imITK, refNifti, refITK, transforms, invert)
  writeTracts(bufferNumTractIndex - 1, numPointsBuffer, seedIndexBuffer, warpedPoints, outputStream)

  rm(warpedPoints)
  gc()
  
  close(tractSource)
  close(outputStream)
 
}

## Warp a bunch of points, returns a matrix in nifti space
transformTractPoints <- function(buffer, numPoints, imNifti, imITK, refNifti, refITK, transforms, invert) {

  ## Buffer will in general be slightly longer than 3 * numPoints
  pointsMat = matrix(buffer, ncol = 3, nrow = length(buffer) / 3, byrow = T)

  pointsMat = pointsMat[1:numPoints,]

  voxels = worldToVoxel(pointsMat, imNifti)

  itkPointsInput = antsTransformIndexToPhysicalPoint(imITK, voxels)

  ## Now warp the points (in ITK physical space)
  itkPointsWarped = as.matrix(antsApplyTransformsToPoints(3, itkPointsInput, transformlist = transforms, whichtoinvert = invert))

  ## Clean up as we go to free memory
  rm(pointsMat, voxels, itkPointsInput)
  gc()
  
  ## Now go back to NIFTI space via the ITK voxel space
  itkVoxelsOutput = antsTransformPhysicalPointToIndex(refITK, itkPointsWarped)
  
  niftiPointsOutputMat = voxelToWorld(itkVoxelsOutput, refNifti)

  rm(itkPointsWarped, itkVoxelsOutput)

  gc()
 

  return(niftiPointsOutputMat)
  
}

## Takes points in matrix format used by RNIfti / ANTsR
## Will transpose then write as array
writeTracts <- function(numTracts, tractNumPoints, seedIndices, points, outputStream) {

  points = t(points)

  dim(points) = c()

  pointIndex = 1

  for (t in 1:numTracts) {
    
    writeBin(tractNumPoints[t], outputStream, endian = "big", size = 4)

    writeBin(seedIndices[t], outputStream, endian = "big", size = 4)
    
    writeBin(points[pointIndex:(pointIndex - 1 + tractNumPoints[t] * 3)], outputStream, endian = "big", size = 4)

    pointIndex = pointIndex + 3 * tractNumPoints[t]
    
  }

  ## Check we've written all the points
  if (pointIndex != length(points) + 1) {
    stop("Did not write all available points")
  }
}


## Takes a mask and warps, and moves seed points to a destination space. Mostly useful to take seeds from T1 to DWI,
## but can go the other way for testing / diagnostic purposes.
## 
warpSeeds <- function(inputSeedImageFile, outputSeedImageFile, transforms, invert) {

  inputSeedImageITK = antsImageRead(inputSeedImageFile)

  outputSeedImageITK = antsImageRead(outputSeedImageFile)

  outputSeedImageNifti = readNifti(outputSeedImageFile)

  maskVoxels = which(as.array(inputSeedImageITK) > 0, arr.ind = T)

  maskPoints = antsTransformIndexToPhysicalPoint(inputSeedImageITK, maskVoxels)
  
  ## Now warp the points (in ITK physical space)
  pointsWarped = as.matrix(antsApplyTransformsToPoints(3, maskPoints, transformlist = transforms, whichtoinvert = invert))

  ## Conserve memory
  rm(maskPoints,maskVoxels)
  gc()
  
  ## Now go back to NIFTI space via the ITK voxel space
  voxelsOutput = antsTransformPhysicalPointToIndex(outputSeedImageITK, pointsWarped)
  
  niftiPointsOutput = voxelToWorld(voxelsOutput, outputSeedImageNifti) 

  rm(pointsWarped)
  gc()
  
  seedDensityImage = antsImageClone(outputSeedImageITK) * 0

  seedDensityArray = as.array(seedDensityImage)
  
  for (i in 1:nrow(voxelsOutput)) {
    seedDensityArray[voxelsOutput[i,1],voxelsOutput[i,2],voxelsOutput[i,3]] = seedDensityArray[voxelsOutput[i,1],voxelsOutput[i,2],voxelsOutput[i,3]] + 1
  }

  seedDensityImage[seedDensityArray > 0] = seedDensityArray[seedDensityArray > 0]
  
  return(list(niftiPointsOutput, seedDensityImage))

}
