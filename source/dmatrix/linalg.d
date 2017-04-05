module dmatrix.linalg;

import dmatrix.init;

DMatrix matrixTranspose(DMatrix input) {
  double[] output = new double[input.shape[0]*input.shape[1]];
  int index = 0;
  for(int i=0; i< input.shape[1]; i++) {
    for(int j=0; j< input.shape[0]; j++) {
      output[index] = input.elements[j*input.shape[1]+i];
      index++;
    }
  }
  int[] resshape = [input.shape[1],input.shape[0]];
  return DMatrix(resshape,output);
}
