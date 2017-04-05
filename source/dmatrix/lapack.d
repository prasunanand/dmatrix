module dmatrix.lapack;

import dmatrix.init;

extern (C) {
  int LAPACKE_dgetrf (int matrix_layout, int m, int n, double* a, int lda, int* ipiv);
  int LAPACKE_dsyevr (int matrix_layout, char jobz, char range, char uplo, int n, double* a, int lda, double vl, double vu, int il, int iu, double abstol, int* m, double* w, double* z, int ldz, int* isuppz);
  int LAPACKE_dgetri (int matrix_layout, int n, double* a, int lda, const(int)* ipiv);
}

struct eighTuple{
  DMatrix kva, kve;

  this(DMatrix kva, DMatrix kve) {
    this.kva = kva;
    this.kve = kve;
  }
}

eighTuple eigh(DMatrix input) {
  double[] z = new double[input.shape[0] * input.shape[1]]; //eigenvalues
  double[] w = new double[input.shape[0]];  // eigenvectors
  double[] elements = input.elements.dup;

  double wi;
  int n = input.shape[0];
  double vu, vl;
  int[] m = new int[input.shape[0]];
  int[] isuppz = new int[2*input.shape[0]];
  int il = 1;
  int iu = input.shape[1];
  int ldz = n;
  double abstol = 0.001; //default value for abstol

  LAPACKE_dsyevr(101, 'V', 'A', 'L', n, elements.ptr, n, vl, vu, il, iu, abstol, m.ptr, w.ptr, z.ptr, ldz, isuppz.ptr);

  DMatrix kva = DMatrix([input.shape[0],1], w);
  DMatrix kve = DMatrix(input.shape, z);
  for(int zq = 0 ; zq < cast(int)kva.elements.length; zq++){
    if(kva.elements[zq]< 1e-6){
      kva.elements[zq] = 0;
    }
  }
  return eighTuple(kva, kve);
}



double det(DMatrix input) {
  double[] narr = input.elements.dup;
  int[] pivot = getrf(narr, input.shape.dup);

  int num_perm = 0;
  int j = 0;
  foreach(swap; pivot) {
    if(swap-1 != j) {num_perm += 1;}
    j++;
  }
  double prod;
  if(num_perm % 2 == 1) {
    prod = 1;
  }else{
   prod = -1; //# odd permutations => negative
  }
  int min = input.shape[0];
  if(input.shape[0]> input.shape[1]) {min = input.shape[1];}
  for(int i =0;i< min; i++) {
    prod *= narr[input.shape[0]*i + i];
  }
  return prod;
}

int[] getrf(double[] arr, int[] shape) {
  auto ipiv = new int[shape[0]+1];
  LAPACKE_dgetrf(101, shape[0],shape[0],arr.ptr,shape[0],ipiv.ptr);
  return ipiv;
}

DMatrix inverse(DMatrix input) {
  double[] elements= input.elements.dup;
  int LWORK = input.shape[0]*input.shape[0];
  double[] WORK = new double[input.shape[0]*input.shape[0]];
  auto ipiv = new int[input.shape[0]+1];
  auto result = new double[input.shape[0]*input.shape[1]];
  int info;
  int output = LAPACKE_dgetrf(101, input.shape[0],input.shape[0],elements.ptr,input.shape[0],ipiv.ptr);
  int[] resshape = [input.shape[0],input.shape[0]];
  LAPACKE_dgetri(101, input.shape[0],elements.ptr, input.shape[0], ipiv.ptr);
  return DMatrix(input.shape, elements);
}
