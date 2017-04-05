module dmatrix.init;

struct DMatrix{
  int[] shape;
  double[] elements;
  bool init = false;

  this(int[] s, double[] e) {
    shape = s;
    elements = e;
    init = true;
  }

  double acc(int row, int col) {
    return this.elements[row*this.shape[1]+col];
  }
}
