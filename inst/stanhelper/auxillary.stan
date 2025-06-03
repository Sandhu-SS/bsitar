int VecSame(vector A, vector B) {
  int N = num_elements(A);
  for (i in 1:N) {
    if (A[i] != B[i]) return 0;
  }
  return 1;
}

int MatSame(matrix A, matrix B) {
  int R = rows(A);
  int C = cols(A);
  for (i in 1:R) {
    for (j in 1:C) {
      if (A[i, j] != B[i, j]) return 0;
    }
  }
  return 1;
}

vector repeat_vector(vector input, int K) {
  int N = num_elements(input);
  vector[N*K] repvec;
  for (k in 1:K)
    for (i in 1:N)
      repvec[i + (k - 1) * N] = input[i];
  return repvec;
}

matrix repeat_matrix(matrix input, int K) {
  int N = rows(input);
  int M = cols(input);
  matrix[N*K, M] repmat;
  for (k in 1:K)
    for (i in 1:N)
      repmat[i + (k - 1) * N] = input[i];
  return repmat;
}
