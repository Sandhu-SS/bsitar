  /////////////////////////////////////////////////////////////////////////
  // Function to calculate H matrix (as defined in GS_ns_getH_stan)
  /////////////////////////////////////////////////////////////////////////
  matrix GS_ns_getH_stan(vector knots, int normalize) {
  int Nintk = num_elements(knots) - 8;
  real C11 = 6 / ((knots[5] - knots[2]) * (knots[5] - knots[3]));
  real C31 = 6 / ((knots[6] - knots[3]) * (knots[5] - knots[3]));
  real C21 = -C11 - C31;
  real Cp22 = 6 / ((knots[Nintk + 6] - knots[Nintk + 3]) * (knots[Nintk + 6] - knots[Nintk + 4]));
  real Cp2 = 6 / ((knots[Nintk + 7] - knots[Nintk + 4]) * (knots[Nintk + 6] - knots[Nintk + 4]));
  real Cp12 = -Cp22 - Cp2;

  if (Nintk == 0) {
    matrix[2, 4] H = to_matrix(transpose([3, 0, 2, 1, 1, 2, 0, 3]), 2, 4);
    if (normalize) {
      vector[2] sumH = H * rep_vector(1.0, 4);
      for (i in 1:2) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else if (Nintk == 1) {
    matrix[3, 5] H = transpose(to_matrix([
      -C21 / C11, 1, 0, 0, 0, 
      0, -C31 / C21, 1, -Cp22 / Cp12, 0, 
      0, 0, 0, 1, -Cp12 / Cp2
    ], 5, 3));
    if (normalize) {
      vector[3] sumH = H * rep_vector(1.0, 5);
      for (i in 1:3) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else if (Nintk == 2) {
    matrix[4, 3] H1 = append_row(
      append_row(
        append_row(rep_matrix(1, 1, 3), [0, 1, -C21 / C31]),
        rep_matrix(0, 0, 3)
      ),
      rep_matrix(0, 2, 3)
    );
    matrix[4, 3] H3 = append_row(
      append_row(
        append_row(rep_matrix(0, 2, 3), rep_matrix(0, 0, 3)),
        [-Cp12 / Cp22, 1, 0]
      ),
      rep_matrix(1, 1, 3)
    );
    matrix[4, 6] H = append_col(H1, H3);
    if (normalize) {
      vector[4] sumH = H * rep_vector(1.0, 6);
      for (i in 1:4) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  } else {
    int K = Nintk - 2;
    int N = Nintk + 2;
    matrix[N, 3] H1 = append_row(
      append_row(
        append_row(rep_matrix(1, 1, 3), [0, 1, -C21 / C31]),
        rep_matrix(0, K, 3)
      ),
      rep_matrix(0, 2, 3)
    );
    matrix[N, K] H2 = append_row(
      append_row(
        rep_matrix(0, 2, K), diag_matrix(rep_vector(1, K))
      ),
      rep_matrix(0, 2, K)
    );
    matrix[N, 3] H3 = append_row(
      append_row(
        append_row(rep_matrix(0, 2, 3), rep_matrix(0, K, 3)),
        [-Cp12 / Cp22, 1, 0]
      ),
      rep_matrix(1, 1, 3)
    );
    matrix[N, N + 2] H = append_col(append_col(H1, H2), H3);
    if (normalize) {
      vector[N] sumH = H * rep_vector(1.0, N + 2);
      for (i in 1:N) if (sumH[i] != 0) H[i, ] /= sumH[i];
    }
    return transpose(H);
  }
} // end matrix GS_ns_getH_stan
