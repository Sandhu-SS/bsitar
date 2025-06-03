

// Efficient vector matching/counting
int num_matches(vector x, real y, real z) {
  int n = 0;
  int N = num_elements(x);
  if (z == 0) {
    for (i in 1:N) if (x[i] == y) n += 1;
  } else if (z == 1.0) {
    for (i in 1:N) if (x[i] > y) n += 1;
  } else if (z == -1.0) {
    for (i in 1:N) if (x[i] < y) n += 1;
  }
  return n;
}

array[] int which_equal(vector x, real y, real z) {
  int N = num_elements(x);
  int n = num_matches(x, y, z);
  array[n] int match_positions;
  int pos = 1;
  if (z == 0) {
    for (i in 1:N) if (x[i] == y) { match_positions[pos] = i; pos += 1; }
  } else if (z == 1.0) {
    for (i in 1:N) if (x[i] > y) { match_positions[pos] = i; pos += 1; }
  } else if (z == -1.0) {
    for (i in 1:N) if (x[i] < y) { match_positions[pos] = i; pos += 1; }
  }
  return match_positions;
}

// B-spline basis and derivative (single function, efficient)
matrix GS_bs_stan(vector x, vector knots, vector bknots, 
                  vector fullknots, vector allknots, 
                  int N, int degree, int ord, 
                  int Nintk, int Nk, int Nintervals,
                  int intercept, int calcderiv) {
  matrix[N, Nintervals] M1 = rep_matrix(0, N, Nintervals);
  for (i in 1:N)
    for (j in 1:Nintervals)
      M1[i, j] = (allknots[j] <= x[i] && x[i] < allknots[j + 1]) ? 1 : 0;

  array[num_matches(x, bknots[2], 0)] int lastknot_index = which_equal(x, bknots[2], 0);
  if (size(lastknot_index) > 0)
    for (ii in 1:size(lastknot_index))
      for (jj in (Nintervals - degree):(Nintervals))
        M1[lastknot_index[ii], jj] = 1;


  matrix[N, Nintervals] M2x;
  matrix[N, Nintervals] M1x = M1;
  for (p in 1:degree) {
    M2x = rep_matrix(0, N, Nintervals);
    for (i in 1:(Nintervals - p)) {
      vector[N] C1 = (allknots[i + p] == allknots[i]) ? rep_vector(0.0, N) : (x - allknots[i]) / (allknots[i + p] - allknots[i]);
      vector[N] C2 = (allknots[i + p + 1] == allknots[i + 1]) ? rep_vector(0.0, N) : (allknots[i + p + 1] - x) / (allknots[i + p + 1] - allknots[i + 1]);
      for (n in 1:N)
        M2x[n, i] = C1[n] * M1x[n, i] + C2[n] * M1x[n, i + 1];
    }
    if (p != degree) M1x = M2x;
  }
  matrix[N, Nintervals - degree] M2 = M2x[, 1:(Nintervals - degree)];
  matrix[N, Nintervals - degree + 1] M1xx = M1x[, 1:(Nintervals - degree + 1)];
  if (calcderiv) {
    matrix[N, Nintervals - degree] deriv;
    for (i in 1:(Nintervals - degree)) {
      real dC1 = (allknots[i + degree] == allknots[i]) ? 0 : degree / (allknots[i + degree] - allknots[i]);
      real dC2 = (allknots[i + degree + 1] == allknots[i + 1]) ? 0 : degree / (allknots[i + degree + 1] - allknots[i + 1]);
      for (n in 1:N)
        deriv[n, i] = dC1 * M1xx[n, i] - dC2 * M1xx[n, i + 1];
    }
    return deriv;
  }
  return M2;
}

// Tuple version for both basis and derivative
tuple(matrix, matrix) GS_bs_stan_tuple(vector x, vector knots, vector bknots, 
                  vector fullknots, vector allknots, 
                  int N, int degree, int ord, 
                  int Nintk, int Nk, int Nintervals,
                  int intercept, int calcderiv) {
  matrix[N, Nintervals] M1 = rep_matrix(0, N, Nintervals);
  for (i in 1:N)
    for (j in 1:Nintervals)
      M1[i, j] = (allknots[j] <= x[i] && x[i] < allknots[j + 1]) ? 1 : 0;

  array[num_matches(x, bknots[2], 0)] int lastknot_index = which_equal(x, bknots[2], 0);
  if (size(lastknot_index) > 0)
    for (ii in 1:size(lastknot_index))
      for (jj in (Nintervals - degree):(Nintervals))
        M1[lastknot_index[ii], jj] = 1;

  matrix[N, Nintervals] M2x;
  matrix[N, Nintervals] M1x = M1;
  for (p in 1:degree) {
    for (i in 1:(Nintervals - p)) {
      vector[N] C1 = (allknots[i + p] == allknots[i]) ? rep_vector(0.0, N) : (x - allknots[i]) / (allknots[i + p] - allknots[i]);
      vector[N] C2 = (allknots[i + p + 1] == allknots[i + 1]) ? rep_vector(0.0, N) : (allknots[i + p + 1] - x) / (allknots[i + p + 1] - allknots[i + 1]);
      M2x[, i] = C1 .* M1x[, i] + C2 .* M1x[, i + 1];
    }
    if (p != degree) M1x = M2x;
  }
  matrix[N, Nintervals - degree] M2 = M2x[, 1:(Nintervals - degree)];
  matrix[N, Nintervals - degree + 1] M1xx = M1x[, 1:(Nintervals - degree + 1)];
  matrix[N, Nintervals - degree] deriv;
  for (i in 1:(Nintervals - degree)) {
    real dC1 = (allknots[i + degree] == allknots[i]) ? 0 : degree / (allknots[i + degree] - allknots[i]);
    real dC2 = (allknots[i + degree + 1] == allknots[i + 1]) ? 0 : degree / (allknots[i + degree + 1] - allknots[i + 1]);
    deriv[, i] = dC1 * M1xx[, i] - dC2 * M1xx[, i + 1];
  }
  return (M2, deriv);
}

// Main B-spline function
matrix GS_ns_stan(
  vector x, vector knots, vector bknots, 
  vector fullknots, vector allknots, 
  int N, int degree, int ord, 
  int Nintk, int Nk, int Nintervals,
  int intercept, int calcderiv, int normalize,
  int preH,
  matrix MatpreH
) {
  int n_basis = Nintervals - degree;
  matrix[N, n_basis] bs = GS_bs_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
  matrix[N, n_basis] bsderiv;
  if (calcderiv) {
    bsderiv = GS_bs_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, 1);
  }

  // Precompute boundary indices only if needed
  int n_below = num_matches(x, bknots[1], -1);
  int n_above = num_matches(x, bknots[2], 1);
  array[n_below] int x_below_boundary;
  array[n_above] int x_above_boundary;
  if (n_below > 0) x_below_boundary = which_equal(x, bknots[1], -1);
  if (n_above > 0) x_above_boundary = which_equal(x, bknots[2], 1);

  if (n_below > 0 || n_above > 0) {
    tuple(matrix[2, n_basis], matrix[2, n_basis]) my_tuple =
      GS_bs_stan_tuple(bknots, knots, bknots, fullknots, allknots, 2, degree, ord, Nintk, Nk, Nintervals, intercept, 0);
    matrix[2, n_basis] bs_bknots = my_tuple.1;
    matrix[2, n_basis] bsderiv_bknots = my_tuple.2;

    if (n_below > 0) {
      vector[n_below] x_below = x[x_below_boundary];
      row_vector[n_basis] bsknot1 = bs_bknots[1, ];
      row_vector[n_basis] bsknot1_deriv = bsderiv_bknots[1, ];
      matrix[n_below, n_basis] rep_bsknot1 = rep_matrix(bsknot1, n_below);
      matrix[n_below, n_basis] rep_bsknot1_deriv = rep_matrix(bsknot1_deriv, n_below);
      matrix[n_below, n_basis] xdiff = rep_matrix(to_row_vector(x_below - bknots[1]), n_basis)';
      bs[x_below_boundary, ] = rep_bsknot1 + rep_bsknot1_deriv .* xdiff;
      if (calcderiv) {
        for (i in 1:n_below)
          bsderiv[x_below_boundary[i], ] = bsknot1_deriv;
      }
    }
    if (n_above > 0) {
      vector[n_above] x_above = x[x_above_boundary];
      row_vector[n_basis] bsknot2 = bs_bknots[2, ];
      row_vector[n_basis] bsknot2_deriv = bsderiv_bknots[2, ];
      matrix[n_above, n_basis] rep_bsknot2 = rep_matrix(bsknot2, n_above);
      matrix[n_above, n_basis] rep_bsknot2_deriv = rep_matrix(bsknot2_deriv, n_above);
      matrix[n_above, n_basis] xdiff = rep_matrix(to_row_vector(x_above - bknots[2]), n_basis)';
      bs[x_above_boundary, ] = rep_bsknot2 + rep_bsknot2_deriv .* xdiff;
      if (calcderiv) {
        for (i in 1:n_above)
          bsderiv[x_above_boundary[i], ] = bsknot2_deriv;
      }
    }
  }

  // Precompute H matrix only once
  matrix[Nk + 2, Nk] H;
  if (preH) {
    H = MatpreH;
  } else {
    H = GS_ns_getH_stan(allknots, normalize);
  }

  int ncolselect = Nk + intercept - 1;
  matrix[N, Nk] out = calcderiv ? (bsderiv * H) : (bs * H);

  // Avoid unnecessary allocation
  if (intercept) {
    return out;
  } else {
    return out[, 2:cols(out)];
  }
}


// Call GS_nsp_call_stan
matrix GS_nsp_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int intercept, int derivs, real centerval, int normalize,
                        int preH, matrix MatpreH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  int degree = 3;
  int ord = degree + 1;
  int Nk = Nintk + 2;
  int df = Nintk + 1 + intercept;
  int calcderiv = (derivs > 1) ? 0 : derivs;
  int ncolselect = Nk + intercept - 1;
  vector[Nintk + 2 * ord] allknots = append_row(append_row(rep_vector(bknots[1], ord), knots), rep_vector(bknots[2], ord));
  int Nintervals = num_elements(allknots) - 1;
  matrix[N, ncolselect] out = GS_ns_stan(x, knots, bknots, fullknots, allknots, N, degree, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH, MatpreH);

  if (centerval != 0) {
    matrix[1, df] cenout = GS_ns_stan(centerval + rep_vector(0.0, 1), knots, bknots, fullknots, allknots, 1, degree, ord, Nintk, Nk, Nintervals, intercept, calcderiv, normalize, preH, MatpreH);
    if (!calcderiv) {
      if (intercept) {
        for (i in 2:cols(cenout)) out[, i] -= cenout[1, i];
      } else {
        for (i in 1:cols(cenout)) out[, i] -= cenout[1, i];
      }
    }
  }
  return out;
}


// Call GS_nsk_call_stan
matrix GS_nsk_call_stan(vector x, vector knotsx, vector bknotsx, 
                        int intercept, int derivs, real centerval, int normalize,
                        int preH, matrix MatpreH) {
  int N = num_elements(x);
  vector[num_elements(knotsx) + 2] fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
  vector[num_elements(fullknots) - 2] knots = segment(fullknots, 2, num_elements(fullknots) - 2);
  vector[2] bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1));
  int Nintk = num_elements(knots);
  int degree = 3;
  int ord = degree + 1;
  int Nk = Nintk + 2;
  int ncolselect = Nk + intercept - 1;
  matrix[N, ncolselect] basis = GS_nsp_call_stan(x, knots, bknots, intercept, derivs, centerval, normalize, preH, MatpreH);
  matrix[Nk, ncolselect] kbasis = GS_nsp_call_stan(fullknots, knots, bknots, intercept, derivs, centerval, normalize, preH, MatpreH);
  matrix[N, ncolselect] out;
  if (intercept) {
    out = basis * inverse(kbasis);
  } else {
    matrix[N, ncolselect + 1] kout = append_col(rep_vector(1, N), basis) * inverse(append_col(rep_vector(1, Nk), kbasis));
    out = kout[, 2:cols(kout)];
  }
  return out;
}
