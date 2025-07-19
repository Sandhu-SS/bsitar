  /////////////////////////////////////////////////////////////////////////
  // Function to calculate truncated power basis based rcs matrix
  /////////////////////////////////////////////////////////////////////////
  /**
   * Calculates the Restricted Cubic Spline (RCS) basis matrix.
   *
   * This Stan function provides the core logic for generating a Restricted Cubic Spline
   * basis matrix, given an input vector `x` and pre-defined fullknots.
   *
   * IMPORTANT: Knot selection logic (e.g., quantiles, `Hmisc::rcspline.eval`)
   * must be performed OUTSIDE of Stan, in your data preparation script (R, Python, etc.).
   * The inclx is now defined within the functiin and not as an argument. 
   * The arguments are exact same as GS_nsp/nsk stan
   *
   * Args:
   * x: A vector of observations for which to calculate the basis.
   * fullknots: A vector of ordered unique fullknots. Must have at least 3 fullknots.
   * derivs: Integer representing the derivative order (0 for value, 1 for 1st derivs, 2 for 2nd derivs).
   * intercept: Integer (0 or 1). If 1, adds an intercept column.
   * inclx: Integer (0 or 1). If 1, includes the 'x' column (for derivs = 0).
   *
   * Returns:
   * A matrix where rows correspond to observations in `x` and columns
   * correspond to the basis functions.
   */
   
  matrix GS_rcs_call_stan(vector x, 
                        vector knotsx, 
                        vector bknotsx,
                        int intercept, 
                        int derivs, 
                        real centerval,
                        int normalize,
                        int preH) {
                        
    int inclx = 1;                  
    int N = num_elements(x);
    
    vector[num_elements(knotsx)+2] fullknots;
    vector[num_elements(fullknots)-2] knots;
    vector[2] bknots;
    fullknots = append_row(append_row(rep_vector(bknotsx[1], 1), knotsx), rep_vector(bknotsx[2], 1));
    knots = segment(fullknots, 2, num_elements(fullknots)-2);
    bknots = append_row(head(bknotsx, 1), tail(bknotsx, 1)); 
    
    int nk = num_elements(fullknots);
    
    if (nk < 3) {
      reject("rcs_matrix: Number of fullknots must be at least 3.");
    }
    
    matrix[N, nk - 1] basis_evals;
    real knot1 = fullknots[1];
    real knot_nk = fullknots[nk];
    real knot_nk_minus_1 = fullknots[nk - 1];
    real denom_common = (knot_nk - knot1)^2;
    real denom_spline_diff = (knot_nk - knot_nk_minus_1);

    matrix[N, nk] Xx_pos; 
    for (n_idx in 1:N) {
      for (k_idx in 1:nk) {
        Xx_pos[n_idx, k_idx] = fmax(0.0, x[n_idx] - fullknots[k_idx]);
      }
    }

    if (derivs == 0) {
      basis_evals[:, 1] = x;
    } else if (derivs == 1) {
      basis_evals[:, 1] = rep_vector(1.0, N);
    } else if (derivs == 2) {
      basis_evals[:, 1] = rep_vector(0.0, N);
    } else {
      reject("rcs_matrix: derivs must be 0, 1, or 2.");
    }

    for (j in 1:(nk - 2)) {
      int jp1 = j + 1; 
      if (derivs == 0) {
        basis_evals[:, jp1] =
          (pow(Xx_pos[:, j], 3) -
           pow(Xx_pos[:, nk - 1], 3) * (fullknots[nk] - fullknots[j]) / denom_spline_diff +
           pow(Xx_pos[:, nk], 3) * (fullknots[nk - 1] - fullknots[j]) / denom_spline_diff) / denom_common;
      } else if (derivs == 1) {
        basis_evals[:, jp1] =
          (3.0 * pow(Xx_pos[:, j], 2)) / denom_common -
          (3.0 * pow(Xx_pos[:, nk - 1], 2)) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (3.0 * pow(Xx_pos[:, nk], 2)) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      } else if (derivs == 2) {
        basis_evals[:, jp1] =
          (6.0 * Xx_pos[:, j]) / denom_common -
          (6.0 * Xx_pos[:, nk - 1]) * (fullknots[nk] - fullknots[j]) / (denom_spline_diff * denom_common) +
          (6.0 * Xx_pos[:, nk]) * (fullknots[nk - 1] - fullknots[j]) / (denom_spline_diff * denom_common);
      }
    }

    matrix[N, 0] final_basis_evals_empty = rep_matrix(0.0, N, 0); 

    if (intercept == 1) {
      if (derivs == 0) {
        if (inclx == 1) {
          matrix[N, nk] temp_basis;
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk)] = basis_evals; 
          return temp_basis;
        } else {
          matrix[N, nk-1] temp_basis; 
          temp_basis[:, 1] = rep_vector(1.0, N);
          temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)];
          return temp_basis;
        }
      } else if (derivs == 1) {
        if (inclx == 1) {
            matrix[N, nk] temp_basis;
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk)] = basis_evals; 
            return temp_basis;
        } else {
            matrix[N, nk-1] temp_basis; 
            temp_basis[:, 1] = rep_vector(0.0, N); 
            temp_basis[:, 2:(nk-1)] = basis_evals[:, 2:(nk-1)]; 
            return temp_basis;
        }
      } else if (derivs == 2) {
        if (inclx == 1) {
             matrix[N, nk] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             print(temp_basis);
             print(basis_evals);
             temp_basis[:, 3:(nk - 0)] = basis_evals[:, 2:(nk - 1)]; 
             return temp_basis;
        } else {
             matrix[N, nk -1] temp_basis; 
             temp_basis[:, 1] = rep_vector(0.0, N); 
             temp_basis[:, 2] = rep_vector(0.0, N); 
             temp_basis[:, 3:(nk-1)] = basis_evals[:, 3:(nk - 1)]; 
             return temp_basis;
        }
      }
    } else { // intercept == 0
      if (inclx == 0) {
        return basis_evals[:, 2:(nk-1)]; 
      } else {
        return basis_evals;
      }
    }
    // Should not reach here
    return final_basis_evals_empty;
  } // end matrix GS_rcs_call_stan
