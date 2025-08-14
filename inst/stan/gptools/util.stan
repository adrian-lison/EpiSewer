// Scalars -----------------------------------------------------------------------------------------

/**
Assert that two integers are equal.

:param actual: Actual value.
:param desired: Desired value.
*/
void assert_equal(int actual, int desired) {
    if (actual != desired) {
        reject(actual, " is not equal to ", desired);
    }
}


/**
Check whether two values are close.

The actual value :code:`x` and desired value :code:`y` may differ by at most
:code:`tol = rtol * y + atol`, where :code:`rtol` is the relative tolerance, and :code:`atol` is the
absolute tolerance. The tolerance :code:`tol` is clipped below at :math:`10^{-15}` to avoid
rejection due to rounding errors.

:param actual: Actual value :code:`x`.
:param desired: Desired value :code:`y`.
:param rtol: Relative tolerance :code:`r`.
:param atol: Absolute tolerance :code:`a`.
:returns: :code:`1` if the values are close, :code:`0` otherwise.
*/
int is_close(real actual, real desired, real rtol, real atol) {
    // We always allow a tolerance of at least 1e-15 in case there are rounding errors.
    real tol = fmax(atol + rtol * abs(desired), 1e-15);
    if (abs(actual - desired) <= tol) {
        return 1;
    }
    return 0;
}

/**
Assert that two values are close. See :stan:func:`is_close(real, real, real, real)` for a
description of the parameters.
*/
void assert_close(real actual, real desired, real rtol, real atol) {
    if (!is_close(actual, desired, rtol, atol)) {
        reject(actual, " is not close to ", desired);
    }
}

/**
Assert that two values are close. See :stan:func:`is_close(real, real, real, real)` for a
description of the parameters.
*/
void assert_close(real actual, real desired) {
    assert_close(actual, desired, 1e-6, 0);
}

/**
Check whether a possibly complex value is finite.

:param x: Value to check.
:returns: :code:`1` if the value is finite, :code:`0` otherwise.
*/
int is_finite(complex x) {
    real rx = get_real(x);
    real ix = get_imag(x);
    if (is_nan(rx) || is_nan(ix) || is_inf(rx) || is_inf(ix)) {
        return 0;
    }
    return 1;
}


// Vectors -----------------------------------------------------------------------------------------

/**
Assert that two vectors are close. See :stan:func:`is_close(real, real, real, real)` for a
description of the parameters.
*/
void assert_close(vector actual, vector desired, real rtol, real atol) {
    int n = size(desired);
    int m = size(actual);
    if (m != n) {
        reject("number of elements are not equal: size(desired)=", n, "; size(actual)=", m);
    }
    for (i in 1:size(actual)) {
        if (!is_close(actual[i], desired[i], rtol, atol)) {
            reject(actual[i], " is not close to ", desired[i], " at position ", i);
        }
    }
}

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(vector actual, vector desired) {
    assert_close(actual, desired, 1e-6, 0);
}

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(vector actual, real desired, real rtol, real atol) {
    assert_close(actual, rep_vector(desired, size(actual)), rtol, atol);
}

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(vector actual, real desired) {
    assert_close(actual, desired, 1e-6, 0);
}

/**
Check whether all elements of a vector are finite.

:param x: Vector to check.
:returns: :code:`1` if all elements of the vector are finite, :code:`0` otherwise.
*/
int is_finite(vector x) {
    for (i in 1:size(x)) {
        if(!is_finite(x[i])) {
            return 0;
        }
    }
    return 1;
}

/**
Assert that all elements of a vector are finite.

:param: Vector to check.
*/
void assert_finite(vector x) {
    int n = size(x);
    for (i in 1:n) {
        if (!is_finite(x[i])) {
            reject(x[i], " at index ", i, " is not finite");
        }
    }
}

// Matrices ----------------------------------------------------------------------------------------

/**
Check whether all elements of a matrix are finite.

:param: Vector to check.
*/
int is_finite(matrix x) {
    for (i in 1:rows(x)) {
        for (j in 1:cols(x)) {
            if (!is_finite(x[i, j])) {
                return 0;
            }
        }
    }
    return 1;
}


/**
Pretty-print a matrix.
*/
void print_matrix(complex_matrix x) {
    print("matrix with ", rows(x), " rows and ", cols(x), " columns");
    for (i in 1:rows(x)) {
        print(x[i]);
    }
}

/**
Assert that two matrices are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(matrix actual, matrix desired, real rtol, real atol) {
    array [2] int nshape = dims(desired);
    array [2] int mshape = dims(actual);
    if (mshape[1] != nshape[1]) {
        reject("number of rows are not equal: dims(desired)[1]=", nshape[1], "; size(actual)[1]=",
               mshape[1]);
    }
    if (mshape[2] != nshape[2]) {
        reject("number of columns are not equal: dims(desired)[2]=", nshape[2],
               "; size(actual)[2]=", mshape[2]);
    }
    for (i in 1:nshape[1]) {
        for (j in 1:nshape[2]) {
            if (!is_close(actual[i, j], desired[i, j], rtol, atol)) {
                reject(actual[i, j], " is not close to ", desired[i, j], " at row ", i, ", column ",
                       j);
            }
        }
    }
}

/**
Assert that two matrices are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(matrix actual, matrix desired) {
    assert_close(actual, desired, 1e-6, 0);
}

/**
Assert that two matrices are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(matrix actual, real desired, real rtol, real atol) {
    array [2] int shape = dims(actual);
    assert_close(actual, rep_matrix(desired, shape[1], shape[2]), rtol, atol);
}

/**
Assert that two matrices are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(matrix actual, real desired) {
    assert_close(actual, desired, 1e-6, 0);
}

// Complex vectors ---------------------------------------------------------------------------------

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(complex_vector actual, complex_vector desired, real rtol, real atol) {
    int n = size(desired);
    int m = size(actual);
    if (m != n) {
        reject("number of elements are not equal: size(desired)=", n, "; size(actual)=", m);
    }
    assert_close(get_real(actual), get_real(desired), rtol, atol);
    assert_close(get_imag(actual), get_imag(desired), rtol, atol);
}


/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(complex_vector actual, complex_vector desired) {
    assert_close(actual, desired, 1e-6, 0);
}

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(complex_vector actual, complex desired, real rtol, real atol) {
    assert_close(actual, rep_vector(desired, size(actual)), rtol, atol);
}

/**
Assert that two vectors are close. See :stan:func:`is_close` for description of parameters.
*/
void assert_close(complex_vector actual, complex desired) {
    assert_close(actual, desired, 1e-6, 0);
}


// Real Fourier transforms -------------------------------------------------------------------------

/**
Compute the one-dimensional discrete Fourier transform for real input.

:param y: Real signal with :code:`n` elements to transform.
:returns: Truncated vector of Fourier coefficients with :code:`n %/% 2 + 1` elements.
*/
complex_vector rfft(vector y) {
    return fft(y)[:size(y) %/% 2 + 1];
}

/**
Expand truncated one-dimensional discrete Fourier transform coefficients for real input to full
Fourier coefficients.
*/
complex_vector expand_rfft(complex_vector y, int n) {
    complex_vector[n] result;
    int nrfft = n %/% 2 + 1;
    if (size(y) != nrfft) {
        reject("expected complex vector with ", nrfft, " elements but got ", size(y));
    }
    int ncomplex = (n - 1) %/% 2;
    result[:nrfft] = y;
    result[nrfft + 1:n] = conj(reverse(y[2:1 + ncomplex]));
    return result;
}

/**
Compute the one-dimensional inverse discrete Fourier transform for real output.

:param z: Truncated vector of Fourier coefficents with :code:`n %/% 2 + 1` elements.
:param n: Length of the signal (required because the length of the signal cannot be determined from
    :code:`z` alone).
:returns: Real signal with :code:`n` elements.
*/
vector inv_rfft(complex_vector z, int n) {
    return get_real(inv_fft(expand_rfft(z, n)));
}

/**
Compute the two-dimensional discrete Fourier transform for real input.

:param y: Real signal with :code:`n` rows and :code:`m` columns to transform.
:returns: Truncated vector of Fourier coefficients with :code:`n` rows and :code:`m %/% 2 + 1`
    elements.
*/
complex_matrix rfft2(matrix y) {
    return fft2(y)[:, :cols(y) %/% 2 + 1];
}

/**
Compute the two-dimensional inverse discrete Fourier transform for real output.

:param z: Truncated vector of Fourier coefficients with :code:`n` rows and :code:`m %/% 2 + 1`
    elements.
:param m: Number of columns of the signal (required because the number of columns cannot be
    determined from :code:`z` alone).
:returns: Real signal with :code:`n` rows and :code:`m` columns.
*/
matrix inv_rfft2(complex_matrix z, int m) {
    int n = rows(z);
    complex_matrix[n, m] x;
    int mrfft = m %/% 2 + 1;
    int mcomplex = (m - 1) %/% 2;
    x[:, 1:mrfft] = z[:, 1:mrfft];
    // Fill redundant values.
    for (i in 1:n) {
        x[i, mrfft + 1:m] = conj(reverse(z[i, 2:1 + mcomplex]));
    }
    // Reverse the order to account for negative frequencies.
    for (i in mrfft + 1:mrfft + mcomplex) {
        x[2:, i] = reverse(x[2:, i]);
    }
    return get_real(inv_fft2(x));
}

// Conditional location and scale parameters for multivariate normal distributions -----------------

/**
Evaluate the conditional location and scale parameter of a univariate normal random variable given
correlated observations from a multivariate normal distribution.

:param y: Observation to condition on.
:param cov11: Marginal variance of the target random variable.
:param cov21: Covariance between :code:`y` and the target random variable.
:param cov22: Covariance amongst the elements of :code:`y`.
:returns: Location and scale as a vector.
*/
vector gp_conditional_loc_scale(vector y, real cov11, vector cov21, matrix cov22) {
    if (size(y) == 0) {
        return [0, sqrt(cov11)]';
    }
    vector[size(y)] v = mdivide_left_spd(cov22, cov21);
    return [dot_product(v, y), sqrt(cov11 - dot_product(v, cov21))]';
}
