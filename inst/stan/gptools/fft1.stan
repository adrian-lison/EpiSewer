// IMPORTANT: stan uses the questionable R indexing which is one-based and inclusive on both ends.
// I.e., x[1:3] includes x[1] to x[3]. More generally, x[i:j] comprises j - i + 1 elements. It could
// at least have been exclusive on the right...

/**
Evaluate the scale of Fourier coefficients.

:param cov_rfft: Precomputed real fast Fourier transform of the kernel with shape
    :code:`(..., n %/% 2 + 1)`.
:param n: Size of the real signal. Necessary because the size cannot be inferred from
    :code:`cov_rfft`.
:returns: Scale of Fourier coefficients with shape :code:`(..., n %/% 2 + 1)`.
*/
vector gp_evaluate_rfft_scale(vector cov_rfft, int n) {
    int nrfft = n %/% 2 + 1;
    vector[nrfft] result = n * cov_rfft / 2;
    // Check positive-definiteness.
    real minval = min(result);
    if (minval < 0) {
        reject("covariance matrix is not positive-definite (minimum eigenvalue is ", minval, ")");
    }
    // The first element has larger scale because it only has a real part but must still have the
    // right variance. The same applies to the last element if the number of elements is even
    // (Nyqvist frequency).
    result[1] *= 2;
    if (n % 2 == 0) {
        result[nrfft] *= 2;
    }
    return sqrt(result);
}


/**
Unpack the Fourier coefficients of a real Fourier transform with :code:`n %/% 2 + 1` elements to a
vector of :code:`n` elements.

:param z: Real Fourier transform coefficients.
:param n: Size of the real signal. Necessary because the size cannot be inferred from :code:`rfft`.

:returns: Unpacked vector of :code:`n` elements comprising the :code:`n %/% 2 + 1` real parts of the
    zero frequency term, complex terms, and Nyqvist frequency term (for even :code:`n`). The
    subsequent :code:`(n - 1) %/% 2` elements are the imaginary parts of complex coefficients.
*/
vector gp_unpack_rfft(complex_vector x, int n) {
    vector[n] z;
    int ncomplex = (n - 1) %/% 2;
    int nrfft = n %/% 2 + 1;
    z[1:nrfft] = get_real(x);
    z[1 + nrfft:n] = get_imag(x[2:1 + ncomplex]);
    return z;
}


/**
Transform a Gaussian process realization to white noise in the Fourier domain.

:param y: Realization of the Gaussian process with shape :code:`(..., n)`.
:param loc: Mean of the Gaussian process with shape :code:`(..., n)`.
:param cov_rfft: Precomputed real fast Fourier transform of the kernel with shape
    :code:`(..., n // 2 + 1)`.
:returns: Fourier-domain white noise with shape :code:`(..., n)`. See :stan:func:`gp_unpack_rfft`
    for details on the data structure.
*/
vector gp_rfft(vector y, vector loc, vector cov_rfft) {
    return gp_unpack_rfft(rfft(y - loc) ./ gp_evaluate_rfft_scale(cov_rfft, size(y)), size(y));
}


/**
Evaluate the log absolute determinant of the Jacobian associated with
:stan:func:`gp_rfft`.

:param cov_rfft: Precomputed real fast Fourier transform of the kernel with shape
    :code:`(..., n %/% 2 + 1)`.
:param n: Size of the real signal. Necessary because the size cannot be inferred from :code:`rfft`.
:returns: Log absolute determinant of the Jacobian.
*/
real gp_rfft_log_abs_det_jac(vector cov_rfft, int n) {
    vector[n %/% 2 + 1] rfft_scale = gp_evaluate_rfft_scale(cov_rfft, n);
    return - sum(log(rfft_scale[1:n %/% 2 + 1])) -sum(log(rfft_scale[2:(n + 1) %/% 2]))
        - log(2) * ((n - 1) %/% 2) + n * log(n) / 2;
}
real gp_rfft_log_abs_det_jacobian(vector cov_rfft, int n) {
    reject("`gp_rfft_log_abs_det_jacobian` has been renamed to `gp_rfft_log_abs_det_jac` to comply with upcoming changes (see https://github.com/stan-dev/stanc3/issues/1470 for details).");
}


/**
Evaluate the log probability of a one-dimensional Gaussian process realization in Fourier space.

:param y: Realization of a Gaussian process with :code:`n` grid points.
:param loc: Mean of the Gaussian process of size :code:`n`.
:param cov_rfft: Precomputed real fast Fourier transform of the kernel of size :code:`n %/% 2 + 1`.
:returns: Log probability of the Gaussian process realization.
*/
real gp_rfft_lpdf(vector y, vector loc, vector cov_rfft) {
    int n = size(y);
    int nrfft = n %/% 2 + 1;
    vector[n] z = gp_rfft(y, loc, cov_rfft);
    return std_normal_lpdf(z) + gp_rfft_log_abs_det_jac(cov_rfft, n);
}


/**
Transform a real vector with :code:`n` elements to a vector of complex Fourier coefficients with
:code:`n` elements ready for inverse real fast Fourier transformation.
*/
complex_vector gp_pack_rfft(vector z) {
    int n = size(z);  // Number of observations.
    int ncomplex = (n - 1) %/% 2;  // Number of complex Fourier coefficients.
    int nrfft = n %/% 2 + 1;  // Number of elements in the real FFT.
    int neg_offset = (n + 1) %/% 2;  // Offset at which the negative frequencies start.
    // Zero frequency, real part of positive frequency coefficients, and Nyqvist frequency.
    complex_vector[nrfft] fft = z[1:nrfft];
    // Imaginary part of positive frequency coefficients.
    fft[2:ncomplex + 1] += 1.0i * z[nrfft + 1:n];
    return fft;
}


/**
Transform white noise in the Fourier domain to a Gaussian process realization, i.e., a
*non-centered* parameterization in the Fourier domain.

The :code:`n` real white noise variables must be assembled into a length-:code:`n %/% 2 + 1` complex
vector with structure expected by the fast Fourier transform. The input vector :code:`z` comprises

- the real zero frequency term,
- :code:`(n - 1) %/% 2` real parts of positive frequency terms,
- the real Nyqvist frequency term if :code:`n` is even,
- and :code:`(n - 1) %/% 2` imaginary parts of positive frequency terms.

:param z: Fourier-domain white noise comprising :code:`n` elements.
:param loc: Mean of the Gaussian process.
:param cov_rfft: Real fast Fourier transform of the covariance kernel.

:returns: Realization of the Gaussian process with :code:`n` elements.
*/
vector gp_inv_rfft(vector z, vector loc, vector cov_rfft) {
    int n = size(z);
    vector[n %/% 2 + 1] rfft_scale = gp_evaluate_rfft_scale(cov_rfft, n);
    return get_real(inv_rfft(rfft_scale .* gp_pack_rfft(z), n)) + loc;
}

/**
Evaluate the real fast Fourier transform of the periodic squared exponential kernel.

:param n: Number of grid points.
:param sigma: Scale of the covariance.
:param length_scale: Correlation length.
:param period: Period for circular boundary conditions.
:returns: Fourier transform of the squared exponential kernel of size :code:`n %/% 2 + 1`.
*/
vector gp_periodic_exp_quad_cov_rfft(int n, real sigma, real length_scale, real period) {
    int nrfft = n %/% 2 + 1;
    return n * sigma ^ 2 * length_scale / period * sqrt(2 * pi())
        * exp(-2 * (pi() * linspaced_vector(nrfft, 0, nrfft - 1) * length_scale / period) ^ 2);
}

/**
Evaluate the real fast Fourier transform of the periodic Matern kernel.

:param nu: Smoothness parameter (1 / 2, 3 / 2, and 5 / 2 are typical values).
:param n: Number of grid points.
:param sigma: Scale of the covariance.
:param length_scale: Correlation length.
:param period: Period for circular boundary conditions.
:returns: Fourier transform of the squared exponential kernel of size :code:`n %/% 2 + 1`.
*/
vector gp_periodic_matern_cov_rfft(real nu, int n, real sigma, real length_scale, real period) {
    int nrfft = n %/% 2 + 1;
    vector[nrfft] k = linspaced_vector(nrfft, 0, nrfft - 1);
    return sigma ^ 2 * n * sqrt(2 * pi() / nu) * tgamma(nu + 0.5) / tgamma(nu)
        * (1 + 2 / nu * (pi() * length_scale / period * k) ^ 2) ^ -(nu + 0.5) * length_scale
        / period;
}
