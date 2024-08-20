import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt

def compute_ar_coefficients(X, order):
    """
    Computes the AR coefficients using the Yule-Walker equations.

    Parameters:
    X (array-like): The time series data.
    order (int): The order of the AR model.

    Returns:
    ar_coefficients (numpy array): The AR coefficients.
    noise_variance (float): The variance of the white noise.
    """
    # Convert X to a numpy array if it isn't already
    X = np.asarray(X)

    # Compute the autocovariance function up to the given order
    autocovariances = [np.cov(X[:-lag], X[lag:])[0, 1] for lag in range(order + 1)]
    
    # Construct the autocovariance matrix (Toeplitz matrix)
    R = np.array([[autocovariances[abs(i - j)] for j in range(order)] for i in range(order)])
    
    # Construct the right-hand side vector
    r = np.array(autocovariances[1:order+1])

    # Solve the Yule-Walker equations to find AR coefficients
    ar_coefficients = np.linalg.solve(R, r)

    # Compute the noise variance
    noise_variance = autocovariances[0] - np.dot(ar_coefficients, r)
    
    return ar_coefficients, noise_variance
    
X = pfits.getdata('data2/disturbance.fits').squeeze()

# Sample autocovariances
gamma_0 = np.var(X)
gamma_1 = np.cov(X[:-1], X[1:])[0, 1]
gamma_2 = np.cov(X[:-2], X[2:])[0, 1]

# Matrix of autocovariances
R = np.array([[gamma_0, gamma_1],
              [gamma_1, gamma_0]])

# Vector of autocovariances
r = np.array([gamma_1, gamma_2])

# Compute the coefficients phi_1 and phi_2
phi = np.linalg.inv(R).dot(r)

phi_1, phi_2 = phi[0], phi[1]

# Compute the variance of the white noise
sigma2_epsilon = gamma_0 - phi_1 * gamma_1 - phi_2 * gamma_2

X_past = X[1:-1]
X_past_past = X[:-2]
X_pred = phi_1*X_past+phi_2*X_past_past

plt.figure()
plt.plot(X)
plt.plot(X[2:])
plt.plot(X_pred)
plt.show()
plt.legend(("X","X_shift","X_pred"))

print(f"phi_1: {phi_1}, phi_2: {phi_2}, variance: {sigma2_epsilon}")