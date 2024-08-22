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
    autocovariances = [np.var(X)]  # Start with the variance at lag 0
    autocovariances += [np.cov(X[:-lag], X[lag:])[0, 1] for lag in range(1, order + 1)]
    
    # Construct the autocovariance matrix (Toeplitz matrix)
    R = np.array([[autocovariances[abs(i - j)] for j in range(order)] for i in range(order)])
    
    # Construct the right-hand side vector
    r = np.array(autocovariances[1:order+1])

    # Solve the Yule-Walker equations to find AR coefficients
    ar_coefficients = np.linalg.solve(R, r)

    # Compute the noise variance
    noise_variance = autocovariances[0] - np.dot(ar_coefficients, r)
    
    return ar_coefficients, noise_variance
    
def forecast_two_steps_ahead(X, ar_coefficients):
    order = len(ar_coefficients)
    
    # Ensure we have enough data to make predictions
    if len(X) < order:
        raise ValueError(f"The length of the input data must be at least equal to the order of the AR model ({order}).")
    
    # Get the last `order` values from the time series
    recent_values = X[-order:]
    
    # One-step ahead forecast
    X_t_plus_1 = np.dot(ar_coefficients, recent_values[::-1])

    # Add the one-step forecast to the recent values for the two-step forecast
    recent_values = np.append(X_t_plus_1, recent_values[:-1])

    # Two-step ahead forecast
    X_t_plus_2 = np.dot(ar_coefficients, recent_values[::-1])

    return X_t_plus_1, X_t_plus_2

def forecast_future_values(X, ar_coefficients, steps):
    """
    Forecast future values using the AR model.

    Parameters:
    X (array-like): The time series data.
    ar_coefficients (numpy array): The AR coefficients.
    steps (int): The number of future steps to forecast.

    Returns:
    forecasted_values (numpy array): The forecasted future values.
    """
    order = len(ar_coefficients)
    forecasted_values = []

    # Use the last 'order' values to start forecasting
    recent_values = list(X[-order:])

    for _ in range(steps):
        # Forecast the next value
        next_value = np.dot(ar_coefficients, recent_values[::-1])
        forecasted_values.append(next_value)

        # Update the recent values list
        recent_values = [next_value] + recent_values[:-1]

    return np.array(forecasted_values)

def predict_two_steps_ahead_series(X, ar_coefficients):
    """
    Predict two steps ahead for the whole time series.

    Parameters:
    X (array-like): The time series data.
    ar_coefficients (numpy array): The AR coefficients.

    Returns:
    predictions (numpy array): Array of two-step ahead predictions.
    """
    order = len(ar_coefficients)
    n = len(X)
    
    # Initialize an array to hold the predictions
    predictions = np.full(n, np.nan)  # Start with NaNs for unavailable predictions

    for i in range(order, n - 1):
        # Get the last 'order' values needed to make predictions
        recent_values = X[i-order:i]
        
        # One-step ahead prediction
        one_step_ahead = np.dot(ar_coefficients, recent_values[::-1])
        
        # Add one-step prediction to recent values to predict two steps ahead
        recent_values = [one_step_ahead] + list(recent_values[:-1])
        
        # Two-step ahead prediction
        two_steps_ahead = np.dot(ar_coefficients, recent_values[::-1])
        
        # Store the prediction
        predictions[i + 1] = two_steps_ahead
    
    return predictions

X = pfits.getdata('data2/disturbance.fits').squeeze()

t = np.arange(0,1000)
X = np.sin(t*0.1)
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
AR_order = 4
X_past = X[1:-1]
X_past_past = X[:-2]
X_pred = phi_1*X_past+phi_2*X_past_past
X_pred = np.zeros(X.shape[0]-AR_order)
ar_coefficients = np.array([phi_1,phi_2])


ar_coefficients,_ = compute_ar_coefficients(X, AR_order)

for i in range(AR_order):
        X_pred += ar_coefficients[i]*X[AR_order-i-1:-i-1]

X_t_plus_2 = predict_two_steps_ahead_series(X, ar_coefficients)

plt.figure()
plt.plot(X[:50])
plt.plot(X[2:52])
# plt.plot(X_pred[:50])
plt.plot(X_t_plus_2[2:52])
plt.legend(("X","X_shift","X_pred"))
plt.show()


print(f"phi_1: {phi_1}, phi_2: {phi_2}, variance: {sigma2_epsilon}")