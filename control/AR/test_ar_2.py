import numpy as np
import matplotlib.pyplot as plt


def compute_ar_coefficients(X, order):
    X = np.asarray(X)
    autocovariances = [np.var(X)]  # Start with the variance at lag 0
    autocovariances += [np.cov(X[:-lag], X[lag:])[0, 1] for lag in range(1, order + 1)]

    R = np.array([[autocovariances[abs(i - j)] for j in range(order)] for i in range(order)])
    r = np.array(autocovariances[1:order + 1])

    ar_coefficients = np.linalg.solve(R, r)
    noise_variance = autocovariances[0] - np.dot(ar_coefficients, r)

    return ar_coefficients, noise_variance


def predict_two_steps_ahead_series(X, ar_coefficients):
    order = len(ar_coefficients)
    n = len(X)

    predictions = np.full(n, np.nan)  # Initialize with NaNs

    for i in range(order, n - 2):
        recent_values = X[i - order:i]  # Get the last 'order' values needed to make predictions
        one_step_ahead = np.dot(ar_coefficients, recent_values[::-1])
        recent_values = np.insert(recent_values[:-1], 0, one_step_ahead)  # Update recent values
        two_steps_ahead = np.dot(ar_coefficients, recent_values[::-1])
        predictions[i + 2] = two_steps_ahead  # Store the two-step-ahead prediction

    return predictions


# Example usage
X = np.array([5, 6, 7, 8, 10, 12, 11, 13, 15, 14, 16, 18, 20, 19, 21, 22], dtype=float)  # Example time series
t = np.arange(0,50)
X = np.sin(t*0.1)
# Specify the order of the AR model
order = 1

# Compute the AR(5) coefficients
ar_coefficients, noise_variance = compute_ar_coefficients(X, order)

# Predict two steps ahead for the whole series
predictions = predict_two_steps_ahead_series(X, ar_coefficients)

# Shift the original X array by 2 steps to align with the predictions
X_shifted = X[2 + order:]  # Shift by 2 plus the order of the AR model

# Align predictions by removing the initial NaN values
predictions_aligned = predictions[2 + order:]

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(X_shifted, label='Original Series (shifted by 2 steps + order)', linestyle='--', marker='o')
plt.plot(predictions_aligned, label='Two-Step Ahead Predictions', marker='x', linestyle='dashed', color='red')
plt.title('Comparison of Two-Step Ahead Predictions and Original Series')
plt.xlabel('Time Step')
plt.ylabel('Value')
plt.legend()
plt.show()
