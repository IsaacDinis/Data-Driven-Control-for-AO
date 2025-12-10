import numpy as np
from scipy.linalg import solve_discrete_lyapunov, solve_lyapunov, eig, expm, inv, sqrtm
import astropy.io.fits as pfits
from scipy.fft import fft
import matplotlib.pyplot as plt
import control as ct

def ajuste_sigmav(Aphi, Gammaphi, Cphi, datavar, option=1):
    """
    Adjust the covariance matrix of the noise driving the identified discrete model.

    Parameters:
        Aphi, Gammaphi, Cphi: State model of the form:
                              xphi(k+1) = Aphi * xphi(k) + Gammaphi * v(k)
                              phi(k) = Cphi * xphi(k)
                              or
                              d(xphi(t)) = Aphi * xphi(t) dt + Gammaphi * dv(t)
                              phi(t) = Cphi * xphi(t)
        datavar: Variance/covariance matrix of phi or trajectory used for identification
        option: 1 for discrete model (default), 2 for continuous model

    Returns:
        Sigmav: Variance of the noise driving the model
    """

    # Define nphi as the number of rows in Cphi
    nphi = Cphi.shape[0]

    # Calculate Sigmaphi
    if datavar.shape[1] == nphi:
        Sigmaphi = datavar
    else:
        ntraj = datavar.shape[1]
        datavar = datavar - np.mean(datavar, axis=1, keepdims=True)
        Sigmaphi = np.zeros((nphi, nphi))
        for j in range(ntraj):
            Sigmaphi += np.outer(datavar[:, j], datavar[:, j])
        Sigmaphi /= ntraj

    # Adjustment using Lyapunov
    Sigmavtestmat = []
    MM = []
    toto = np.zeros((nphi, nphi))

    for i in range(nphi):
        for j in range(i + 1):
            Sigmavtestk = toto.copy()
            Sigmavtestk[i, i] = 1
            if j < i:
                Sigmavtestk[j, j] = 1
                Sigmavtestk[i, j] = 1
                Sigmavtestk[j, i] = 1
            Sigmavtestmat.append(Sigmavtestk.flatten())
            if option == 1:
                Sigmaxtestk = solve_discrete_lyapunov(Aphi, Gammaphi @ Sigmavtestk @ Gammaphi.T)
            else:
                Sigmaxtestk = solve_lyapunov(Aphi, Gammaphi @ Sigmavtestk @ Gammaphi.T)
            varphitestk = Cphi @ Sigmaxtestk @ Cphi.T
            MM.append(varphitestk.flatten())

    Sigmavtestmat = np.array(Sigmavtestmat).T
    MM = np.array(MM).T
    coefsSigmav = np.linalg.lstsq(MM, Sigmaphi.flatten(), rcond=None)[0]
    SigmavL = Sigmavtestmat @ coefsSigmav
    Sigmav = SigmavL.reshape((nphi, nphi))

    # Check positivity
    minvp = np.min(np.abs(eig(Sigmav)[0]))
    if minvp < 0:
        print("** Warning: Sigmav has negative eigenvalues **")

    return Sigmav


def rediscret(A_c, gamma_c, C_c, var, fs):
    """
    Rediscretize continuous-time shaping filter.

    Parameters:
        A_c, gamma_c, C_c: Continuous-time state-space matrices.
        var: Variance of output.
        fs: Sampling frequency (Hertz).

    Returns:
        A_d, gamma_d, C_d: Discrete-time state-space matrices.
        rac_sigma_x: Square root of steady-state variance matrix for filter initialization.
    """

    # Discretization of state matrices
    A_d = expm(A_c / fs)

    XX = inv(A_c) @ (A_d - np.eye(A_d.shape[0]))
    C_d = C_c @ XX * fs
    gamma_d0 = XX @ gamma_c

    # Adjust discrete filter gain
    sigmav = ajuste_sigmav(A_d, gamma_d0, C_d, var)
    gamma_d = gamma_d0 @ sqrtm(sigmav)

    # Square root of steady-state variance matrix (for initialization)
    sigma_x = solve_discrete_lyapunov(A_d, gamma_d @ gamma_d.T)
    rac_sigma_x = sqrtm(sigma_x)
    rac_sigma_x = np.real(rac_sigma_x)

    return A_d, gamma_d, C_d, rac_sigma_x


class VibrationGenerator:

    def __init__(self, model_c, fs, seed = 0):

        self.fs = fs
        self.seed = seed
        self.A_c = model_c[1].data
        self.C_c = model_c[2].data
        self.gamma_c = model_c[3].data
        self.var = model_c[4].data
        self.A_d, self.gamma_d, self.C_d, self.rac_sigma_x = rediscret(self.A_c, self.gamma_c, self.C_c, self.var, self.fs)
        self.nxvib = self.A_d.shape[0]

        self.wn_gen = np.random.default_rng(seed=seed) # white noise generator
        self.xvib = self.rac_sigma_x @ self.wn_gen.standard_normal((self.nxvib, 1))

    def step(self):
        vk = self.wn_gen.standard_normal((self.C_d.shape[0], 1))  # White noise input
        self.xvib = self.A_d @ self.xvib + self.gamma_d @ vk
        out = (self.C_d @ self.xvib).squeeze()
        return out


    def set_seed(self, seed):
        self.wn_gen = np.random.default_rng(seed=seed) # white noise generator
        self.xvib = self.rac_sigma_x @ self.wn_gen.standard_normal((self.nxvib, 1))

if __name__ == "__main__":

    model_path= '../studies/vibration/models/modelevib0606tip.fits'
    model_c =  pfits.open(model_path)
    fs = 3000
    seed = 0

    gen = VibrationGenerator(model_c, fs, seed)

    print(gen.step())