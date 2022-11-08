import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
def plot_DM_shape(n_act,n_act_eff):
    x_pos,y_pos = compute_xy_pos()
    DM_shape = np.zeros((n_act, n_act))
    for i in range(n_act_eff):
        DM_shape[y_pos[i],x_pos[i]] = 1
    plt.figure(figsize=(7, 7))
    plt.imshow(DM_shape)
    plt.show()
def plot_modes(M,n,n_act,n_act_eff):
    DM_shape = np.zeros((n_act,n_act,n**2))
    x_pos, y_pos = compute_xy_pos()
    plt.figure(figsize=(7, 7))
    for i in range(n):
        for j in range(n):
            for k in range(n_act_eff):
                DM_shape[y_pos[k], x_pos[k],i*n+j] = M[k,i*n+j]
            plt.subplot(n,n,i*n+j+1)
            plt.imshow(DM_shape[:,:,i*n+j])
            # plt.magma()

def plot_PCA(PCA, n, pupil_diam):
    plt.figure(figsize=(7, 7))
    for i in range(n):
        for j in range(n):
            phase_shape = PCA[:,i*n+j].reshape((pupil_diam,pupil_diam))
            plt.subplot(n, n, i * n + j + 1)
            plt.imshow(phase_shape)
            plt.magma()
def compute_xy_pos():
    x_pos = np.load('x_pos.npy')
    y_pos = np.load('y_pos.npy')
    x_pos -= np.min(x_pos)
    y_pos -= np.min(y_pos)
    x_pos /= 4
    y_pos /= 4
    return x_pos.astype(int),y_pos.astype(int)

def compute_modal_basis(inf_mat,n_lambda):
    Gamma = inf_mat.T @ inf_mat
    D, M, dumm = np.linalg.svd(Gamma)

    sort_perm = np.argsort(-D)
    D = D[sort_perm]
    M = M[:, sort_perm]

    D = D[:n_lambda]
    M = M[:,:n_lambda]

    M /= np.sqrt(D)

    cond = D[0]/D[-1]
    noise_prop = np.sum(1/D)

    print('Cond = {0:.2f}, Noise prop = {1:.2f}'.format(cond,noise_prop))
    return M

def phase_PCA(dist_mat):
    Gamma = dist_mat @ dist_mat.T
    D, M = np.linalg.eigh(Gamma)
    sort_perm = np.argsort(-D)
    D = D[sort_perm]
    M = M[:, sort_perm]

    return M

def statistical_diag(M,dist_mat):

    dist_mean = np.mean(dist_mat,axis=1)
    dist_mat -= dist_mean[:,None]
    dist_std = np.std(dist_mat,axis=1)
    dist_mat /= dist_std[:,None]

    # dist_mat = gaussian_filter(dist_mat, sigma=30)
    dist_cov = np.cov(dist_mat)

    # dist_cov = gaussian_filter(dist_cov, sigma=7)

    C = np.linalg.inv(M) @ dist_cov @ np.linalg.inv(M.T)
    Sigma, A = np.linalg.eigh(C)

    sort_perm = np.argsort(-Sigma)
    Sigma = Sigma[sort_perm]
    A = A[:, sort_perm]

    B = M @ A
    return B

if __name__ == '__main__':
    n_act = 41
    n_act_eff = 1276
    pupil_diam = 41

    inf_mat = np.load('inf_mat.npy')
    dist_mat = np.load('dist_matrix_act.npy')
    # PCA = phase_PCA(dist_mat)
    # M = compute_modal_basis(dist_mat.T,200)
    # np.save('M.npy', M)
    M = compute_modal_basis(inf_mat, n_act_eff)
    np.save('M.npy', M[:,:])
    # B = statistical_diag(M, dist_mat)
    # np.save('B.npy', B)
    # plot_PCA(PCA, 4, pupil_diam)
    plot_modes(M, 3, n_act, n_act_eff)
    # plot_modes(B, 3, n_act, n_act_eff)


    basis = np.load('inf_mat_KL2V.npy')
    plot_modes(basis, 3, n_act, n_act_eff)
    plt.show()