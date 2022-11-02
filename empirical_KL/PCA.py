import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors



if __name__ == '__main__':

    # S2A = np.load('S2A_I.npy')
    # S2A/= np.max(S2A)
    # S2A = S2A[:1284,:]
    # U,S,V = np.linalg.svd(S2A)
    # np.save('U.npy', U)

    # A2S = np.load('inf_mat.npy')
    # # A2S = A2S[:,:1284]
    # # # A2S/= np.max(A2S)
    # # U2,S2,V2 = np.linalg.svd(A2S)
    # #
    # # np.save('U.npy',U2)
    #
    # Sigma = A2S.T @ A2S
    # # Sigma -= np.mean(Sigma)
    # # Sigma[np.abs(Sigma)<0.1] = 0
    # # Sigma /= np.mean(Sigma)
    #
    # D,M = np.linalg.eigh(Sigma)
    # #
    # sort_perm = np.argsort(-D)
    # D = D[sort_perm]
    # M = M[:, sort_perm]
    # #
    # # # D = np.sqrt(D)
    # np.save('M.npy', M)
    M = np.load('M.npy')
    # plt.imshow(Sigma)
    plt.plot(M[:,1])
    plt.show()