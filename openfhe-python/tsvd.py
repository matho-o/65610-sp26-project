from openfhe import *
from matrix import matrix_multiply, transpose
import random

def truncated_svd(A:Ciphertext, n:int, cc:CryptoContext, pk:PublicKey, alpha:float=0.01, iters:int=4):
    """
    Halko randomized algorithm for TSVD
    A is an n by n matrix to find the TSVD of
    """

    # Generate random Gaussian matrix and encrypt
    rand_ptx = [random.gauss(0, 1) for i in range(n ** 2)]
    rand_ctx = cc.Encrypt(pk, cc.MakeCKKSPackedPlaintext(rand_ptx))

    # Compute random sample matrix from A and its transpose
    Y = matrix_multiply(A, rand_ctx, n, cc, pk)
    Y_T = transpose(Y, n, cc, pk)

    # Compute dense covariance matrix
    M = matrix_multiply(Y_T, Y, n, cc, pk)

    print("Cov M computed")

    # Inverse square root by Newton-Schulz iteration
    # Initialize guess to scaled identity matrix
    X0_ptx = [(alpha if (i%n) == (i//n) else 0) for i in range(n ** 2)]
    X = cc.Encrypt(pk, cc.MakeCKKSPackedPlaintext(X0_ptx))

    arr_3I = [(3 if (i%n) == (i//n) else 0) for i in range(n ** 2)]
    ptx_3I = cc.MakeCKKSPackedPlaintext(arr_3I)

    ctx_3I = cc.Encrypt(pk, ptx_3I)

    for i in range(iters):

        # X_sq = matrix_multiply(X, X, n, cc, pk)
        # M_X_sq = matrix_multiply(M, X_sq, n, cc, pk)
        # diff = cc.EvalSub(ptx_3I, M_X_sq)
        # X = matrix_multiply(X, diff, n, cc, pk)
        
        # optimizing for mult depth
        X_sq = matrix_multiply(X, X, n, cc, pk)
        X_M = matrix_multiply(X, M, n, cc, pk)

        t1 = matrix_multiply(X, ctx_3I, n, cc, pk, c=0.5)
        t2 = matrix_multiply(X_M, X_sq, n, cc, pk, c=0.5)

        X = cc.EvalSub(t1, t2)

        # X = cc.EvalMult(X, 0.5)

        print("Newton-Schulz iter", i)

    # Compute orthonormal basis of random sample matrix
    Q = matrix_multiply(Y, X, n, cc, pk)
    Q_T = transpose(Q, n, cc, pk)

    # project into subspace and compute standard SVD of projection B
    B = matrix_multiply(Q_T, A, n, cc, pk)

    # TODO update
    return B
