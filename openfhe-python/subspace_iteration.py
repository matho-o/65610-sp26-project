import numpy as np
from openfhe import *
import matrix

# ── plaintext helpers ──────────────────────────────────────────────────────────

def subspace_iteration_plaintext(A, k, num_iterations=2, seed=42):
    """
    Plaintext reference implementation of subspace iteration.
    A: 2D numpy array (m x n)
    k: number of singular vectors to find
    num_iterations: more = more accurate, use this to study the tradeoff
    """
    
    np.random.seed(seed)  # add this
    m, n = A.shape

    # initialise with random n x k matrix, orthonormalized
    Q, _ = np.linalg.qr(np.random.randn(n, k))

    for _ in range(num_iterations):
        Z = A @ Q          # m x k
        Z, _ = np.linalg.qr(Z)
        Q = A.T @ Z        # n x k
        Q, _ = np.linalg.qr(Q)

    # recover singular values and vectors from final subspace
    Z = A @ Q              # m x k
    U_small, s, Vt_small = np.linalg.svd(Z, full_matrices=False)
    U = U_small[:, :k]
    V = Q @ Vt_small.T     # n x k
    return U, s[:k], V[:, :k].T


# ── polynomial approximation of 1/sqrt(x) ─────────────────────────────────────

def poly_invsqrt_coeffs(degree, x_min, x_max):
    """
    Compute Chebyshev approximation coefficients for 1/sqrt(x)
    on the interval [x_min, x_max].
    degree: higher = more accurate, more depth. start with 3 or 5.
    """
    # sample points in Chebyshev nodes on [x_min, x_max]
    nodes = np.cos(np.pi * (2 * np.arange(1, degree + 2) - 1) / (2 * (degree + 1)))
    x_nodes = 0.5 * (x_max - x_min) * nodes + 0.5 * (x_max + x_min)
    y_nodes = 1.0 / np.sqrt(x_nodes)
    # fit a polynomial through those points
    coeffs = np.polyfit(x_nodes, y_nodes, degree)
    return coeffs  # numpy poly1d convention, highest degree first


def eval_poly(coeffs, x):
    """Evaluate polynomial with given coeffs at x (scalar or array)."""
    return np.polyval(coeffs, x)


# ── FHE normalization ──────────────────────────────────────────────────────────

def fhe_inner_product_to_scalar(cc, ct, n):
    """
    Rotate-and-sum to accumulate the sum of first n slots into every slot.
    This gives <v,v> (the squared norm) broadcast across all slots.
    Costs log2(n) rotations and additions — no multiplications, so depth-free.
    """
    result = ct
    step = 1
    while step < n:
        rotated = cc.EvalRotate(result, step)
        result = cc.EvalAdd(result, rotated)
        step *= 2
    return result  # every slot now holds sum of all slots


def fhe_normalize(cc, pub_key, ct_v, n, poly_degree=3, x_min=0.01, x_max=10.0):
    """
    Normalize an encrypted vector v of length n using a polynomial
    approximation of 1/sqrt(||v||^2).

    poly_degree: degree of the approximating polynomial (higher = more depth)
    x_min, x_max: expected range of ||v||^2 — tune these to your data
                  for better approximation accuracy

    Depth cost: 1 (for squaring via EvalMult to get v*v)
              + ceil(log2(poly_degree+1)) for polynomial evaluation
              + 1 for the final scaling multiply
              ≈ poly_degree 3 costs roughly 3-4 levels total
    """

    slot_size = cc.GetRingDimension() // 2

    # step 1: compute v * v elementwise, then sum → ||v||^2 in every slot
    ct_sq = cc.EvalMult(ct_v, ct_v)                    # depth +1
    ct_norm_sq = fhe_inner_product_to_scalar(cc, ct_sq, n)  # depth +0

    # step 2: evaluate polynomial approximation of 1/sqrt(x) at ct_norm_sq
    coeffs = poly_invsqrt_coeffs(poly_degree, x_min, x_max)

    # Horner evaluation: coeffs are [c0, c1, ..., cd] highest degree first
    pt = lambda c: cc.MakeCKKSPackedPlaintext([float(c)] * slot_size)
    
    ct_poly = pt(coeffs[0])  # start with highest degree coefficient as plaintext
    for c in coeffs[1:]:
        ct_poly = cc.EvalMult(ct_poly, ct_norm_sq)  # multiply first
        ct_poly = cc.EvalAdd(ct_poly, pt(c))         # then add next coefficient

    # step 3: scale v by the approximated 1/||v||
    ct_normalized = cc.EvalMult(ct_v, ct_poly)         # depth +1
    return ct_normalized


# ── FHE gram-schmidt for k columns ────────────────────────────────────────────

def fhe_orthonormalize(cc, pub_key, ct_columns, n, poly_degree=3,
                        x_min=1.0, x_max=5.0, use_gram_schmidt=False):
    slot_size = cc.GetRingDimension() // 2
    pt_one = cc.MakeCKKSPackedPlaintext([1.0] * slot_size)
    
    result = []
    for i, ct_qi in enumerate(ct_columns):
        if use_gram_schmidt:
            for ct_qj in result:
                while ct_qj.GetLevel() > ct_qi.GetLevel():
                    ct_qi = cc.EvalMult(ct_qi, pt_one)
                while ct_qi.GetLevel() > ct_qj.GetLevel():
                    ct_qj = cc.EvalMult(ct_qj, pt_one)
                ct_dot_elems = cc.EvalMult(ct_qi, ct_qj)
                ct_dot = fhe_inner_product_to_scalar(cc, ct_dot_elems, n)
                ct_proj = cc.EvalMult(ct_qj, ct_dot)
                ct_qi = cc.EvalSub(ct_qi, ct_proj)
        ct_qi = fhe_normalize(cc, pub_key, ct_qi, n, poly_degree=poly_degree,
                               x_min=x_min, x_max=x_max)
        result.append(ct_qi)
    return result

# ── main FHE subspace iteration ────────────────────────────────────────────────

def subspace_iteration_fhe(cc, pub_key, ct_A, ct_At, n, k,
                            num_iterations=1, poly_degree=3,
                            x_min_z=0.1, x_max_z=2.0,
                            x_min_q=1.0, x_max_q=5.0,
                            seed=42, use_gram_schmidt=False):
    """
    FHE subspace iteration for truncated SVD.

    ct_A  : encrypted n x n matrix A (row-major, periodically repeated)
    ct_At : encrypted transpose of A
    n     : matrix dimension (square for now)
    k     : number of singular vectors
    num_iterations: 1 is cheapest; each iteration costs roughly:
                    2 matrix multiplies + k * fhe_normalize depth
    poly_degree: degree of 1/sqrt approximation (3 is cheapest reasonable choice)

    Returns: (ct_U_cols, ct_s, ct_Vt_rows) — lists of k ciphertexts
             each encoding one singular vector
    """
    
    np.random.seed(seed)  # add this

    slot_size = cc.GetRingDimension() // 2

    Q_plain = np.random.randn(n, k)
    Q_plain, _ = np.linalg.qr(Q_plain)

    ct_Q_cols = []
    for j in range(k):
        col_data = Q_plain[:, j].tolist() * (slot_size // n)
        ct_Q_cols.append(cc.Encrypt(pub_key, cc.MakeCKKSPackedPlaintext(col_data)))

    print(f"Q initial level: {ct_Q_cols[0].GetLevel()}")

    for iteration in range(num_iterations):
        print(f"\n  subspace iteration {iteration + 1}/{num_iterations}")

        ct_Z_cols = []
        for idx, ct_q in enumerate(ct_Q_cols):
            ct_z = matrix.matrix_multiply(ct_A, ct_q, n, cc, pub_key)
            print(f"    after A@q[{idx}], level: {ct_z.GetLevel()}")
            ct_Z_cols.append(ct_z)

        # Z vectors have smaller norms — use z range
        ct_Z_cols = fhe_orthonormalize(cc, pub_key, ct_Z_cols, n,
                                        poly_degree=poly_degree,
                                        x_min=x_min_z, x_max=x_max_z, use_gram_schmidt=use_gram_schmidt)
        print(f"    after ortho Z, level: {ct_Z_cols[0].GetLevel()}")

        ct_Q_cols = []
        for idx, ct_z in enumerate(ct_Z_cols):
            ct_q = matrix.matrix_multiply(ct_At, ct_z, n, cc, pub_key)
            print(f"    after At@z[{idx}], level: {ct_q.GetLevel()}")
            ct_Q_cols.append(ct_q)

        # Q vectors have larger norms — use q range
        ct_Q_cols = fhe_orthonormalize(cc, pub_key, ct_Q_cols, n,
                                        poly_degree=poly_degree,
                                        x_min=x_min_q, x_max=x_max_q)
        print(f"    after ortho Q, level: {ct_Q_cols[0].GetLevel()}")

    # final Z = A @ Q
    ct_Z_cols = []
    for ct_q in ct_Q_cols:
        ct_z = matrix.matrix_multiply(ct_A, ct_q, n, cc, pub_key)
        ct_Z_cols.append(ct_z)
    print(f"\n  final Z level: {ct_Z_cols[0].GetLevel()}")

    ct_U_cols = fhe_orthonormalize(cc, pub_key, ct_Z_cols, n,
                                    poly_degree=poly_degree,
                                    x_min=x_min_z, x_max=x_max_z)
    ct_Vt_rows = ct_Q_cols
    ct_s = None

    return ct_U_cols, ct_s, ct_Vt_rows