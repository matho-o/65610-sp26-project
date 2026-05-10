import numpy as np
from openfhe import *
import matrix


def power_iteration_plaintext(A, num_iterations=1, seed=42):
    np.random.seed(seed)

    m, n = A.shape
    v = np.random.randn(n)
    v = v / np.linalg.norm(v)

    for _ in range(num_iterations):
        z = A @ v
        z = z / np.linalg.norm(z)

        v = A.T @ z
        v = v / np.linalg.norm(v)

    Av = A @ v
    sigma = np.linalg.norm(Av)
    u = Av / sigma

    return u, sigma, v

# delete 3 funcs later after merging delete and import from subspace_iterations
def poly_invsqrt_coeffs(degree, x_min, x_max):
    nodes = np.cos(
        np.pi * (2 * np.arange(1, degree + 2) - 1) / (2 * (degree + 1))
    )
    x_nodes = 0.5 * (x_max - x_min) * nodes + 0.5 * (x_max + x_min)
    y_nodes = 1.0 / np.sqrt(x_nodes)
    return np.polyfit(x_nodes, y_nodes, degree)


def pack_vector_as_repeated_columns(v):
    n = len(v)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = v[i]
    return M.flatten().tolist()


def extract_vector_from_repeated_columns(values, n):
    M = np.array(values[: n * n]).reshape(n, n)
    return M[:, 0]


def fhe_sum_first_slots(cc, ct, total_slots_to_sum):
    result = ct
    step = 1

    while step < total_slots_to_sum:
        rotated = cc.EvalRotate(result, step)
        result = cc.EvalAdd(result, rotated)
        step *= 2

    return result


def fhe_inv_sqrt_from_norm_sq(
    cc,
    ct_norm_sq,
    poly_degree=3,
    x_min=0.01,
    x_max=5.0,
):
    slot_size = cc.GetRingDimension() // 2
    coeffs = poly_invsqrt_coeffs(poly_degree, x_min, x_max)

    def pt(c):
        return cc.MakeCKKSPackedPlaintext([float(c)] * slot_size)

    ct_poly = pt(coeffs[0])
    for c in coeffs[1:]:
        ct_poly = cc.EvalMult(ct_poly, ct_norm_sq)
        ct_poly = cc.EvalAdd(ct_poly, pt(c))

    return ct_poly


def fhe_norm_sq_repeated_columns(cc, ct_v_mat, n):

    ct_sq = cc.EvalMult(ct_v_mat, ct_v_mat)
    ct_sum = fhe_sum_first_slots(cc, ct_sq, n * n)
    return cc.EvalMult(ct_sum, 1.0 / n)


def fhe_normalize_repeated_columns(
    cc,
    pub_key,
    ct_v_mat,
    n,
    poly_degree=3,
    x_min=0.01,
    x_max=5.0,
):
    ct_norm_sq = fhe_norm_sq_repeated_columns(cc, ct_v_mat, n)

    ct_inv_norm = fhe_inv_sqrt_from_norm_sq(
        cc,
        ct_norm_sq,
        poly_degree=poly_degree,
        x_min=x_min,
        x_max=x_max,
    )

    return cc.EvalMult(ct_v_mat, ct_inv_norm)


def fhe_inner_product_repeated_columns(cc, ct_x_mat, ct_y_mat, n):

    ct_prod = cc.EvalMult(ct_x_mat, ct_y_mat)
    ct_sum = fhe_sum_first_slots(cc, ct_prod, n * n)
    return cc.EvalMult(ct_sum, 1.0 / n)


def encrypted_matvec_repeated_columns(ct_A, ct_v_mat, n, cc, pub_key):
    return matrix.matrix_multiply(ct_A, ct_v_mat, n, cc, pub_key)


def power_iteration_fhe(
    cc,
    pub_key,
    ct_A,
    ct_At,
    n,
    num_iterations=1,
    poly_degree=3,
    x_min_z=0.01,
    x_max_z=5.0,
    x_min_v=0.01,
    x_max_v=5.0,
    seed=42,
):
    """
    FHE rank 1 SVD using encrypted A.

    Returns:
        ct_u, ct_sigma, ct_v
    """
    np.random.seed(seed)
    slot_size = cc.GetRingDimension() // 2

    v0 = np.random.randn(n)
    v0 = v0 / np.linalg.norm(v0)

    packed_v0_mat = pack_vector_as_repeated_columns(v0)

    ct_v = cc.Encrypt(
        pub_key,
        cc.MakeCKKSPackedPlaintext(
            packed_v0_mat * (slot_size // (n * n))
        ),
    )

    print(f"initial v-matrix level: {ct_v.GetLevel()}")

    for it in range(num_iterations):
        print(f"\n  power iteration {it + 1}/{num_iterations}")

        ct_z = encrypted_matvec_repeated_columns(
            ct_A, ct_v, n, cc, pub_key
        )
        print(f"    after A@v, level: {ct_z.GetLevel()}")

        ct_z = fhe_normalize_repeated_columns(
            cc,
            pub_key,
            ct_z,
            n,
            poly_degree=poly_degree,
            x_min=x_min_z,
            x_max=x_max_z,
        )
        print(f"    after normalize z, level: {ct_z.GetLevel()}")

        ct_v = encrypted_matvec_repeated_columns(
            ct_At, ct_z, n, cc, pub_key
        )
        print(f"    after At@z, level: {ct_v.GetLevel()}")

        ct_v = fhe_normalize_repeated_columns(
            cc,
            pub_key,
            ct_v,
            n,
            poly_degree=poly_degree,
            x_min=x_min_v,
            x_max=x_max_v,
        )
        print(f"    after normalize v, level: {ct_v.GetLevel()}")

    # Av
    ct_Av = encrypted_matvec_repeated_columns(
        ct_A, ct_v, n, cc, pub_key
    )
    print(f"    final after A@v, level: {ct_Av.GetLevel()}")

    # u = Av / ||Av||
    ct_u = fhe_normalize_repeated_columns(
        cc,
        pub_key,
        ct_Av,
        n,
        poly_degree=poly_degree,
        x_min=x_min_z,
        x_max=x_max_z,
    )
    print(f"    final u level: {ct_u.GetLevel()}")

    # sigma = u^T Av
    # This avoids computing sqrt(||Av||^2), which was inaccurate.
    ct_sigma = fhe_inner_product_repeated_columns(
        cc,
        ct_u,
        ct_Av,
        n,
    )
    print(f"    final sigma level: {ct_sigma.GetLevel()}")

    return ct_u, ct_sigma, ct_v


def estimate_power_iteration_depth(num_iterations=1, poly_degree=3):
    matmul_depth = 2
    norm_depth = poly_degree + 2
    dot_depth = 2  # multiply + scalar divide/mask-ish

    per_iter = 2 * matmul_depth + 2 * norm_depth
    final_step = matmul_depth + norm_depth + dot_depth

    return num_iterations * per_iter + final_step
