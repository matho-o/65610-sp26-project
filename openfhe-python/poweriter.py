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
    x_max=9.0,
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
    x_max=9.0,
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
    x_max_z=9.0,
    x_min_v=0.01,
    x_max_v=9.0,
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


    ct_Av = encrypted_matvec_repeated_columns(
        ct_A, ct_v, n, cc, pub_key
    )
    print(f"    final after A@v, level: {ct_Av.GetLevel()}")

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
    dot_depth = 2
    per_iter = 2 * matmul_depth + 2 * norm_depth
    final_step = matmul_depth + norm_depth + dot_depth

    return num_iterations * per_iter + final_step


def make_plain_one(cc):
    slot_size = cc.GetRingDimension() // 2
    return cc.MakeCKKSPackedPlaintext([1.0] * slot_size)


def raise_to_level(cc, ct, target_level):
    """
    Increase ct.GetLevel() until it matches target_level by multiplying by 1.
    """
    pt_one = make_plain_one(cc)

    while ct.GetLevel() < target_level:
        ct = cc.EvalMult(ct, pt_one)

    return ct


def match_levels(cc, ct_a, ct_b):
    """
    Make two ciphertexts have the same level.
    """
    la = ct_a.GetLevel()
    lb = ct_b.GetLevel()

    if la < lb:
        ct_a = raise_to_level(cc, ct_a, lb)
    elif lb < la:
        ct_b = raise_to_level(cc, ct_b, la)

    return ct_a, ct_b


def encrypted_outer_product_from_repeated_columns(cc, pub_key, ct_u, ct_v, n):
    """
    ct_u stores u as repeated columns:
        slot(i,j) = u_i

    ct_v stores v as repeated columns:
        slot(i,j) = v_i

    To form u v^T, we need:
        slot(i,j) = u_i v_j

    So transpose ct_v to get repeated rows:
        slot(i,j) = v_j

    Then multiply slotwise.
    """
    ct_v_rows = matrix.transpose(ct_v, n, cc, pub_key)

    ct_u, ct_v_rows = match_levels(cc, ct_u, ct_v_rows)

    return cc.EvalMult(ct_u, ct_v_rows)


def deflate_encrypted_matrix(cc, pub_key, ct_A, ct_u, ct_sigma, ct_v, n):
    """
    Deflate:
        A <- A - sigma * u v^T

    All encrypted.
    """
    ct_outer = encrypted_outer_product_from_repeated_columns(
        cc,
        pub_key,
        ct_u,
        ct_v,
        n,
    )

    ct_outer, ct_sigma = match_levels(cc, ct_outer, ct_sigma)

    ct_rank1 = cc.EvalMult(ct_outer, ct_sigma)

    ct_A, ct_rank1 = match_levels(cc, ct_A, ct_rank1)

    ct_A_new = cc.EvalSub(ct_A, ct_rank1)

    return ct_A_new


def rank_k_svd_fhe(
    cc,
    pub_key,
    ct_A,
    n,
    k=2,
    num_iterations=1,
    poly_degree=3,
    x_min_z=0.01,
    x_max_z=5.0,
    x_min_v=0.01,
    x_max_v=5.0,
    seed=42,
):
    """
    Encrypted rank-k SVD using repeated rank-1 power iteration + deflation.

    For i = 1..k:
        compute sigma_i, u_i, v_i
        A <- A - sigma_i u_i v_i^T

    Returns:
        ct_U_list: list of encrypted u_i
        ct_S_list: list of encrypted sigma_i
        ct_V_list: list of encrypted v_i
        ct_A_work: final deflated encrypted matrix
    """
    ct_A_work = ct_A

    ct_U_list = []
    ct_S_list = []
    ct_V_list = []

    for i in range(k):
        print(f"\n========== rank component {i + 1}/{k} ==========")

        ct_At_work = matrix.transpose(ct_A_work, n, cc, pub_key)

        ct_u, ct_sigma, ct_v = power_iteration_fhe(
            cc,
            pub_key,
            ct_A_work,
            ct_At_work,
            n=n,
            num_iterations=num_iterations,
            poly_degree=poly_degree,
            x_min_z=x_min_z,
            x_max_z=x_max_z,
            x_min_v=x_min_v,
            x_max_v=x_max_v,
            seed=seed + i,
        )

        ct_U_list.append(ct_u)
        ct_S_list.append(ct_sigma)
        ct_V_list.append(ct_v)

        print(f"  deflating A by component {i + 1}")
        ct_A_work = deflate_encrypted_matrix(
            cc,
            pub_key,
            ct_A_work,
            ct_u,
            ct_sigma,
            ct_v,
            n,
        )
        print(f"  after deflation, A level: {ct_A_work.GetLevel()}")

    return ct_U_list, ct_S_list, ct_V_list, ct_A_work
