from openfhe import *
import matrix
import os
import numpy as np
import random
import poweriter
from time import perf_counter
# --- YOUR ORIGINAL PARAMETER SETUP ---
parameters = CCParamsCKKSRNS()

secret_key_dist = SecretKeyDist.UNIFORM_TERNARY
parameters.SetSecretKeyDist(secret_key_dist)

parameters.SetSecurityLevel(SecurityLevel.HEStd_NotSet)
parameters.SetRingDim(1 << 16)

rescale_tech = ScalingTechnique.FLEXIBLEAUTO
dcrt_bits = 59
first_mod = 60

parameters.SetScalingModSize(dcrt_bits)
parameters.SetScalingTechnique(rescale_tech)
parameters.SetFirstModSize(first_mod)

level_budget = [4, 4]
#levels_available_after_bootstrap = 3
#depth = levels_available_after_bootstrap + FHECKKSRNS.GetBootstrapDepth(level_budget, secret_key_dist)

parameters.SetMultiplicativeDepth(130)

cryptocontext = GenCryptoContext(parameters)
cryptocontext.Enable(PKESchemeFeature.PKE)
cryptocontext.Enable(PKESchemeFeature.KEYSWITCH)
cryptocontext.Enable(PKESchemeFeature.LEVELEDSHE)
cryptocontext.Enable(PKESchemeFeature.ADVANCEDSHE)
cryptocontext.Enable(PKESchemeFeature.FHE)

ring_dim = cryptocontext.GetRingDimension()
num_slots = int(ring_dim / 2)
print(f"CKKS is using ring dimension {ring_dim} with {num_slots} slots.")

# --- KEY GENERATION ---
cryptocontext.EvalBootstrapSetup(level_budget)
key_pair = cryptocontext.KeyGen()
cryptocontext.EvalMultKeyGen(key_pair.secretKey)
# cryptocontext.EvalBootstrapKeyGen(key_pair.secretKey, num_slots)

# --- ROTATION KEY GENERATION (Using your requested range) ---
rot_indices = [-128, -64, -32, -16, -8, -4, -2, -1, 1, 2, 4, 8, 16, 32, 64, 128]
cryptocontext.EvalAtIndexKeyGen(key_pair.secretKey, rot_indices)
print(f"Rotation keys generated for: {rot_indices}")

# --- FULL TEST SUITE ---

# 1. Prepare Test Data
size = 4
matrix_a_data = [random.random() for i in range(size*size)]
matrix_b_data = [random.random() for i in range(size*size)]

np_a = np.array(matrix_a_data).reshape(size, size)
np_b = np.array(matrix_b_data).reshape(size, size)

# 2. Encrypt
cipher_a = cryptocontext.Encrypt(key_pair.publicKey, cryptocontext.MakeCKKSPackedPlaintext(matrix_a_data*(num_slots//(size*size))))
cipher_b = cryptocontext.Encrypt(key_pair.publicKey, cryptocontext.MakeCKKSPackedPlaintext(matrix_b_data*(num_slots//(size*size))))

# 4. Matrix Multiplication Test
print("\n--- Step 2: Testing matrix.matrix_multiply ---")
try:
    # Format: (ct, ct, cc, key_pair, size)
    res_mult = matrix.matrix_multiply(cipher_a, cipher_b, size, cryptocontext, key_pair.publicKey)

    dec_mult = cryptocontext.Decrypt(res_mult, key_pair.secretKey)
    dec_mult.SetLength(size * size)
    actual_mult = [x.real for x in dec_mult.GetCKKSPackedValue()]
    # Calculate Expected
    np_a = np.array(matrix_a_data).reshape(size, size)
    np_b = np.array(matrix_b_data).reshape(size, size)
    expected_mult = np.matmul(np_a, np_b).flatten().tolist()

    print(f"Expected: {expected_mult}")
    print(f"Actual:   {actual_mult}")

    if np.allclose(actual_mult, expected_mult, atol=0.001):
        print(" Multiplication Success!")
except Exception as e:
    print(f"Multiplication test failed: {e}")

# 5. Transpose Test
print("\n--- Step 3: Testing matrix.transpose ---")
try:
    # Format: (ct, cc, key_pair, size)
    res_trans = matrix.transpose(cipher_a, size, cryptocontext, key_pair.publicKey)

    dec_trans = cryptocontext.Decrypt(res_trans, key_pair.secretKey)
    dec_trans.SetLength(size * size)
    actual_trans = [x.real for x in dec_trans.GetCKKSPackedValue()]

    expected_trans = np_a.T.flatten().tolist()

    print(f"Expected: {expected_trans}")
    print(f"Actual:   {actual_trans}")

    if np.allclose(actual_trans, expected_trans, atol=0.001):
        print(" Transpose Success!")
except Exception as e:
    print(f"Transpose test failed: {e}")




# 6. power iteration test
print("\n--- Step 4: Testing SVD by power iteration ---")
try:
    # Normalize A for stability of polynomial inverse-sqrt
    frob_norm = np.linalg.norm(np_a, "fro")
    A_normalized = np_a / frob_norm

    cipher_a_svd = cryptocontext.Encrypt(
        key_pair.publicKey,
        cryptocontext.MakeCKKSPackedPlaintext(
            A_normalized.flatten().tolist() * (num_slots // (size * size))
        ),
    )

    ct_At = matrix.transpose(
        cipher_a_svd,
        size,
        cryptocontext,
        key_pair.publicKey,
    )

    fhe_iters = 3
    poly_deg = 3



    # Plaintext ref
    u_pt, sigma_pt, v_pt = poweriter.power_iteration_plaintext(
        A_normalized,
        num_iterations=fhe_iters,
        seed=42,
    )

    U_ref, s_ref, Vt_ref = np.linalg.svd(A_normalized, full_matrices=False)

    u_ref = U_ref[:, 0]
    sigma_ref = s_ref[0]
    v_ref = Vt_ref[0, :]

    print("\nPlaintext power iteration:")
    print("sigma PI:", sigma_pt)
    print("sigma SVD:", sigma_ref)
    print("|<u_pt,u_ref>|:", abs(np.dot(u_pt, u_ref)))
    print("|<v_pt,v_ref>|:", abs(np.dot(v_pt, v_ref)))

    X_MIN_Z, X_MAX_Z = 0.01, 2.0
    X_MIN_V, X_MAX_V = 0.01, 7.0

    start = perf_counter()

    rank_k = 4

    ct_U_list, ct_S_list, ct_V_list, ct_A_deflated = poweriter.rank_k_svd_fhe(
        cryptocontext,
        key_pair.publicKey,
        cipher_a_svd,
        n=size,
        k=rank_k,
        num_iterations=fhe_iters,
        poly_degree=poly_deg,
        x_min_z=X_MIN_Z,
        x_max_z=X_MAX_Z,
        x_min_v=X_MIN_V,
        x_max_v=X_MAX_V,
        seed=42,
    )

    elapsed = perf_counter() - start

    # Decrypt u: repeated-column matrix, extract first column
    dec_u = cryptocontext.Decrypt(ct_u, key_pair.secretKey)
    dec_u.SetLength(size * size)

    u_fhe = poweriter.extract_vector_from_repeated_columns(
        dec_u.GetRealPackedValue(),
        size,
    )
    u_fhe = u_fhe / np.linalg.norm(u_fhe)

    # Decrypt v: repeated-column matrix, extract first column
    dec_v = cryptocontext.Decrypt(ct_v, key_pair.secretKey)
    dec_v.SetLength(size * size)

    v_fhe = poweriter.extract_vector_from_repeated_columns(
        dec_v.GetRealPackedValue(),
        size,
    )
    v_fhe = v_fhe / np.linalg.norm(v_fhe)

    # Decrypt sigma: scalar broadcast
    dec_s = cryptocontext.Decrypt(ct_s, key_pair.secretKey)
    dec_s.SetLength(1)
    sigma_fhe = dec_s.GetRealPackedValue()[0]

    print("\nFHE rank-1 SVD:")
    print("sigma_fhe:", sigma_fhe)
    print("sigma_ref:", sigma_ref)
    print("sigma_pt:", sigma_pt)
    print("|<u_fhe,u_ref>|:", abs(np.dot(u_fhe, u_ref)))
    print("|<v_fhe,v_ref>|:", abs(np.dot(v_fhe, v_ref)))
    print("|<u_fhe,u_pt>|:", abs(np.dot(u_fhe, u_pt)))
    print("|<v_fhe,v_pt>|:", abs(np.dot(v_fhe, v_pt)))
    print("runtime seconds:", elapsed)

except Exception as e:
    print(f"Power iteration test failed: {e}")


# 7. power iteration test: rank-k SVD
print("\n--- Step 4: Testing rank-k SVD by deflation ---")
try:
    rank_k = 4
    fhe_iters = 3
    poly_deg = 3

    # Normalize A for  stability
    A_np = np.array(matrix_a_data).reshape(size, size)
    A_normalized = A_np / np.linalg.norm(A_np, "fro")

    cipher_a_svd = cryptocontext.Encrypt(
        key_pair.publicKey,
        cryptocontext.MakeCKKSPackedPlaintext(
            A_normalized.flatten().tolist() * (num_slots // (size * size))
        ),
    )

    U_ref, s_ref, Vt_ref = np.linalg.svd(A_normalized, full_matrices=False)

    # Plaintext  ref
    A_work_plain = A_normalized.copy()
    U_pt_list, S_pt_list, V_pt_list = [], [], []

    for i in range(rank_k):
        u_pt, sigma_pt, v_pt = poweriter.power_iteration_plaintext(
            A_work_plain,
            num_iterations=fhe_iters,
            seed=42 + i,
        )

        U_pt_list.append(u_pt)
        S_pt_list.append(sigma_pt)
        V_pt_list.append(v_pt)

        A_work_plain = A_work_plain - sigma_pt * np.outer(u_pt, v_pt)

    print("\nPlaintext deflation reference:")
    for i in range(rank_k):
        print(f"component {i}")
        print(f"  sigma_pt : {S_pt_list[i]:.6f}")
        print(f"  sigma_ref: {s_ref[i]:.6f}")
        print(f"  |<u_pt,u_ref>|: {abs(np.dot(U_pt_list[i], U_ref[:, i])):.6f}")
        print(f"  |<v_pt,v_ref>|: {abs(np.dot(V_pt_list[i], Vt_ref[i, :])):.6f}")

    X_MIN_Z, X_MAX_Z = 0.01, 2.0
    X_MIN_V, X_MAX_V = 0.01, 7.0

    start = perf_counter()

    ct_U_list, ct_S_list, ct_V_list, ct_A_deflated = poweriter.rank_k_svd_fhe(
        cryptocontext,
        key_pair.publicKey,
        cipher_a_svd,
        n=size,
        k=rank_k,
        num_iterations=fhe_iters,
        poly_degree=poly_deg,
        x_min_z=X_MIN_Z,
        x_max_z=X_MAX_Z,
        x_min_v=X_MIN_V,
        x_max_v=X_MAX_V,
        seed=42,
    )

    elapsed = perf_counter() - start

    print(f"\nFHE deflation complete in {elapsed:.2f}s\n")

    for i in range(rank_k):
        dec_u = cryptocontext.Decrypt(ct_U_list[i], key_pair.secretKey)
        dec_u.SetLength(size * size)

        u_fhe = poweriter.extract_vector_from_repeated_columns(
            dec_u.GetRealPackedValue(),
            size,
        )
        u_fhe = u_fhe / np.linalg.norm(u_fhe)

        dec_v = cryptocontext.Decrypt(ct_V_list[i], key_pair.secretKey)
        dec_v.SetLength(size * size)

        v_fhe = poweriter.extract_vector_from_repeated_columns(
            dec_v.GetRealPackedValue(),
            size,
        )
        v_fhe = v_fhe / np.linalg.norm(v_fhe)

        dec_s = cryptocontext.Decrypt(ct_S_list[i], key_pair.secretKey)
        dec_s.SetLength(1)
        sigma_fhe = dec_s.GetRealPackedValue()[0]

        print(f"Component {i}:")
        print(f"  sigma_fhe: {sigma_fhe:.6f}")
        print(f"  sigma_pt : {S_pt_list[i]:.6f}")
        print(f"  sigma_ref: {s_ref[i]:.6f}")
        print(f"  |<u_fhe,u_ref>|: {abs(np.dot(u_fhe, U_ref[:, i])):.6f}")
        print(f"  |<v_fhe,v_ref>|: {abs(np.dot(v_fhe, Vt_ref[i, :])):.6f}")
        print(f"  |<u_fhe,u_pt>| : {abs(np.dot(u_fhe, U_pt_list[i])):.6f}")
        print(f"  |<v_fhe,v_pt>| : {abs(np.dot(v_fhe, V_pt_list[i])):.6f}")
        print()

except Exception as e:
    print(f"Deflation SVD test failed: {e}")
