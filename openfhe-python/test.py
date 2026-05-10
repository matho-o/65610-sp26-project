from openfhe import *
import matrix
import subspace_iteration
import os
import numpy as np
import random

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
# levels_available_after_bootstrap = 3
# depth = levels_available_after_bootstrap + FHECKKSRNS.GetBootstrapDepth(level_budget, secret_key_dist)

parameters.SetMultiplicativeDepth(70)

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
# cryptocontext.EvalBootstrapSetup(level_budget)
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
random.seed(99)
np.random.seed(99)
matrix_a_data = [random.random() for i in range(size*size)]
matrix_b_data = [random.random() for i in range(size*size)]

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
        print("✅ Multiplication Success!")
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
        print("✅ Transpose Success!")
except Exception as e:
    print(f"Transpose test failed: {e}")

# 6. Subspace Iteration Test
print("\n--- Step 4: Testing subspace iteration ---")

# reshape to 2D for numpy
A = np.array(matrix_a_data).reshape(size, size)
A_np = np.array(matrix_a_data).reshape(size, size)
frob_norm = np.linalg.norm(A_np, 'fro')
A_normalized = A_np / frob_norm
matrix_a_data_normalized = A_normalized.flatten().tolist()

# encrypt the normalized version
cipher_a_svd = cryptocontext.Encrypt(key_pair.publicKey, 
    cryptocontext.MakeCKKSPackedPlaintext(
        matrix_a_data_normalized * (num_slots // (size * size))))

# ── plaintext reference ───────────────────────────────────────────────────────
fhe_iters = 1
k=2
print(f"\n  [Plaintext reference — {fhe_iters} iteration(s)]")
U_pt, s_pt, Vt_pt = subspace_iteration.subspace_iteration_plaintext(
    A, k, num_iterations=fhe_iters, seed=42)
U_ref, s_ref, Vt_ref = np.linalg.svd(A, full_matrices=False)

print(f"  {'':12} {'Plaintext SI':>14} {'Numpy SVD':>12} {'Rel error':>10}")
print(f"  {'-'*50}")
for i in range(k):
    rel_err = abs(s_pt[i] - s_ref[i]) / s_ref[i]
    print(f"  sigma_{i}        {s_pt[i]:>14.6f} {s_ref[i]:>12.6f} {rel_err:>9.4%}")

# ── norm diagnostic ───────────────────────────────────────────────────────────
print("\n  [Norm diagnostic]")
np.random.seed(42)
Q_check = np.random.randn(size, k)
Q_check, _ = np.linalg.qr(Q_check)

X_MIN_Z, X_MAX_Z = 0.05, 2.0
X_MIN_Q, X_MAX_Q = 0.05, 7.0

norm_ok = True
for j in range(k):
    z = A_normalized @ Q_check[:, j]
    nz = np.dot(z, z)
    z_unit = z / np.linalg.norm(z)
    q2 = A_normalized.T @ z_unit
    nq = np.dot(q2, q2)
    z_ok = X_MIN_Z <= nz <= X_MAX_Z
    q_ok = X_MIN_Q <= nq <= X_MAX_Q
    print(f"    col {j}: {'✅' if z_ok else '⚠️ '} ||A@q||²={nz:.4f} [{X_MIN_Z},{X_MAX_Z}]  "
          f"{'✅' if q_ok else '⚠️ '} ||Aᵀz||²={nq:.4f} [{X_MIN_Q},{X_MAX_Q}]")
    if not z_ok or not q_ok:
        norm_ok = False

if not norm_ok:
    print("\n  ⚠️  WARNING: norms outside range — accuracy may be degraded.\n")
else:
    print("\n  ✅ All norms in range — proceeding.\n")

# ── FHE run ───────────────────────────────────────────────────────────────────
import time

poly_deg = 7
print(f"  Config: size={size}, k={k}, iters={fhe_iters}, poly_degree={poly_deg}, GS=False")

ct_At = matrix.transpose(cipher_a_svd, size, cryptocontext, key_pair.publicKey)

t0 = time.time()
ct_U_cols, ct_s, ct_Vt_rows = subspace_iteration.subspace_iteration_fhe(
    cryptocontext, key_pair.publicKey,
    cipher_a_svd, ct_At,
    n=size, k=k,
    num_iterations=fhe_iters,
    poly_degree=poly_deg,
    x_min_z=X_MIN_Z, x_max_z=X_MAX_Z,
    x_min_q=X_MIN_Q, x_max_q=X_MAX_Q,
    seed=42, use_gram_schmidt=False
)
fhe_elapsed = time.time() - t0

# ── results ───────────────────────────────────────────────────────────────────
print(f"\n  FHE complete in {fhe_elapsed:.1f}s\n")
print(f"  {'':4} {'|<u_fhe,u_ref>|':>17} {'|<u_fhe,u_pt>|':>16} {'||u_fhe||':>11}  note")
print(f"  {'-'*65}")

dots_ref, dots_pt, norms_fhe = [], [], []
for i, ct_u in enumerate(ct_U_cols):
    dec = cryptocontext.Decrypt(ct_u, key_pair.secretKey)
    dec.SetLength(size)
    u_fhe = np.array(dec.GetRealPackedValue()[:size])
    norm = np.linalg.norm(u_fhe)
    u_fhe_n = u_fhe / norm
    dot_ref = abs(np.dot(u_fhe_n, U_ref[:, i]))
    dot_pt  = abs(np.dot(u_fhe_n, U_pt[:, i]))
    norms_fhe.append(norm)
    dots_ref.append(dot_ref)
    dots_pt.append(dot_pt)
    note = "✅ good" if dot_ref > 0.8 else ("⚠️  moderate" if dot_ref > 0.5 else "❌ poor")
    print(f"  u{i}  {dot_ref:>17.4f} {dot_pt:>16.4f} {norm:>11.4f}  {note}")

avg_ref = np.mean(dots_ref)
avg_pt  = np.mean(dots_pt)
print(f"\n  Error decomposition (avg over {k} vectors):")
print(f"    Total accuracy vs numpy SVD   : {avg_ref:.4f}  (ground truth)")
print(f"    Accuracy vs plaintext SI      : {avg_pt:.4f}  (same algorithm, no FHE noise)")
print(f"    Algorithmic error  (1-pt)     : {1-avg_pt:.4f}  (finite iterations)")
print(f"    FHE overhead  (pt-ref)        : {avg_pt-avg_ref:.4f}  (CKKS approximation)")
print(f"    Avg ||u_fhe||                 : {np.mean(norms_fhe):.4f}  (should be ~1.0)")
print(f"    Runtime                       : {fhe_elapsed:.1f}s")
