from openfhe import *
import matrix
import os
import numpy as np
import random

# --- YOUR ORIGINAL PARAMETER SETUP ---
parameters = CCParamsCKKSRNS()

secret_key_dist = SecretKeyDist.UNIFORM_TERNARY
parameters.SetSecretKeyDist(secret_key_dist)

parameters.SetSecurityLevel(SecurityLevel.HEStd_NotSet)
parameters.SetRingDim(1 << 9)

rescale_tech = ScalingTechnique.FLEXIBLEAUTO
dcrt_bits = 29
first_mod = 30

parameters.SetScalingModSize(dcrt_bits)
parameters.SetScalingTechnique(rescale_tech)
parameters.SetFirstModSize(first_mod)

level_budget = [4, 4]
levels_available_after_bootstrap = 3
depth = levels_available_after_bootstrap + FHECKKSRNS.GetBootstrapDepth(level_budget, secret_key_dist)

parameters.SetMultiplicativeDepth(depth)

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
cryptocontext.EvalBootstrapKeyGen(key_pair.secretKey, num_slots)

# --- ROTATION KEY GENERATION (Using your requested range) ---
rot_indices = [-16, -8, -4, -2, -1, 1, 2, 4, 8, 16]
cryptocontext.EvalAtIndexKeyGen(key_pair.secretKey, rot_indices)
print(f"Rotation keys generated for: {rot_indices}")

# --- FULL TEST SUITE ---

# 1. Prepare Test Data
size = 16
matrix_a_data = [random.random() for i in range(size*size)]
matrix_b_data = [random.random() for i in range(size*size)]

# 2. Encrypt
cipher_a = cryptocontext.Encrypt(key_pair.publicKey, cryptocontext.MakeCKKSPackedPlaintext(matrix_a_data*(num_slots//(size*size))))
cipher_b = cryptocontext.Encrypt(key_pair.publicKey, cryptocontext.MakeCKKSPackedPlaintext(matrix_b_data*(num_slots//(size*size))))

# 4. Matrix Multiplication Test
print("\n--- Step 2: Testing matrix.matrix_multiply ---")
try:
    # Format: (ct, ct, cc, key_pair, size)
    res_mult = matrix.matrix_multiply(cipher_a, cipher_b, size, cryptocontext, key_pair)
    
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
    res_trans = matrix.transpose(cipher_a, size, cryptocontext, key_pair)
    
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