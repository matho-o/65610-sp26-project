from openfhe import *
from matrix import matrix_multiply, transpose

if __name__ == 'main':
    params = CCParamsCKKSRNS()

    params.SetSecretKeyDist(SecretKeyDist.UNIFORM_TERNARY)
    params.SetSecurityLevel(SecurityLevel.HEStd_128_classic)
    params.SetRingDim(512)

    params.SetScalingTechnique(ScalingTechnique.FLEXIBLEAUTO)
    params.SetScalingModSize(59)
    params.SetFirstModSize(60)
    params.SetMultiplicativeDepth(32)

    cc = GenCryptoContext(params)
    cc.Enable(PKESchemeFeature.PKE)
    cc.Enable(PKESchemeFeature.KEYSWITCH)
    cc.Enable(PKESchemeFeature.LEVELEDSHE)
    cc.Enable(PKESchemeFeature.ADVANCEDSHE)
    cc.Enable(PKESchemeFeature.FHE)

    level_budget = [4, 4]
    cc.EvalBootstrapSetup(level_budget)

    num_slots = cc.GetRingDimension() // 2

    rot_indices = [-16, -8, -4, -2, -1, 1, 2, 4, 8, 16]

    key_pair = cc.KeyGen()
    print("key pair generated")
    
    cc.EvalMultKeyGen(key_pair.secretKey)
    print("mult key generated")
    
    cc.EvalBootstrapKeyGen(key_pair.secretKey, num_slots)
    print("bootstrap key generated")
    
    cc.EvalRotateKeyGen(key_pair.secretKey, rot_indices)
    print("rotate key generated")

    ## test

    n = 16

    slot_size = cc.GetRingDimension()//2
    # can precompute this before any matrix multiplications if n known
    col = [1 if j%n == 0 else 0 for j in range(slot_size)]
    ct_col = cc.Encrypt(key_pair.publicKey, cc.MakeCKKSPackedPlaintext(col))
    row = [1 if j%(n*n)<n else 0 for j in range(slot_size)]
    ct_row = cc.Encrypt(key_pair.publicKey, cc.MakeCKKSPackedPlaintext(row))

    print(matrix_multiply(ct_col, ct_row, n, cc, key_pair))

    pass