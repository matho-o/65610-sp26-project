from openfhe import *
from collections import deque

# This is a sample, should be the ones created in crypotcontext
rotation_keys = [-16, -8, -4, -2, -1, 1, 2, 4, 8, 16]
optimal_rotation = None

def optimize_rotation(mx):
    assert 1 in rotation_keys or -1 in rotation_keys, "impossible to rotate by 1"
    dp = [None for _ in range(mx)]
    global optimal_rotation
    optimal_rotation = [None for _ in range(mx)]
    queue = deque([0])
    dp[0] = 0
    while len(queue)!= 0:
        i = queue.popleft()
        for rot in rotation_keys:
            nxt_i = (i+rot)%mx
            if dp[nxt_i] is None:
                dp[nxt_i] = dp[i] + 1
                optimal_rotation[nxt_i] = rot
                queue.append(nxt_i)


def rotate(cc:CryptoContext, A:Ciphertext, ind:int):
    slot_size = cc.GetRingDimension()//2
    if optimal_rotation is None:
        optimize_rotation(slot_size)
    ind = ind % slot_size
    result = A
    while ind != 0:
        result = cc.EvalAtIndex(result, optimal_rotation[ind])
        ind = (ind - optimal_rotation[ind])%slot_size
    return result


def matrix_multiply(A:Ciphertext, B:Ciphertext, n:int, cc:CryptoContext, key_pair: KeyPair):
    """
    A is n by n, B is n by n
    They should be in row order, and periodically repeating to fill the cipher text
    This assumes n^2 divides the number of slots, and is a power of 2
    """
    slot_size = cc.GetRingDimension()//2
    # can precompute this before any matrix multiplications if n known
    col = [1 if j%n == 0 else 0 for j in range(slot_size)]
    ct_col = cc.Encrypt(key_pair.publicKey, cc.MakeCKKSPackedPlaintext(col))
    row = [1 if j%(n*n)<n else 0 for j in range(slot_size)]
    ct_row = cc.Encrypt(key_pair.publicKey, cc.MakeCKKSPackedPlaintext(row))
    result = None
    log2 = (n).bit_length() - 1
    for i in range(n):
        # column mask at i
        rot_A = rotate(cc, A, i)
        masked_A = cc.EvalMult(ct_col, rot_A)
        for j in range(log2):
            cc.EvalAddInPlace(masked_A, rotate(cc, masked_A, -pow(2, j)))
        # row mask at i
        rot_B = rotate(cc, B, n*i)
        masked_B = cc.EvalMult(ct_row, rot_B)
        for j in range(log2):
            cc.EvalAddInPlace(masked_B, rotate(cc, masked_B, -n * pow(2, j)))
        if result is None:
            result = cc.EvalMult(masked_A, masked_B)
        else:
            cc.EvalAddInPlace(result, cc.EvalMult(masked_A, masked_B))
    return result

def transpose(A: Ciphertext, n:int, cc:CryptoContext, key_pair: KeyPair):
    """
    A is n by n
    It should be in row order, and periodically repeating to fill the cipher text
    This assumes n^2 divides the number of slots
    """
    slot_size = cc.GetRingDimension()//2
    result = None
    for i in range(-n+1, n):
        pt_ti = [0 for _ in range(n*n)]
        for j in range(n - abs(i)):
            pt_ti[((n+1)*j+(i if i > 0 else abs(i)*n))%(n*n)] = 1
        pt_ti = pt_ti*(slot_size//(n*n))
        ti = cc.Encrypt(key_pair.publicKey, cc.MakeCKKSPackedPlaintext(pt_ti))
        if result is None:
            result = cc.EvalMult(ti, rotate(cc, A, (n-1)*i))
        else:
            cc.EvalAddInPlace(result, cc.EvalMult(ti, rotate(cc, A, (n-1)*i)))
    return result  