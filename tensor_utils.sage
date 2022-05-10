import numpy as np

def dyck_word_to_pairing(w):
    out = []
    for i in range(len(w)):
        j = w.associated_parenthesis(i)
        if j > i:
            out.append((i,j))
    return out

def pairing_to_dyck_word(L):
    w = [None]*(len(L)*2)
    for i,j in L:
        w[i] = 1
        w[j] = 0
    w = DyckWord(w)
    assert dyck_word_to_pairing(w) == L,'diagram not uncrossing'
    return w

def pairing_to_tensor(L):
    n = len(L)*2
    dimpermute = [k for ij in L for k in ij]
    dimpermute = [dimpermute.index(i) for i in range(n)]
    arr = np.zeros((2,2),dtype='object')
    arr[0,1] = 1
    arr[1,0] = -1
    T = reduce(lambda a,b: np.tensordot(a,b,0),[arr]*(n//2))
    return T.transpose(dimpermute)

# currently does not check sign of T
def tensor_to_pairing(T):
    n = len(T.shape)
    supp = T.nonzero()
    pairs = []
    check = set(range(n))
    while len(check) > 0:
        i = check.pop()
        j = next(j for j in range(n) if (supp[j] == 1-supp[i]).all())
        if j < i:
            i,j = j,i
        check.remove(j)
        pairs.append((i,j))
    pairs.sort()
    return pairs

def basis_vectors(n):
    basis = {}
    for w in DyckWords(n//2):
        basis[w] = pairing_to_tensor(dyck_word_to_pairing(w))
    return basis

def as_linear_combination(T):
    B = basis_vectors(len(T.shape)).items()
    M = matrix([list(S.ravel()) for _,S in B]).T
    return [(e,B[i][0]) for i,e in M.solve_right(vector(T.ravel())).dict().items()]

# vim: ft=python
