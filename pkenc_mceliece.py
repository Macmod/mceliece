'''
McEliece Public Key Encryption Scheme

| Available from: https://en.wikipedia.org/wiki/McEliece_cryptosystem

* type:          encryption (public key)

:Authors: Artur Henrique Marzano Gonzaga
:Date: 15/10/2020
'''

from charm.toolbox.PKEnc import PKEnc
from sage.all import *
from sage.coding.goppa_code import GoppaCode, GoppaCodeEncoder
from numpy.random import choice, permutation
from itertools import combinations

class McElieceCipher(dict):
    def __init__(self, ct):
        dict.__init__(self, {'ct': ct})

    def __repr__(self):
        return hex(int(''.join(list(map(str, list(self['ct'][0])))), 2))[2:]

    def __str__(self):
        return hex(int(''.join(list(map(str, list(self['ct'][0])))), 2))[2:]

class McEliece(PKEnc):
    def __init__(self):
        PKEnc.__init__(self)
        self.k = None

    def _patterson_decode(self, MSG_e, sk):
        print('(Patterson decoding)')
        print()

        t = sk['g'].degree()
        X = sk['g'].parent().gen()

        # Get syndrome polynomial
        syndrome_poly = sum(list(
            (
                (X - sk['alphas'][i]).inverse_mod(sk['g'])
            ) * MSG_e[0,i] for i in range(MSG_e.ncols())
        ))

        print('Syndrome Polynomial:')
        print(syndrome_poly)
        print()

        s_inv_g = syndrome_poly.inverse_mod(sk['g'])
        print('s_inv_g:')
        print(s_inv_g)
        print()

        # Calculate sqrt(s_inv_g - X) mod g
        h = s_inv_g - X
        h_list = h.list()
        g_list = sk['g'].list()

        gs_even = self.PR(list(map(sqrt,
            [g_list[x] for x in range(len(g_list)) if x % 2 == 0]
        )))
        gs_odd = self.PR(list(map(sqrt,
            [g_list[x] for x in range(len(g_list)) if x % 2 != 0]
        )))

        hs_even = self.PR(list(map(sqrt,
            [h_list[x] for x in range(len(h_list)) if x % 2 == 0]
        )))
        hs_odd = self.PR(list(map(sqrt,
            [h_list[x] for x in range(len(h_list)) if x % 2 != 0]
        )))

        gs_odd_inv = gs_odd.inverse_mod(sk['g'])

        h_root = (hs_even + hs_odd * gs_even * gs_odd_inv)
        assert (h_root**2).mod(sk['g']) == h.mod(sk['g'])

        # Basis reduce
        (q, r) = sk['g'].quo_rem(h_root)
        z = simplify((sk['g'] - q*h_root, 0 - q))
        a = [z[0]]
        b = [z[1]]

        i = 0
        while ((a[i]**2 + X*b[i]**2).degree()) > t:
            f = a[i-1] if i > 0 else h_root
            k = b[i-1] if i > 0 else 1
            (q, r) = f.quo_rem(a[i])
            a.append(r)
            b.append(k - q*b[i])

            i += 1

        # Decode resulting code
        error_locator = a[i]**2 + X*(b[i]**2)
        print('Error Locator:')
        print(error_locator)
        print()

        error = [
            i for i in range(len(sk['alphas']))
                if error_locator(sk['alphas'][i]) == 0
        ]

        print('Error bits (after Pinv):')
        print(error)
        print()

        print('Error bits (before Pinv):')
        error_vec = matrix(GF(2), [(1 if i in error else 0) for i in range(len(sk['alphas']))])
        error_vec_fixed = error_vec * sk['P']
        error_fixed = [i for i in range(len(sk['alphas'])) if error_vec_fixed[0,i] == 1]
        print(error_fixed)
        print()

        MSG = copy(MSG_e)
        for ebit in error:
            MSG[0,ebit] = MSG[0,ebit]+1

        MS = sk['G'].solve_left(MSG)

        return MS

    def _rand_poly(self, PR, t):
        poly = PR.random_element(t)
        poly_list = poly.list()
        poly_list[-1] = 1

        return PR(poly_list)

    """ Generates parameters for McEliece cryptosystem """
    def keygen(self, secparams=(6,64,4)):
        print('============ KEY GENERATION ===========')
        print('m, n, t = %s' % str(secparams))
        print()
        pk = {}
        sk = {}

        m, n, t = secparams
        if n > 2**m:
            raise Exception("Invalid parameters (n should be <= 2**m)")

        F2m = GF(2**m)
        PR = PolynomialRing(F2m, 'X')
        self.PR = PR

        # Get irreducible polynomial of degree t
        # I could also explicitly use rabin test for this step
        g = self._rand_poly(PR, t)
        while g.degree() != t or not g.is_irreducible():
            g = self._rand_poly(PR, t)

        sk['g'] = g

        print('Polynomial:')
        print(g)
        print()

        # Get random code locators
        alphas = []
        while len(alphas) < n:
            elem = F2m.random_element()
            if g(elem) != 0 and elem not in alphas:
                alphas.append(elem)

        sk['alphas'] = alphas

        print('Alphas (code locators):')
        print(alphas)
        print()

        # G: Code's generator matrix
        # H: Code's parity-check matrix
        # g_poly: Polynomial of degree t that defines the code
        goppa_obj = GoppaCode(g, alphas)
        sk["H"] = goppa_obj.parity_check_matrix()

        encoder = GoppaCodeEncoder(goppa_obj)
        sk["G"] = encoder.generator_matrix()

        k = self.k = sk["G"].nrows()
        print('k = %d' % k)

        # P: Random permutation matrix
        # S: Random invertible matrix
        # G~: S*G*P (public code)

        P = permutation(list(range(1,n+1)))
        sk["P"] = matrix(GF(2), Permutation(P).to_matrix())

        sk["S"] = random_matrix(GF(2), k, k)
        while det(sk["S"]) != 1:
            sk["S"] = random_matrix(GF(2), k, k)

        pk["G~"] = (sk["S"] * sk["G"]) * sk["P"]

        print("S Matrix: %d x %d" % (sk["S"].nrows(), sk["S"].ncols()))
        print("G Matrix: %d x %d" % (sk["G"].nrows(), sk["G"].ncols()))
        print("P Matrix: %d x %d" % (sk["P"].nrows(), sk["P"].ncols()))
        print("G~ Matrix: %d x %d" % (pk["G~"].nrows(), pk["G~"].ncols()))

        pk["t"] = t

        return (pk, sk)
    
    def encrypt(self, pk, M):
        print('============ ENCRYPTION ===========')

        M_bin = ''.join(format(x, '08b') for x in M)
        M_bin = '0'*((8-len(M_bin)%8)%8) + M_bin
        M_bin = list(map(int, list(M_bin)))
        print('Plaintext:', M.hex())
        print(M_bin)
        print()

        assert len(M_bin) == self.k

        MG = matrix(GF(2), [M_bin]) * pk["G~"]

        error = list(choice(MG.ncols(), size=pk["t"], replace=False))
        print('Error bits:')
        print(error)
        print()

        for ebit in error:
            MG[0,ebit] = MG[0,ebit]+1
        ciphertext = McElieceCipher(MG)

        print('Ciphertext:', ciphertext)
        print(list(MG[0]))

        return ciphertext

    def decrypt(self, pk, sk, C):
        print('============ DECRYPTION ===========')
        if C.ncols() != sk["H"].ncols():
            raise Exception(f"Wrong message length.")

        # Remove effect of "P" from the code
        # This operation also permutates the error vector,
        # that is, simply transforms it into another error vector with the same weight.
        MSG_e = C * sk["P"].inverse()

        MS = self._patterson_decode(MSG_e, sk)

        M = MS * sk['S'].inverse()
        M = ''.join(list(map(str, list(M[0]))))

        M_bytes = bytes(int(M[i:i+8], 2) for i in range(0, len(M), 8))

        print('Decrypted message: %s' % M_bytes.hex())

        return M_bytes
