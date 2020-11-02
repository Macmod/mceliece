from pkenc_mceliece import McEliece
from secrets import token_bytes

if __name__ == '__main__':
    # Test key generation, encryption/decryption
    mcel = McEliece()
    (public_key, secret_key) = mcel.keygen(secparams=(6,64,4))

    plaintext = token_bytes(5)
    ciphertext = mcel.encrypt(public_key, plaintext)

    decrypted_msg = mcel.decrypt(public_key, secret_key, ciphertext["ct"])

    assert decrypted_msg == plaintext
