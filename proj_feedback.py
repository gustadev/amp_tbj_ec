import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import control as ct
import matplotlib.pyplot as plt
import numpy as np

# Primeiro Estagio
#           2.313e-45 s^2 + 5.562e-36 s - 1.821e-25
#   -------------------------------------------------------
#   7.772e-54 s^3 + 8.702e-44 s^2 + 8.219e-35 s + 1.211e-27
# Segundo Estagio
#           3.704e-90 s^5 - 4.3e-80 s^4 - 8.725e-73 s^3 - 3.559e-66 s^2 - 4.824e-60 s - 1.672e-54
#   ------------------------------------------------------------------------------------------------------
#   1.651e-94 s^6 + 6.096e-87 s^5 - 1.432e-78 s^4 - 3.018e-71 s^3 - 1.237e-64 s^2 - 1.68e-58 s - 5.828e-53


s = ct.TransferFunction.s

G1 = (2.313e-45*s*s + 5.562e-36*s - 1.821e-25)/(7.772e-54*s*s*s + 8.702e-44*s*s + 8.219e-35*s + 1.211e-27)
#G2 = (8.038e-31*s*s + 6.666e-2*s - 2.16e-13)/(1.178e-28*s*s + 9.395e-20*s + 9.922e-14)
G2 = (
    3.704e-90*s**5 - 4.3e-80*s**4 - 8.725e-73*s**3 - 3.559e-66*s**2 
    - 4.824e-60*s - 1.672e-543
) / (
    1.651e-94*s**6 + 6.096e-87*s**5 - 1.432e-78*s**4 
    - 3.018e-71*s**3 - 1.237e-64*s**2 - 1.68e-58*s - 5.828e-53
)
G = G1*G2

def main():
    sistema = ct.TransferFunction(G2)

    print("Funcao de transferencia:")
    print(sistema)

    ct.bode_plot(
        sistema,
       dB=True,
       deg=True,
       grid=True,
    )

    plt.show()


if __name__ == "__main__":
    main()
