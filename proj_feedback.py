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


#                     -2.046e-60 s^6 + 1.583e-50 s^5 + 2.535e-43 s^4 + 1.232e-36 s^3 + 1.636e-30 s^2
#   --------------------------------------------------------------------------------------------------------------------------------------
#   2.372e-74 s^8 + 1.191e-65 s^7 + 1.173e-57 s^6 + 1.689e-50 s^5 + 8.001e-44 s^4 + 1.075e-37 s^3 + 4.236e-33 s^2 + 4.663e-30 s + 4.02e-28


s = ct.TransferFunction.s

zeros = -2.046e-60*s**6 + 1.583e-50*s**5 + 2.535e-43*s**4 + 1.232e-36*s**3 + 1.636e-30*s**2
polos = 2.372e-74*s**8 + 1.191e-65*s**7 + 1.173e-57*s**6 + 1.689e-50*s**5 + 8.001e-44*s**4 + 1.075e-37*s**3 + 4.236e-33*s**2 + 4.663e-30*s + 4.02e-28


G = zeros/polos

def main():
    sistema = ct.TransferFunction(G)

    print("Funcao de transferencia:")
    print(sistema)
    
    plt.figure(figsize=(11, 7))
    frequencias_hz = np.logspace(0, 8, 1000000)
    ct.bode_plot(
        sistema,
        2*np.pi*frequencias_hz,
        label="Av Total",
        dB=True,
        deg=True,
        Hz=True,
        grid=True,
        display_margins="overlay"
        )

    plt.show()


if __name__ == "__main__":
    main()
