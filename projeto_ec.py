#%%
import numpy as np


def prl(a):
    rt = 0 
    for i in a:
        rt += 1/i
    return 1/rt


def projeto_EC(
    Vss=12,
    hfe=150,
    Rs=50,
    Vce=6,
    Ic=1e-3,
    Ve=0.4,
    fator_divisor=10,
    fL = 150,
    Rl = 3.3e3,
    fT = 120e6,
    Cmu = 5e-12,
    fH = 170e3
):

    # Constantes
    VT = 25.8e-3
    Vbe = 0.65

    # Correntes
    Ib = Ic / hfe
    Ie = Ic + Ib

    # Pequenos sinais
    gm = Ic / VT
    hie = hfe / gm      # = beta * re
    re = 1 / gm

    # Resistores de polarização
    Re = Ve / Ie
    Rc = (Vss - Vce - Ve) / Ic

    # Tensão da base
    Vb = Ve + Vbe

    # Dimensionamento do divisor
    IR2 = fator_divisor * Ib
    IR1 = IR2 + Ib

    R2 = Vb / IR2
    R1 = (Vss - Vb) / IR1
    Rb = prl([R1,R2])

    # Ganho aproximado (Re bypassado)
    Av = -gm * Rc

    Zin = prl([Rb,hie])
    #Rce = prl([Re,(hie+prl([Rs,R1,R2])/hfe)])#prl([(Rb+hie)/hfe,Re])
    Rce = prl([
        Re,
        hie/hfe + prl([Rs,Rb])/(hfe+1)
    ])
    Rci = Zin + Rs
    Rco = Rc + Rl

    Ci = 1/(2*np.pi*Rci*fL/10)
    Co = 1/(2*np.pi*Rco*fL/10)
    Ce = 1/(2*np.pi*Rce*fL)

    #Alta frequencia

    #---------------------------------
    # Miller
    #---------------------------------
    Ctotal = gm / (2*np.pi*fT)
    Cpi = Ctotal - Cmu
    Cin = Cpi + Cmu*(1-Av)
    Cout = Cmu*(1-1/Av)
    Rin = prl([Rs,R1,R2,hie/hfe])
    Rout = prl([Rc,Rl])
    fH_in = 1/(2*np.pi*Rin*Cin)
    fH_out = 1/(2*np.pi*Rout*Cout)
    CL = 1/(2*np.pi*Rout*fH)


    print("--------------------------------")
    print("Projeto DC")
    print("--------------------------------")
    print(f"Ic   = {Ic*1e3:.2f}mA")
    print(f"Ib   = {Ib*1e6:.2f}uA")
    print(f"Ve   = {Ve:.3f}V")
    print(f"Vb   = {Vb:.3f}V")
    print(f"Vce  = {Vce:.3f}V")
    print()
    print(f"Rc   = {Rc:.1f}Ω")
    print(f"Re   = {Re:.1f}Ω")
    print(f"R1   = {R1:.1f}Ω")
    print(f"R2   = {R2:.1f}Ω")
    print(f"Rb   = {Rb:.1f}Ω")
    print()
    print(f"IR2  = {IR2*1e6:.2f}uA")
    print(f"IR1  = {IR1*1e6:.2f}uA")
    print()
    print(f"gm   = {gm:.4f} S")
    print(f"re   = {re:.2f}Ω")
    print(f"hie  = {hie:.2f}Ω")
    print(f"Av(Max)   = {Av:.2f}V/V")
    print()
    print(f"Zin = {Zin:.3e}Ω")
    print(f"Rce = {Rce:.3e}Ω")
    print(f"Rci = {Rci:.3e}Ω")
    print(f"Rco = {Rco:.3e}Ω")
    print()
    print(f"Ci = {Ci:.3e}F")
    print(f"Co = {Co:.3e}F")
    print(f"Ce = {Ce:.3e}F")
    print()
    print(f"Cpi   = {Cpi*1e12:.2f}pF")
    print(f"Cmu   = {Cmu*1e12:.2f}pF")
    print(f"Cin   = {Cin*1e12:.2f}pF")
    print(f"Cout  = {Cout*1e12:.2f}pF")
    print()
    print(f"RCin   = {Rin:.2f}Ω")
    print(f"RCout  = {Rout:.2f}Ω")
    print()
    print(f"CL = {CL:.3e}F")
    print()
    print(f"fH_in  = {fH_in*1e-6:.2f}MHz")
    print(f"fH_out = {fH_out*1e-6:.2f}MHz")
    print()
    print(f"fH ≈ {min(fH_in,fH_out)*1e-6:.2f}MHz")
    print("--------------------------------")

#%%
if __name__ == "__main__":
    print("----- Primeiro Estagio ------")
    projeto_EC(
        Ic=1e-3,
        Ve=0.5,
        Vce=6
    )
    print("----- Segundo Estagio ------")
    projeto_EC(
        Ic=10e-3,
        Ve=3,
        Vce=7
    )

# %%
