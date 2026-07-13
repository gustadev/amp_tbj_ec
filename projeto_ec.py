#%%
import numpy as np
import control as ct
import matplotlib.pyplot as plt


def formatar_frequencia(freq_rad_s):
    if not np.isfinite(freq_rad_s):
        return "infinita"

    freq_hz = freq_rad_s / (2 * np.pi)

    if freq_hz < 1e3:
        return f"{freq_hz:.4f} Hz"
    if freq_hz < 1e6:
        return f"{freq_hz * 1e-3:.4f} kHz"
    if freq_hz < 1e9:
        return f"{freq_hz * 1e-6:.4f} MHz"
    return f"{freq_hz * 1e-9:.4f} GHz"


def imprimir_margens_estabilidade(sistema, method="best"):
    gm, pm, sm, wpc, wgc, wms = ct.stability_margins(sistema, method=method)

    gm_db = 20 * np.log10(gm) if np.isfinite(gm) and gm > 0 else np.inf

    print()
    print("Margens de estabilidade")
    print("--------------------------------")
    if np.isfinite(gm):
        print(f"Margem de ganho (GM): {gm:.4f} V/V ({gm_db:.2f} dB)")
    else:
        print("Margem de ganho (GM): infinita")

    if np.isfinite(pm):
        print(f"Margem de fase (PM): {pm:.2f}°")
    else:
        print("Margem de fase (PM): infinita")

    if np.isfinite(sm):
        print(f"Margem de estabilidade (SM): {sm:.4f}")
    else:
        print("Margem de estabilidade (SM): infinita")

    print(f"Fase = -180° em: {formatar_frequencia(wpc)}")
    print(f"|G(jw)| = 1 (0 dB) em: {formatar_frequencia(wgc)}")
    print(f"Nyquist mais proximo de -1 em: {formatar_frequencia(wms)}")
    print("--------------------------------")

def prl(elementos):
    rt = 0

    for z in elementos:
        if z is None:
            continue
        rt += 1/z

    return 1/rt

def projeto_EC(
    Vss=12,
    hfe=160,
    Rs=50,
    Vce=6,
    Ic=1e-3,
    Ve=0.4,
    fator_divisor=10,
    fL = 150,
    Rl = 3.3e3,
    fT = 170e6,
    Cmu = 5e-12,
    fH = 170e3,
    Re_dc = 0
):

    # Constantes
    VT = 25.8e-3
    Vbe = 0.65

    # Correntes
    Ib = Ic / hfe
    Ie = Ic + Ib

    # Pequenos sinais
    gm = Ic / VT
    hie = hfe / gm
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
    Zin = prl([
        Rb,
        hie + (hfe+1)*Re_dc
    ])
    Av = -hfe * prl([Rc,Rl]) / (hie + (hfe + 1)*Re_dc)

    Rce = prl([
        Re-Re_dc,
        hie/hfe + Re_dc + prl([Rs,Rb])/(hfe+1)
    ])
    Rci = Zin + Rs
    Rco = Rc + Rl

    Ci = 1/(2*np.pi*Rci*(fL/10))
    Co = 1/(2*np.pi*Rco*(fL/10))
    Ce = 1/(2*np.pi*Rce*(fL))

    #Alta frequencia

    #---------------------------------
    # Miller
    #---------------------------------
    Ctotal = gm / (2*np.pi*fT)
    Cpi = Ctotal - Cmu
    Cin = Cpi + Cmu*(1-Av)
    Cout = Cmu*(1-1/Av)
    Rin = prl([
        Rs,
        R1,
        R2,
        hie + (hfe+1)*Re_dc
    ])
    Rout = prl([Rc,Rl])
    fH_in = 1/(2*np.pi*Rin*Cin)
    fH_out = 1/(2*np.pi*Rout*Cout)
    CL = 1/(2*np.pi*Rout*fH)

    #-----------------------------------
    # Modelo Altas Frequencias
    #------------------------------------

    s = ct.TransferFunction.s
    rx = 25
    gx = 1/rx
    Gs = 1/Rs
    Gb = 1/Rb
    RL = prl([Rc,Rl,1/(s*CL)])

    if Re_dc > 0:
        a = prl([hie,1/(s*Cpi)])/(prl([hie,1/(s*Cpi)])+Re_dc*hfe)
        b = hfe/prl([hie,1/(s*Cpi)])
        Avs = -a*(s*Cmu + b)/(s*Cmu + 1/Rc + 1/Rl)
        Avs = Avs * (s**2 / ((s + 1/(Rce*Ce) + 1/(Rci*Ci))*(s + 1/(Rco*Co))))
    else:
        Ra = prl([RL,hie,(rx+prl([Rs,Rb]))])
        b = (1/(RL*Cmu)) + ((1 + gm*Ra)/(Ra*Cpi))
        c = 1/(prl([hie,(rx + prl([Rs,Rb]))])*RL*Cmu*Cpi)
        Avs = (-gx*Gs/((gx+Gs+Gb)*Cmu*Cpi)) * ((gm - s*Cmu)/(s*s + b*s + c))
        Avs = Avs * (s**2 / ((s + 1/(Rce*Ce) + 1/(Rci*Ci))*(s + 1/(Rco*Co))))

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
    if Re_dc > 0:
        print(f"Re(AC)= {Re-Re_dc:.1f}Ω")
        print(f"Re(DC)= {Re_dc:.1f}Ω")
    else :   
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
    print(f"Av   = {Av:.2f}V/V")
    print()
    print(f"Zin = {Zin:.3e}Ω")
    print(f"Rce = {Rce:.3e}Ω")
    print(f"Rci = {Rci:.3e}Ω")
    print(f"Rco = {Rco:.3e}Ω")
    print()
    print(f"Ci = {Ci*1e6:.3e}uF")
    print(f"Co = {Co*1e6:.3e}uF")
    print(f"Ce = {Ce*1e6:.3e}uF")
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

    # print("Funcao de transferencia:")
    # print(Avs)
    # print()

    return Avs

#%%
if __name__ == "__main__":
    # Projeto funcionando sem realimentação
    # print("----- Primeiro Estagio ------")
    # Avs1 = projeto_EC(
    #     Ic=1e-3,
    #     Ve=0.5,
    #     Vce=6,
    #     Rl = 1.383e4,
    #     hfe=170,
    #     fH = 6.3e3#20e6
    # )
    # print("----- Segundo Estagio ------")
    # Avs2 = projeto_EC(
    #     Ic=1.5e-3,
    #     Ve=1,
    #     Vce=5,
    #     Rs = 0.000000001,
    #     fH= 6.3e3,#6.3e3,
    #     Re_dc=600,
    #     hfe=170
    # )
    print("----- Primeiro Estagio ------")
    # Avs1 = projeto_EC(
    #     Ic=1e-3,
    #     Ve=0.51,
    #     Vce=1.5,
    #     Rl = 2470,
    #     Rs = 50,
    #     hfe=170,
    #     Re_dc = 500,
    #     fH = 6.3e3,#20e6
    #     fL=100
    # )
    Avs1 = projeto_EC(
        Ic=10e-3,
        Ve=0.750,
        Vce=6,
        Rl=2470,
        Rs=50,
        hfe=170,
        Re_dc=70,
        fH=5e3,
        fL=15
    )
    print("----- Segundo Estagio ------")
    Avs2 = projeto_EC(
        Ic=1.5e-3,
        Ve=1,
        Vce=5,
        Rs = 525,#0.000000001,
        fH= 6.3e3,#6.3e3,
        fL=150,
        Re_dc=0,
        hfe=170
    )
    s = ct.TransferFunction.s

    AvsT = Avs1 * Avs2 #* AvLP

    sistema = ct.TransferFunction(AvsT)
    FT1 = ct.TransferFunction(Avs1)
    FT2 = ct.TransferFunction(Avs2)
    

    # print()
    # print("Funcao de transferencia:")
    # print(sistema)

    print()
    # polos = sistema.poles()
    # zeros = sistema.zeros()
    # print(f"Polos = {polos}")
    # print("Lista de polos:")
    # for i, polo in enumerate(polos, start=1):
    #     modulo = np.abs(polo)
    #     angulo_deg = np.angle(polo, deg=True)
    #     frequencia_hz = modulo / (2*np.pi)
    #     if frequencia_hz < 1e3 :
    #         print(f"p{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz:.4e} Hz")
    #     elif frequencia_hz < 1e6 :
    #         print(f"p{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz*1e-3:.4e} kHz")
    #     elif frequencia_hz < 1e9 :
    #         print(f"p{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz*1e-6:.4e} MHz")
    
    # for i, zero in enumerate(zeros, start=1):
    #     modulo = np.abs(zero)
    #     angulo_deg = np.angle(zero, deg=True)
    #     frequencia_hz = modulo / (2*np.pi)
    #     if frequencia_hz < 1e3 :
    #         print(f"z{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz:.4e} Hz")
    #     elif frequencia_hz < 1e6 :
    #         print(f"z{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz*1e-3:.4e} kHz")
    #     elif frequencia_hz < 1e9 :
    #         print(f"z{i} : {modulo:.4e} ∠ {angulo_deg:.2f}° : {frequencia_hz*1e-6:.4e} MHz")

    plt.figure(figsize=(11, 7))

    imprimir_margens_estabilidade(sistema, method="best")

    frequencias_hz = np.logspace(0, 8, 1000000)

    ct.bode_plot(
        sistema,
        2*np.pi*frequencias_hz,
        dB= True,
        deg=True,
        Hz=True,
        grid=True
    )

    plt.show()
# %%
