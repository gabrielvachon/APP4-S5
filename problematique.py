import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import scipy.signal
from scipy import signal
import zplane

# Définition variables
e = np.e
pi = np.pi

def aberration():
    img = np.load("goldhill_aberrations.npy")   # Chargement de l'image
    plt.gray()                                  # Image en niveau de gris

    # Numérateur de la transformée em Z
    H_z_num = np.poly([0.9 * e ** (1j * pi / 2), 0.9 * e ** (-1j * pi / 2), 0.95 * e ** (1j * pi / 8), 0.95 * e ** (-1j * pi / 8)])
    # Démumérateur de la transformée en Z
    H_z_denum = np.poly([0, -0.99, -0.99, 0.8])

    zplane.zplane(H_z_num, H_z_denum)   # Plan en Z de la fonction H(Z)
    zplane.zplane(H_z_denum, H_z_num)   # Plan en Z de la fonction H^-1(Z)

    # Obtention de la fonction en fréquence normalisée H(w)
    _, H_w = signal.freqz(H_z_num, H_z_denum)

    # Filtrage de l'image pour retirer les aberrations
    img_no_abberation = signal.lfilter(H_z_denum, H_z_num, img)

    # Affichage de l'amplitude et de la phase
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(np.log10(np.abs(H_w)), 'red')
    ax2.plot(np.unwrap(np.abs(H_w)), 'blue')
    ax1.set_xlabel('Freq [rad/ech.]')
    ax1.set_ylabel('Amplitude [dB]', color='red')
    ax2.set_ylabel('Phase [rad]', color='blue')
    plt.show()

    # Affichage de l'image filtrée
    plt.figure()
    plt.imshow(img)
    plt.figure()
    plt.imshow(img_no_abberation)
    plt.show()

def debruitage_comp_filtre(pFilter=None):
    # Définir spécifications
    Fe = 1600       # Fréquence d'échantillonnage
    wPass = 500     # Fréquence passante 500Hz
    wStop = 750     # Fréquence coupure à 750Hz
    gPass = 0.2     # Variation de gain max de la bande passante (dB)
    gStop = 60      # Maximum d'atténuation de la bande coupée (dB)

    if pFilter == 'b' or pFilter == None:
        # Obtenir le filtre qui respecte les critères spécifiés
        N, _ = signal.buttord(wPass, wStop, gPass, gStop, fs=Fe)
        print("Ordre du filtre Butterworth : " + str(N))

        # Conception du filtre Butterworth
        num, denum = signal.butter(N, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

        # Chercher la réponse en fréquence du filtre
        w, H_w = signal.freqz(num, denum)
        frequencies = w * Fe / (2 * np.pi)  # Conversion de rad/s en Hz

        # Afficher, en fréquence, le module et la phase du filtre
        fig, ax1 = plt.subplots()
        ax1.plot(frequencies, 20 * np.log10(np.abs(H_w)), 'red')
        ax1.set_xlabel('Freq [Hz]')
        ax1.set_ylabel('Amplitude [dB]', color='red')
        plt.show()

        # Afficher le zplane
        zplane.zplane(num, denum)

    if pFilter == 'c1' or pFilter == None:
        # Obtenir le filtre qui respecte les critères spécifiés
        N, _ = signal.cheb1ord(wPass, wStop, gPass, gStop, fs=Fe)
        print("Ordre du filtre Chebyshev de type 1 : " + str(N))

        # Conception du filtre Chebyshev type 1
        num, denum = signal.cheby1(N, rp=gPass, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

        # Chercher la réponse en fréquence du filtre
        _, H_w = signal.freqz(num, denum)
        frequencies = w * Fe / (2 * np.pi)  # Conversion de rad/s en Hz

        # Afficher, en fréquence, le module et la phase du filtre
        fig, ax1 = plt.subplots()
        ax1.plot(frequencies, 20 * np.log10(np.abs(H_w)), 'red')
        ax1.set_xlabel('Freq [Hz]')
        ax1.set_ylabel('Amplitude [dB]', color='red')
        plt.show()

        # Afficher le zplane
        zplane.zplane(num, denum)

    if pFilter == 'c2' or pFilter == None:
        # Obtenir le filtre qui respecte les critères spécifiés
        N, _ = signal.cheb2ord(wPass, wStop, gPass, gStop, fs=Fe)
        print("Ordre du filtre Chebyshev de type 2 : " + str(N))

        # Conception du filtre Chebyshev type 2
        num, denum = signal.cheby2(N, rs=gPass, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

        # Chercher la réponse en fréquence du filtre
        _, H_w = signal.freqz(num, denum)
        frequencies = w * Fe / (2 * np.pi)  # Conversion de rad/s en Hz

        # Afficher, en fréquence, le module et la phase du filtre
        fig, ax1 = plt.subplots()
        ax1.plot(frequencies, 20 * np.log10(np.abs(H_w)), 'red')
        ax1.set_xlabel('Freq [Hz]')
        ax1.set_ylabel('Amplitude [dB]', color='red')
        plt.show()

        # Afficher le zplane
        zplane.zplane(num, denum)

    if pFilter == 'e' or pFilter == None:
        # Obtenir le filtre qui respecte les critères spécifiés
        N, _ = signal.ellipord(wPass, wStop, gPass, gStop, fs=Fe)
        print("Ordre du filtre elliptique : " + str(N))

        # Conception du filtre Butterworth
        num, denum = signal.ellip(N, rp=gPass, rs=gStop, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

        # Chercher la réponse en fréquence du filtre
        _, H_w = signal.freqz(num, denum)
        frequencies = w * Fe / (2 * np.pi)  # Conversion de rad/s en Hz

        # Afficher, en fréquence, le module et la phase du filtre
        fig, ax1 = plt.subplots()
        ax1.plot(frequencies, 20 * np.log10(np.abs(H_w)), 'red')
        ax1.set_xlabel('Freq [Hz]')
        ax1.set_ylabel('Amplitude [dB]', color='red')
        plt.show()

        # Afficher le zplane
        zplane.zplane(num, denum)

#aberration()
debruitage_comp_filtre()