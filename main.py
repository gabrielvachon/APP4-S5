import test.future_test1

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from scipy import signal

import zplane

e = np.e
pi = np.pi

def num1():
    # NUM1
    # a) afficher poles et zeros
    # Declaration des variables
    K = 1
    z1 = 0.8j
    z2 = -0.8j
    p1 = 0.95 * e ** (1j * pi / 8)
    p2 = 0.95 * e ** (-1j * pi / 8)

    # Ecrire numerateur/denominateur sous forme z^2 + z^1 + z^0
    # on peut aussi le faire avec la fonction poly
    num = np.array([1, -z1 - z2, -z1 * -z2])
    denum = np.array([1, -p1 - p2, -p1 * -p2])

    # Utiliser fonction zplane
    zplane.zplane(num, denum)

    # b) oui le filtre est stable car les poles sont dans le cercle unitaire

    # c) trouver/tracer reponse impulsionnelle en frequence H(w) du filtre
    _, H_w = signal.freqz(num, denum)

    plt.figure()
    plt.plot(np.log10(np.abs(H_w)), color="green")
    #plt.plot(w, np.angle(H_w), color="red")
    plt.savefig("rep_imp_H_w_num1c.png")
    plt.show()

    # d) generer une impulsion
    dirac = np.zeros(1000)
    dirac[int(len(dirac)/2 - 1)] = 1

    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(dirac)

    # filtre non inversé
    h = signal.lfilter(num, denum, dirac)

    plt.subplot(3, 1, 2)
    plt.plot(h)

    # e) filtre en cascade

    # filtre inverse
    h = signal.lfilter(denum, num, h)

    plt.subplot(3, 1, 3)
    plt.plot(h)
    plt.show()

def num2():
    # NUM2
    # Band-reject sur la fréquence pi/16
    w_freq = pi/16

    # Créer des zéros très proche du cercle unitaire (pas a 1 pcq le filtre va être trop pointu) et à w = pi/16
    # Créer des poles très proche des zéros mais avec une magnitude un tout petit peu moins forte pour adoucir la coupure,
    # ils doivent être aussi à w = pi/16
    num = np.poly([0.999 * e ** (1j*w_freq), 0.999 * e ** (-1j*w_freq)])
    denum = np.poly([0.95 * e ** (1j*w_freq), 0.95 * e ** (-1j*w_freq)])

    # Plot avec la fonction zplane
    zplane.zplane(num, denum)

    # Obtenir le H(z) du filtre
    _, H_w = signal.freqz(num, denum)

    # Filtre H(z) montré graphiquement
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(np.log10(np.abs(H_w)), 'red')
    ax2.plot(np.angle(H_w), 'blue')

    ax1.set_xlabel('Freq [rad/ech.]')
    ax1.set_ylabel('Amplitude [dB]', color='red')
    ax2.set_ylabel('Phase [rad]', color='blue')

    plt.show()

    # Créer la fonction x(n) à filtrer contenant deux sinus
    n = np.arange(0, 500)
    x_n = np.sin(n * pi / 16) + np.sin(n * pi / 32)

    # On passe le signal dans le filtre pour avoir la réponse du filtre
    x_n_out = signal.lfilter(num, denum, x_n)

    # Afficher la réponse du filtre avec x(n)
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(n, x_n)
    plt.subplot(2,1,2)
    plt.plot(n, x_n_out)
    plt.show()

def num3a():
    # NUM3 a)
    # Définir spécifications
    Fe = 48000      # Fréquence d'échantillonnage
    wPass = 2500    # Fréquence passante 0 - 2.5kHz
    wStop = 3500    # Fréquence coupure à 3.5kHz
    gPass = 0.2     # Variation de gain max de la bande passante (dB)
    gStop = 40      # Maximum d'atténuation de la bande coupée (dB)

    # Obtenir le filtre qui respecte les critères spécifiés
    N, _ = signal.buttord(wPass, wStop, gPass, gStop, fs=Fe)
    print("Ordre du filtre Butterworth de type passe-bas : " + str(N))

    # Conception du filtre Butterworth
    num, denum = signal.butter(N, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

    # Chercher la réponse en fréquence du filtre
    _, H_w = signal.freqz(num, denum)

    # Afficher, en fréquence, le module et la phase du filtre
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(np.log10(np.abs(H_w)), 'red')
    ax2.plot(np.angle(H_w), 'blue')
    ax1.set_xlabel('Freq [rad/ech.]')
    ax1.set_ylabel('Amplitude [dB]', color='red')
    ax2.set_ylabel('Phase [rad]', color='blue')
    plt.show()

    # Afficher le zplane
    zplane.zplane(num, denum)

def num3b():
    # NUM3 b)
    # Définir spécifications
    Fe = 48000      # Fréquence d'échantillonnage
    wPass = 2500    # Fréquence passante 0 - 2.5kHz
    wStop = 3500    # Fréquence coupure à 3.5kHz
    gPass = 0.2     # Variation de gain max de la bande passante (dB)
    gStop = 40      # Maximum d'atténuation de la bande coupée (dB)

    # Obtenir le filtre qui respecte les critères spécifiés
    N, _ = signal.cheb1ord(wPass, wStop, gPass, gStop, fs=Fe)
    print("Ordre du filtre Chebyshev de type 1 : " + str(N))

    # Conception du filtre Chebyshev type 1
    num, denum = signal.cheby1(N, rp=gPass, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

    # Chercher la réponse en fréquence du filtre
    _, H_w = signal.freqz(num, denum)

    # Afficher, en fréquence, le module et la phase du filtre
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(np.log10(np.abs(H_w)), 'red')
    ax2.plot(np.angle(H_w), 'blue')
    ax1.set_xlabel('Freq [rad/ech.]')
    ax1.set_ylabel('Amplitude [dB]', color='red')
    ax2.set_ylabel('Phase [rad]', color='blue')
    plt.show()

    # Afficher le zplane
    zplane.zplane(num, denum)

def num3c():
    # NUM3 c)
    # Définir spécifications
    Fe = 48000      # Fréquence d'échantillonnage
    wPass = 2500    # Fréquence passante 0 - 2.5kHz
    wStop = 3500    # Fréquence coupure à 3.5kHz
    gPass = 0.2     # Variation de gain max de la bande passante (dB)
    gStop = 40      # Maximum d'atténuation de la bande coupée (dB)

    # Obtenir le filtre qui respecte les critères spécifiés
    N, _ = signal.cheb2ord(wPass, wStop, gPass, gStop, fs=Fe)
    print("Ordre du filtre Butterworth de type passe-bas : " + str(N))

    # Conception du filtre Chebyshev type 2
    num, denum = signal.cheby2(N, Wn=wStop, btype='lowpass', output='ba', fs=Fe)

    # Chercher la réponse en fréquence du filtre
    _, H_w = signal.freqz(num, denum)

    # Afficher, en fréquence, le module et la phase du filtre
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(np.log10(np.abs(H_w)), 'red')
    ax2.plot(np.angle(H_w), 'blue')
    ax1.set_xlabel('Freq [rad/ech.]')
    ax1.set_ylabel('Amplitude [dB]', color='red')
    ax2.set_ylabel('Phase [rad]', color='blue')
    plt.show()

    # Afficher le zplane
    zplane.zplane(num, denum)

def formatif_q1():
    Fs = 8000
    wPass = [1000/(Fs/2), 2000/(Fs/2)]
    wStop = [750/(Fs/2), 2250/(Fs/2)]
    gPass = 0.5
    gStop = 40

    N, Wn = signal.cheb2ord(wPass, wStop, gPass, gStop)
    print("Ordre du filtre : " + str(N))

    b, a = signal.cheby2(N, gStop, Wn, 'bandpass')

    w, H_w = signal.freqz(b, a)

    plt.figure()
    plt.plot(w, 20 * np.log10(np.abs(H_w)))
    plt.xlabel("Fréquence normalisée (rad/ech)")
    plt.ylabel("Amplitude (dB)")
    plt.show()

def formatif_q2():
    #   y[n] = x[n] + 0.7 x[n-4] - 0.5 y[n-4]
    # a)
    num = [1, 0, 0, 0, 0.7]
    denum = [1, 0, 0, 0, 0.5]

    zeros = np.roots(num)
    poles = np.roots(denum)

    print("Zeros : " + str(zeros))
    print("Poles : " + str(poles))

    zplane.zplane(num, denum, 'z_plane_q2.png')

    # b)
    w, H_w = signal.freqz(num, denum)

    plt.figure()
    plt.plot(w, 20 * np.log10(abs(H_w)))
    plt.xlabel("Fréquence normalisée (rad/ech)")
    plt.ylabel("Amplitude (dB)")
    plt.show()

    # c)
    dirac = np.zeros(100)
    dirac[1] = 1
    h = signal.lfilter(num, denum, dirac)

    plt.figure()
    plt.stem(h)
    plt.xlabel("Fréquence normalisée (rad/ech)")
    plt.ylabel("Amplitude (dB)")
    plt.show()

#num3a()
formatif_q1()
formatif_q2()


