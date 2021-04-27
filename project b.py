# %%

import matplotlib.pyplot as plt  ### plotting things
import numpy as np  ## one of python's main maths packages
import pandas as pd  ## for reading in our data
from scipy.optimize import curve_fit  ## for fitting a line to our data
import matplotlib.ticker as ticker  ## this one lets us change some parameters in our plots.
%matplotlib inline
%config InlineBackend.figure_format = 'retina'
from IPython.display import display, Markdown
from scipy import fftpack

repeats = ["a", "b", "c", "d", "e", "f"]
pres = int(input('What is the pressure (hPa)'))
pressure = pres * 100
volumes = ["10", "20", "30", "40", "50", "60", "70", "80", "90", "100"]

for VOL in volumes:
    for REP in repeats:
        g = "N2 " + VOL + "ml-" + REP + ".csv"
        if g == "N2 90ml-f.csv":
            pass
        else:
            data = pd.read_csv(g,
                               header=1,
                               names=('time', 'emf'))

        # print(data.time)

        time = []
        emf = []
        count = -1
        for each in data.time:  # Loading in data into lists time and emf
            count = count + 1
            if str(each) != 'second':
                if each > 0:
                    time.append(each)
                    emf.append(data.emf[count])

        timestep = round((time[1] - time[0]), 5)
        # print(timestep)

        Signal = fftpack.fft(emf)  # Fourier Transform the signal to inverse domain, frequency

        Amplitude = np.abs(Signal)  # Find the absolute value of the signal
        sample_freq = fftpack.fftfreq(Signal.size, d=timestep)  # Return Discrete fourier transform sampling frequencies

        Amp_freq = np.array([Amplitude, sample_freq])  # Gives frequencies of all amplitudes
        Amp_pos = Amp_freq[0, :].argmax()  # Finds Max amplitude
        Fd = Amp_freq[1, Amp_pos]  # Finds the corresponding frequency for the max amplitude, dampening frequency

        # Logarithic Decrement to get the Dampening Ratio, Zeta

        count = -1
        repeat = 0
        Vpeaks = []
        tpeaks = []
        peak = 0
        for each in time:
            count = count + 1
            if emf[count] > peak:
                peak = float(emf[count])  # Code to find the peaks of each oscillation and store the co-ordinates
                tpeak = float(each)
            if peak != 0:  # Voltage Values are stored in Vpeaks and correspondin time values in tpeaks
                if emf[count] < 0:
                    Vpeaks.append(peak)
                    tpeaks.append(tpeak)
                    peak = 0

        # print(Vpeaks)

        Vp = []
        tp = []

        x = []  # stores integer values 1 - range
        lg = []  # stores logged peak values

        # for n in range(len(Vpeaks)-1):
        for n in range(5):
            lg.append(np.log(
                Vpeaks[0] / Vpeaks[n + 1]))  # Code finds the log of peak n over peak n+1, as required in the formula
            x.append(n + 1)  # Note, only first 10 log values taken, change range value to adjust

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, lg, marker='o', linestyle='none', color='black')
        ax.set_xlabel('Peak Number n')
        ax.set_ylabel('Log(An/An+1)')


        # Plots the logs against peak number n to get a gradient
        # gradient should equal (2 x PI x Zeta) / SQRT(1-Zeta^2)
        # gradient can be rearranged to get Zeta, which is the dampening RATIO, equal to Gamma over the natural frequency
        # Gamma is gamma from damped harmonic motion formula (damping coefficient over 2 x mass) (gamma = c/2m)
        # Zeta should be less than one as the system is underdamped

        def line(x, slope, intercept):
            return slope * x + intercept


        popt, pcov = curve_fit(line, x, lg)
        slope = popt[0]
        intercept = popt[1]
        err_slope = np.sqrt(float(pcov[0][0]))
        err_intercept = np.sqrt(float(pcov[1][1]))  # Code to plot Line Of Best Fit and give the lines equation
        # print('gradient = ', slope,'+/-',err_slope)
        # print('intercept = ', intercept,' +/- ',err_intercept)

        y = []

        for each in x:
            point = each * slope + intercept
            y.append(point)
        ax.plot(x, y, linestyle='--', color='orange')
        plt.show()

        zeta = slope / (np.sqrt(((2 * np.pi) ** 2) + (slope ** 2)))

        fig.savefig('LogGraph.png', dpi=200)

        # THIS IS WORKING FULLY
        # Note not enough values for an accurate dampening constant of 20ml nitrogen.

        # Can overide all constants here

        Wd = 2 * np.pi * Fd
        zeta = zeta
        Wn = Wd / (np.sqrt(1 - (zeta ** 2)))
        Fn = Wn / (np.pi * 2)
        # print('Dampening Ratio = ', zeta)
        # print('Dampening Angular Frequency = ', Wd)
        # print('Natural Angular Frequency = ', Wn)
        Amp = 0.4
        phi = 1.364 + np.pi

        x = []
        for each in time:
            dis = Amp * np.exp(-zeta * Wn * each) * np.cos(Wd * each + phi)
            x.append(dis)

        '''fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
        ax.plot(time,x)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Displacement (m)')'''

        # Quick plot of modelled displacement against time

        v = []
        for each in time:
            vel = -Amp * np.exp(-zeta * Wn * each) * (
            (Wd * np.sin(Wd * each + phi) - Wn * zeta * np.cos(Wd * each + phi)))
            v.append(vel)
        '''
        fig2 = plt.figure(figsize=(8,8))
        ax = fig2.add_subplot(1,1,1)
        ax.plot(time,v)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Velocity (ms^-1)')'''

        # FINDING GAMMA

        volume = (float(VOL) + 6.3) * 0.000001
        print(volume)
        area = 0.00091649
        mass = 0.10668
        # print("\n\n\n\n")
        gamma = ((Wn ** 2) * mass * volume) / (pressure * (area ** 2))
        # print("Gamma = ", gamma)

        a = 0.00091649
        da = 0.00000054
        m = 0.10668
        dm = 0.01 / 106.68
        v = volume
        dv = (0.05 * 0.000001) / v
        p = pressure
        dp = 100 / pressure

        dx = dp + (2 * da) + dm + dv

        print(dx)

        filedata = VOL + "," + str(gamma) + "," + str(Wn) + "," + str(dx)

        file_object = open('N2 values.txt', 'a')
        file_object.write(filedata)
        file_object.write("\n")
        file_object.close()
        gamma = 0

# %%


# %%


