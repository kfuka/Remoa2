import numpy as np


def wave_analysis(dicom_data):
    wave = dicom_data[0x5400, 0x0100][0][0x5400, 0x1010].value
    wave = [wave[i] for i in range(len(wave))]
    wave = np.array(wave[::2] - np.min(wave[::2])) / np.max(wave[::2] - np.min(wave[::2])) * 100
    freq = dicom_data[0x5400, 0x0100][0][0x003a, 0x001a].value
    wave_samples = dicom_data[0x5400, 0x0100][0][0x003a, 0x0010].value
    wave_time = [i / freq for i in range(wave_samples)]

    wave = np.delete(wave,len(wave)-1)
    wave_time = np.delete(wave_time,0)
    return wave, wave_time


if __name__ == "__main__":
    import pydicom
    import matplotlib.pyplot as plt

    dicom_file = "../710509.dcm"
    dicom = pydicom.dcmread(dicom_file)
    wave, wave_time = wave_analysis(dicom)
    print(wave)
    print(wave_time)
    plt.plot(wave_time, wave, "o")
    plt.show()
