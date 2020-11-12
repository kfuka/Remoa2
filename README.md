# Remoa2
REspiratory MOtion Analyzer(Remoa)

This program was coded for the automated marker motion analysis on x-ray fluoroscopic image using the CSRT tracker in open-cv in python. 

By checking the marker position manually in the first image, then the CSRT tracker will track the markers for the following images. 

The input file should be one pair of DICOM files that have multiple continuous fluoroscopic images. Of two Dicom files, one is from horizontal, and another is from vertical exposure.


Reference
CSRT: Lukezic, A., Vojir, T., ˇCehovin Zajc, L., Matas, J., & Kristan, M. (2017). Discriminative correlation filter with channel and spatial reliability. In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (pp. 6309-6318).

Open-CV: https://opencv.org/
