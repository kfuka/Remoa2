import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def track_MIL(roi,array):
    history=[]
    global groi
    groi = roi
    # roi = tuple(roi)
    history.append(list(roi))
    before = array[0, :, :]
    before = from_f32_to_uint8(before)
    # before = self.gamma_correction(before,3.2)

    # tracker = cv2.TrackerMIL_create()
    # tracker = cv2.Tracker_create("KCF")
    # Tracker CSRT was employed
    tracker = cv2.TrackerCSRT_create()
    ok = tracker.init(before, roi)

    for i in range(len(array)-1):
        after = array[i + 1, :, :]
        after = from_f32_to_uint8(after)
        ok, bbox = tracker.update(after)
        if ok:
            history.append(list(bbox))
        else:
            history.append([0, 0, 0, 0])
            print("Tracking failed")
            print(ok)
            print(bbox)


    # print(history)
    return history


def from_f32_to_uint8(img):

    img3 = cv2.convertScaleAbs(img, alpha=(255.0/np.max(img)))
    shrinked_img = cv2.adaptiveThreshold(img3, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11,2)

    img2 = shrinked_img.astype(np.uint8)


    return img2