import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from PIL import Image

def track_MIL(roi,array):
    history=[]
    history.append(list(roi))
    before = array[0, :, :]
    before = from_f32_to_uint8(before, roi)
    # Tracker CSRT was employed
    tracker = cv2.TrackerCSRT_create()
    tracker2 = cv2.TrackerMIL_create()
    ok = tracker.init(before, roi)
    ok2=tracker2.init(before, roi)

    for i in range(len(array)-1):
        after = array[i + 1, :, :]
        if i == 0:
            bbox = roi
        before_bbox = bbox
        after = from_f32_to_uint8(after, bbox)
        ok, bbox = tracker.update(after)
        if ok:
            history.append(list(bbox))
        elif i == 0:
            print("Tracking failed")
            history.append(before_bbox)
            bbox = before_bbox
        else:
            print("Tracking failed")
            moving_vector = history[-1] - history[-2]
            history.append(history[-1] + moving_vector)
            bbox = history[-1] + moving_vector



    # print(history)
    return history


def from_f32_to_uint8(img, groi):
    cropped_img = img[int(groi[1]):int(groi[1])+30, int(groi[0]):int(groi[0])+30]
    # plt.figure()
    # plt.imshow(cropped_img,"gray")
    # plt.show()
    # plt.close()
    img = img - 65535
    ret, thresh1 = cv2.threshold(img, np.min([np.max(cropped_img)+1000, 65535]), np.max(cropped_img), cv2.THRESH_TRUNC)
    ret, thresh2 = cv2.threshold(thresh1, np.max(np.min(cropped_img)-1000, 0), np.min(cropped_img), cv2.THRESH_TOZERO)
    # plt.imshow(thresh2, "gray")
    # plt.show()
    thresh2 = thresh2 - np.min(thresh2)
    img3 = cv2.convertScaleAbs(thresh2, alpha=(255.0/np.max(thresh2)))
    #img3 = cv2.adaptiveThreshold(img3, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2)
    # img3 = cv2.adaptiveThreshold(img3, 65535, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_TOZERO,11,2)
    # plt.figure()
    # plt.imshow(img3, "gray")
    # plt.show()
    # plt.close()


    img2 = img3.astype(np.uint8)
    #img2 = img3/np.max(img3)
    #img2 = img2.astype("float32")


    return img2