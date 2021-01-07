import configparser
import os

import cv2
import numpy as np

config = configparser.ConfigParser()
config.read("./config.ini")
save_deep = config.get("settings", "save_image_for_deep_learning")
data_base_folder = config.get("settings", "database_path")
roi_size = int(config.get("settings", "roi_size"))


def track_MIL(roi, array, sopuid, the_id, id):
    """
    CSRT tracker, the name MIL though.

    :param roi: ROI position
    :param array: 2d arrays for tracking
    :param sopuid: SOP UID of patient
    :param the_id: Original id for folder
    :param id: real id
    :return: tracking history. ROI positions.
    """
    directory_path = data_base_folder + the_id
    savedir = directory_path + "/deepimages/"
    if save_deep == "yes":
        os.makedirs(directory_path + "/deepimages", exist_ok=True)
        cv2.imwrite(savedir + str(sopuid) + str(id) + ".png",
                    array[0, int(roi[1]):int(roi[1]) + int(roi[3]), int(roi[0]):int(roi[0]) + int(roi[2])].astype(
                        np.uint16))
    history = []
    history.append(list(roi))
    before = array[0, :, :]
    before = from_uint16_to_uint8(before, roi)
    # Tracker CSRT was employed
    tracker = cv2.TrackerCSRT_create()
    tracker2 = cv2.TrackerMIL_create()
    ok = tracker.init(before, roi)
    ok2 = tracker2.init(before, roi)

    for i in range(len(array) - 1):
        after = array[i + 1, :, :]
        if i == 0:
            bbox = roi
        before_bbox = bbox
        after = from_uint16_to_uint8(after, bbox)
        ok, bbox = tracker.update(after)
        if ok:
            history.append(list(bbox))
            if save_deep == "yes":
                os.makedirs(directory_path + "/deepimages", exist_ok=True)
                cv2.imwrite(savedir + str(sopuid) + str(id) + "_" + str(i) + ".png",
                            array[i + 1, int(bbox[1]):int(bbox[1]) + int(bbox[3]),
                            int(bbox[0]):int(bbox[0]) + int(bbox[2])].astype(np.uint16))
        elif i == 0:
            print("Tracking failed")
            history.append(before_bbox)
            bbox = before_bbox
        else:
            print("Tracking failed")
            moving_vector = history[i - 1] - history[i - 2]
            history.append(history[i - 1] + moving_vector)
            bbox = history[i - 1] + moving_vector

    # print(history)
    return history



def from_uint16_to_uint8(img, groi):
    """
    Convert array from uint16 to uint8. groi holds ROI coordinate information.
    When convert the array, following procedure is also performed.
    1. Get max, min of pixel value in ROI.
    2. Median blur 3x3
    3. threshold trunc
    4. threshold tozero
    5. cv2.convertscaleAbs
    6. fastNlMeansDenoising
    :param img: image array
    :param groi: roi
    :return:
    """
    cropped_img = img[int(groi[1]) - 10:int(groi[1]) + 40, int(groi[0]) - 10:int(groi[0]) + 40]
    img = cv2.medianBlur(img, 3)
    img = img - 65535
    ret, thresh1 = cv2.threshold(img, np.min([np.max(cropped_img) + 1000, 65535]), np.max(cropped_img),
                                 cv2.THRESH_TRUNC)
    ret, thresh2 = cv2.threshold(thresh1, np.max(np.min(cropped_img) - 1000, 0), np.min(cropped_img), cv2.THRESH_TOZERO)
    thresh2 = thresh2 - np.min(thresh2)
    img3 = cv2.convertScaleAbs(thresh2, alpha=(255.0 / np.max(thresh2)))
    img2 = img3.astype(np.uint8)
    img2 = cv2.fastNlMeansDenoising(img2, None, 5, 7, 7)

    return img2
