import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def track_MIL(roi,array):
    history=[]
    # roi = tuple(roi)
    history.append(list(roi))
    before = array[0, :, :]
    before = from_f32_to_uint8(before)
    # before = self.gamma_correction(before,3.2)

    tracker = cv2.TrackerMIL_create()
    ok = tracker.init(before, roi)

    for i in range(len(array)-1):
        after = array[i + 1, :, :]
        after = from_f32_to_uint8(after)
        ok, bbox = tracker.update(after)
        if ok:
            history.append(list(bbox))

        """
        if ok:
            print(i, bbox)
            p1 = (int(bbox[0]), int(bbox[1]))
            p2 = (int(bbox[0] + bbox[2]), int(bbox[1] + bbox[3]))
            c = patches.Rectangle(xy=(bbox[0], bbox[1]), width=bbox[2], height=bbox[3], ec='white', fill=False)
            fig = plt.figure()
            ax = plt.axes()
            ax.add_patch(c)
            ax.imshow(array[i + 1, :, :], cmap='Greys', vmax=5000)
            plt.savefig("./app2_fig_3/" + str(i) + ".png")
            plt.close()
        else:
            fig = plt.figure()
            ax = plt.axes()
            ax.imshow(array[i + 1, :, :], cmap='Greys', vmax=5000)
            plt.text(2, 3, "Target Lost")
            plt.savefig("./app2_fig_3/" + str(i) + ".png")
            plt.close()
        """

    # print(history)
    return history


def from_f32_to_uint8(img):
    img2 = cv2.convertScaleAbs(img, alpha=(255.0 / 65535.0))
    img2 = img2.astype(np.uint8)
    return img2