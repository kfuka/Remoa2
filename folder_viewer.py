"""
This module is intended to list dicom file in a folder.
Listed dicom files are visualized on gui created by tkinter with SOPInstance UID.
User be asked to select two dicom files.
"""

import glob
import tkinter as tk

import pydicom

wave_colors = ['magenta', 'blue', 'red', 'green', 'cyan', 'yellow', 'black']

def get_wave_dicoms(folder_name):
    """
    get dicom with wave data
    :param folder_name: folder path
    :return: dictionary of dicom file path and acquisition time sorted with the time
    """
    dicom_list = glob.glob(folder_name + "*.dcm")
    time_and_dicom = {}
    for a_dicom in dicom_list:
        dicom_data = pydicom.dcmread(a_dicom)
        if len(dicom_data[0x5400, 0x0100][0][0x5400, 0x1010].value) > 10:
            # print(dicom_data[0x0008, 0x0018].value)
            if dicom_data[0x0008, 0x1010].value == "H-SIM1":
                direction = "H"
            else:
                direction = "V"
            time_and_dicom[a_dicom] = [dicom_data.AcquisitionTime, dicom_data[0x0008, 0x0018].value, direction]

    sorted_t_d = sorted(time_and_dicom.items(), key=lambda x: x[1], reverse=True)
    return sorted_t_d


class DicomSelectGui(tk.Frame):
    """
    Open dicom selecting GUI
    """

    def __init__(self, dicoms, master):
        """
        initialize gui
        :param dicoms:list of dicom files
        :param master: tk frame
        """
        super().__init__(master)
        self.dicom = dicoms
        # self.master = master
        self.master.title("Select dicom")
        self.master.geometry("800x400")
        self.grid()
        self.checkval = []
        self.init_draw()

    def init_draw(self):
        """
        initial drawings for gui
        :return: none
        """
        self.frame = tk.Frame(self.master)
        self.frame.grid(row=0, column=0)
        self.file_label = tk.Label(self.frame, text="File name")
        self.time_label = tk.Label(self.frame, text="Aqc. Time")
        self.direction_label = tk.Label(self.frame, text="Direction")
        self.uid_label = tk.Label(self.frame, text="SOP Instance UID")
        self.title_label = tk.Label(self.frame, text="Select ONE PAIR of Dicom")
        self.title_label.grid(row=0, column=0, columnspan=3)
        self.file_label.grid(row=1, column=3, pady=10)
        self.uid_label.grid(row=1, column=4)
        self.time_label.grid(row=1, column=2)
        self.direction_label.grid(row=1, column=1)

        for i in range(len(self.dicom)):
            bl = tk.BooleanVar()
            bl.set(False)
            b = tk.Checkbutton(self.frame, variable=bl)
            b.grid(row=i + 2, column=0)
            self.checkval.append(bl)

            label = tk.Label(self.frame, text=str(self.dicom[i][0]).split("/")[-1], fg=wave_colors[int(i/2.0)])
            label.grid(row=i + 2, column=3, padx=10)
            label2 = tk.Label(self.frame, text=str(
                self.dicom[i][1][0][:2] + ":" + self.dicom[i][1][0][2:4] + ":" + self.dicom[i][1][0][4:6]), fg=wave_colors[int(i/2.0)])
            label2.grid(row=i + 2, column=2, padx=10)
            label3 = tk.Label(self.frame, text=str(self.dicom[i][1][1]), fg=wave_colors[int(i/2.0)])
            label3.grid(row=i + 2, column=4, padx=10)
            label4 = tk.Label(self.frame, text=str(self.dicom[i][1][2]), fg=wave_colors[int(i/2.0)])
            label4.grid(row=i + 2, column=1, padx=10)
        exe_button = tk.Button(self.frame, text="open", command=self.dicom_open)
        exe_button.grid(row=len(self.dicom) + 2, column=0)

    def dicom_open(self):
        """
        open selected dicom pair.
        :return:
        """
        dicom_selection = []
        for i in range(len(self.checkval)):
            if self.checkval[i].get():
                dicom_selection.append(self.dicom[i][0])
        if len(dicom_selection) == 2:
            self.return_dicom = dicom_selection
            self.master.quit()
        else:
            label = tk.Label(self.frame, text="Choose ONE PAIR means select TWO!", fg="red")
            label.grid(row=len(self.dicom) + 3, column=0, columnspan=3)

    def get_chosen_dicoms(self):
        return self.return_dicom


if __name__ == "__main__":
    folder = "../"
    wave_dicoms = get_wave_dicoms(folder)
    dicom_window = tk.Tk()
    dicom_select_app = DicomSelectGui(wave_dicoms, master=dicom_window)
    dicom_select_app.mainloop()
    # print(dicom_select_app.return_dicom)
