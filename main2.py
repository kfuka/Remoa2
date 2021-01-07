import configparser
import os
import subprocess
import time
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import tkinter.font as font
import tkinter.ttk as ttk

import matplotlib.animation as animation
import matplotlib.patches as pat
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import statsmodels.api as sm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from scipy import signal
from statsmodels.sandbox.regression.predstd import wls_prediction_std

import csrt_tracker
import folder_viewer
import unite_class
import wave_analysis

config = configparser.ConfigParser()
config.read('./config.ini')
roi_size = int(config.get("settings", "roi_size"))
save_mpeg4 = config.get("settings", "save_mpeg4")
wave_colors = ["g", "b", "c", "m", "y", "k", "w"]
wave_colors_full = ["green", "blue", "cyan", "magenta", "yellow", "black", "white"]
data_base_folder = config.get("settings", "database_path")
spacing_correction_factor = 1.5 * 1550 / 2111.63
ct_val_min = -2000
ct_val_max = 10000
ct_val_scale = 0.5
preset_list = ["Bone", "Lung", "Soft"]



class Application(tk.Frame):
    """
    Main application
    """

    def __init__(self, master=None):
        """
        initialize
        :param master:
        """
        self.roi_center1 = []
        self.roi_center2 = []
        self.for_wave1 = []
        self.for_wave1_raw = []
        self.for_wave2 = []
        self.for_wave2_raw = []
        super().__init__(master)
        self.master = master
        self.master.bind("<Key>", self.on_key_press)
        self.master.bind("<KeyRelease>", self.on_key_release)
        self.master.title('Respiratory motion analyzer')
        self.grid()
        self.init_draw()
        self.create_widgets()
        self.start_up()
        self.ctrlpress1 = None
        self.ctrlpress2 = None
        self.ctrl_is_held = False
        self.shift_is_held = False
        self.rpress1 = None
        self.rpress2 = None
        self.calculated = False
        self.closest_point = []
        self.closest_point2 = []
        self.current_slice = 0
        self.imshow1 = None
        self.imshow2 = None
        self.timer = ""
        self.max_x1, self.max_x2, self.max_y1, self.max_y2 = [], [], [], []
        self.current_patient_dicom_folder = []
        self.id = []
        self.is_in_loop = False


    def create_widgets(self):
        """
        create tk wigdets
        :return: nothing
        """
        self.canvas_frame1 = tk.Frame(self.master)
        self.canvas_frame1.grid(row=0, column=0, rowspan=2)
        self.canvas_frame2 = tk.Frame(self.master)
        self.canvas_frame2.grid(row=0, column=2, rowspan=2)
        self.control_frame = tk.Frame(self.master)
        self.control_frame.grid(row=0, column=1, rowspan=2)
        self.control_frame2 = tk.Frame(self.master)
        self.control_frame2.grid(row=0, column=3, rowspan=2)
        self.wave_frame = tk.Frame(self.master)
        self.scroll_frame = tk.Frame(self.master)
        self.scroll_frame.grid(row=2, column=0, columnspan=4)
        self.wave_frame.grid(row=3, column=0, columnspan=1)
        self.wave_frame2 = tk.Frame(self.master)
        self.wave_frame2.grid(row=3, column=2, columnspan=1)
        self.main_frame = tk.Frame(self.master)
        self.main_frame.grid(row=0, column=4, rowspan=4, sticky=tk.N)

        self.canvas1 = FigureCanvasTkAgg(self.fig1, self.canvas_frame1)
        self.canvas1.draw()
        # self.canvas1.get_tk_widget().grid(row=0, column=0)
        self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas1.mpl_connect('button_press_event', self.onclick1)
        self.canvas1.mpl_connect('scroll_event', self.mouse_scrolled1)
        self.canvas1.mpl_connect('motion_notify_event', self.on_motion1)
        self.canvas1.mpl_connect('button_release_event', self.on_release1)
        self.canvas1.mpl_connect('key_press_event', self.on_key1)

        self.canvas2 = FigureCanvasTkAgg(self.fig2, self.canvas_frame2)
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas2.mpl_connect('button_press_event', self.onclick2)
        self.canvas2.mpl_connect('scroll_event', self.mouse_scrolled2)
        self.canvas2.mpl_connect('motion_notify_event', self.on_motion2)
        self.canvas2.mpl_connect('button_release_event', self.on_release2)

        self.wave_canvas = FigureCanvasTkAgg(self.wave_fig, self.wave_frame)
        self.wave_canvas.draw()
        self.wave_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.wave_canvas2 = FigureCanvasTkAgg(self.wave_fig2, self.wave_frame2)
        self.wave_canvas2.draw()
        self.wave_canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.x_v = tk.IntVar()
        self.x_v2 = tk.IntVar()
        self.slice_num = tk.IntVar()
        self.x_scale = tk.Scale(self.control_frame,
                                variable=self.x_v,
                                from_=ct_val_max,
                                to=ct_val_min,
                                resolution=100,
                                orient=tk.VERTICAL,
                                length=400,
                                showvalu=0,
                                command=self.draw_plot)
        self.x_scale.grid(row=0, column=0, padx=0, pady=0, sticky=tk.W + tk.N + tk.S)

        # self.ct1max = tk.Entry(self.control_frame, text="100", width=10, bg="red")
        # self.ct1max.grid(row=0, column=0, columnspan=2)

        self.x_scale2 = tk.Scale(self.control_frame2,
                                 variable=self.x_v2,
                                 from_=ct_val_max,
                                 to=ct_val_min,
                                 resolution=100,
                                 orient=tk.VERTICAL,
                                 length=400,
                                 showvalu=0,
                                 command=self.draw_plot)
        self.x_scale2.grid(row=0, column=0, padx=0, pady=0, sticky=tk.W + tk.N + tk.S)

        self.slice_scroll = tk.Scale(self.scroll_frame,
                                     variable=self.slice_num,
                                     from_=0,
                                     to=29,
                                     resolution=1,
                                     orient=tk.HORIZONTAL,
                                     length=500,
                                     showvalu=0,
                                     command=self.draw_plot)
        self.slice_scroll.grid(row=0, column=2, padx=1, pady=1, sticky=tk.N + tk.S)
        self.left_edge_button = tk.Button(self.scroll_frame, text="<<", height=2, width=4, command=self.go_left_edge)
        self.left_edge_button.grid(row=0, column=0, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.left_button = tk.Button(self.scroll_frame, text="<", height=2, width=4, command=self.go_left)
        self.left_button.grid(row=0, column=1, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.right_button = tk.Button(self.scroll_frame, text=">", height=2, width=4, command=self.go_right)
        self.right_button.grid(row=0, column=3, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.right_edge_button = tk.Button(self.scroll_frame, text=">>", height=2, width=4, command=self.go_right_edge)
        self.right_edge_button.grid(row=0, column=4, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.loop_on_button = tk.Button(self.scroll_frame, text="∞", height=2, width=4, command=self.loop_on)
        self.loop_on_button.grid(row=0, column=5, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.loop_off_button = tk.Button(self.scroll_frame, text="II", height=2, width=4, command=self.loop_off)
        self.loop_off_button.grid(row=0, column=6, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)
        self.preset_txt = tk.StringVar()
        self.preset_txt.set("preset")
        # self.preset = tk.OptionMenu(self.scroll_frame, preset_txt, "Bone", "Lung", "Soft", command=self.preset)
        # preslist = self.preset_list
        self.preset = tk.OptionMenu(self.scroll_frame, self.preset_txt, *preset_list, command=self.preset_loop)
        self.preset.grid(row=0, column=7, padx=1, pady=1, sticky=tk.W + tk.E + tk.N + tk.S)


        self.y_v = tk.IntVar()
        self.y_v2 = tk.IntVar()
        self.y_scale = tk.Scale(self.control_frame,
                                variable=self.y_v,
                                from_=ct_val_max,
                                to=ct_val_min,
                                resolution=100,
                                orient=tk.VERTICAL,
                                showvalue=0,
                                command=self.draw_plot)
        self.y_scale.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S)

        self.y_scale2 = tk.Scale(self.control_frame2,
                                 variable=self.y_v2,
                                 from_=ct_val_max,
                                 to=ct_val_min,
                                 resolution=100,
                                 orient=tk.VERTICAL,
                                 showvalue=0,
                                 command=self.draw_plot)
        self.y_scale2.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S)

        plus_val = 5
        self.open_button = tk.Button(self.main_frame, text="Open DICOM", height=3, command=self.open_dicom)
        self.open_button.grid(row=1+plus_val, column=0, padx=2, pady=2, rowspan=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.show_rois = tk.Text(self.main_frame, height=30, width=10, wrap=tk.CHAR)
        self.show_rois.grid(row=0, column=0, padx=0, pady=0, columnspan=6, sticky=tk.N + tk.EW)
        self.update_show_rois("1. Open DICOM\n")
        self.update_show_rois("2. Click markers\n")
        self.update_show_rois("3. Calculate!\n")
        self.update_show_rois("---------------\n")

        my_font = font.Font(self.main_frame, family="Arial", size=12, weight="bold")

        self.label1 = tk.Label(self.main_frame, text="#:", font=my_font)
        self.label2 = tk.Label(self.main_frame, text="L-R", font=my_font)
        self.label3 = tk.Label(self.main_frame, text="A-P", font=my_font)
        self.label4 = tk.Label(self.main_frame, text="S-I", font=my_font)
        self.label5 = tk.Label(self.main_frame, text="3D (mm)", font=my_font)
        self.label6 = tk.Label(self.main_frame, text="", font=my_font)

        self.label1.grid(row=1, column=0, padx=2, pady=2)
        self.label2.grid(row=1, column=1, padx=2, pady=2)
        self.label3.grid(row=1, column=2, padx=2, pady=2)
        self.label4.grid(row=1, column=3, padx=2, pady=2)
        self.label5.grid(row=1, column=4, padx=2, pady=2)
        self.label6.grid(row=1, column=5, padx=2, pady=2)

        self.sv11 = tk.StringVar()
        self.sv12 = tk.StringVar()
        self.sv13 = tk.StringVar()
        self.sv14 = tk.StringVar()
        self.sv15 = tk.StringVar()

        self.sv21 = tk.StringVar()
        self.sv22 = tk.StringVar()
        self.sv23 = tk.StringVar()
        self.sv24 = tk.StringVar()
        self.sv25 = tk.StringVar()

        self.sv31 = tk.StringVar()
        self.sv32 = tk.StringVar()
        self.sv33 = tk.StringVar()
        self.sv34 = tk.StringVar()
        self.sv35 = tk.StringVar()

        self.sv41 = tk.StringVar()
        self.sv42 = tk.StringVar()
        self.sv43 = tk.StringVar()
        self.sv44 = tk.StringVar()
        self.sv45 = tk.StringVar()

        self.label11 = tk.Label(self.main_frame, textvariable=self.sv11, font=my_font)
        self.label12 = tk.Label(self.main_frame, textvariable=self.sv12, font=my_font)
        self.label13 = tk.Label(self.main_frame, textvariable=self.sv13, font=my_font)
        self.label14 = tk.Label(self.main_frame, textvariable=self.sv14, font=my_font)
        self.label15 = tk.Label(self.main_frame, textvariable=self.sv15, font=my_font)

        self.label21 = tk.Label(self.main_frame, textvariable=self.sv21, font=my_font)
        self.label22 = tk.Label(self.main_frame, textvariable=self.sv22, font=my_font)
        self.label23 = tk.Label(self.main_frame, textvariable=self.sv23, font=my_font)
        self.label24 = tk.Label(self.main_frame, textvariable=self.sv24, font=my_font)
        self.label25 = tk.Label(self.main_frame, textvariable=self.sv25, font=my_font)

        self.label31 = tk.Label(self.main_frame, textvariable=self.sv31, font=my_font)
        self.label32 = tk.Label(self.main_frame, textvariable=self.sv32, font=my_font)
        self.label33 = tk.Label(self.main_frame, textvariable=self.sv33, font=my_font)
        self.label34 = tk.Label(self.main_frame, textvariable=self.sv34, font=my_font)
        self.label35 = tk.Label(self.main_frame, textvariable=self.sv35, font=my_font)

        self.label41 = tk.Label(self.main_frame, textvariable=self.sv41, font=my_font)
        self.label42 = tk.Label(self.main_frame, textvariable=self.sv42, font=my_font)
        self.label43 = tk.Label(self.main_frame, textvariable=self.sv43, font=my_font)
        self.label44 = tk.Label(self.main_frame, textvariable=self.sv44, font=my_font)
        self.label45 = tk.Label(self.main_frame, textvariable=self.sv45, font=my_font)

        self.label11.grid(row=2, column=0, padx=2, pady=2)
        self.label12.grid(row=2, column=1, padx=2, pady=2)
        self.label13.grid(row=2, column=2, padx=2, pady=2)
        self.label14.grid(row=2, column=3, padx=2, pady=2)
        self.label15.grid(row=2, column=4, padx=2, pady=2)

        self.label21.grid(row=3, column=0, padx=2, pady=2)
        self.label22.grid(row=3, column=1, padx=2, pady=2)
        self.label23.grid(row=3, column=2, padx=2, pady=2)
        self.label24.grid(row=3, column=3, padx=2, pady=2)
        self.label25.grid(row=3, column=4, padx=2, pady=2)

        self.label31.grid(row=4, column=0, padx=2, pady=2)
        self.label32.grid(row=4, column=1, padx=2, pady=2)
        self.label33.grid(row=4, column=2, padx=2, pady=2)
        self.label34.grid(row=4, column=3, padx=2, pady=2)
        self.label35.grid(row=4, column=4, padx=2, pady=2)

        self.label41.grid(row=5, column=0, padx=2, pady=2)
        self.label42.grid(row=5, column=1, padx=2, pady=2)
        self.label43.grid(row=5, column=2, padx=2, pady=2)
        self.label44.grid(row=5, column=3, padx=2, pady=2)
        self.label45.grid(row=5, column=4, padx=2, pady=2)

        self.roi_clear = tk.Button(self.main_frame, text="ROI Clear", height=1, command=self.clear_roi)
        self.roi_clear.grid(row=5+plus_val, column=3, padx=2, pady=2, columnspan=3, sticky=tk.EW)

        self.calc_exe = tk.Button(self.main_frame, text="Calculate!", height=3, command=self.calc)
        self.calc_exe.grid(row=3+plus_val, column=0, padx=2, pady=2, rowspan=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.print_report = tk.Button(self.main_frame, text="REPORT", command=self.report)
        self.print_report.grid(row=3+plus_val, column=3, padx=2, pady=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.patient_folder = tk.Button(self.main_frame, text="Patient folder", command=self.open_patient_folder)
        self.patient_folder.grid(row=4+plus_val, column=3, padx=2, pady=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.save_data_b = tk.Button(self.main_frame, text="Save Data", height=3, command=self.save_data)
        self.save_data_b.grid(row=5+plus_val, column=0, padx=2, pady=2, rowspan=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.show_animation = tk.Button(self.main_frame, text="Save Animation", height=3, command=self.animate)
        self.show_animation.grid(row=1+plus_val, column=3, padx=2, pady=2, rowspan=2, columnspan=3, sticky=tk.EW)

        self.quit_button = tk.Button(self.main_frame, text="quit", command=self.go_quit)
        self.quit_button.grid(row=6+plus_val, column=3, padx=2, pady=2, columnspan=3, sticky=tk.EW + tk.NS)

        self.save_data_b["state"] = "disable"
        self.calc_exe["state"] = "disable"
        self.show_animation["state"] = "disable"
        self.patient_folder["state"] = "disable"
        self.loop_off_button["state"] = "disable"

    def preset_loop(self, value):
        if value == "Bone":
            self.y_v.set(2200)  # vmax
            self.x_v.set(200)  # vmin
            self.y_v2.set(2200)  # vmax2
            self.x_v2.set(200)  # vmin2
        elif value == "Lung":
            self.y_v.set(3000)  # vmax
            self.x_v.set(0)  # vmin
            self.y_v2.set(5000)  # vmax2
            self.x_v2.set(1000)  # vmin2
        else:
            self.y_v.set(1000)  # vmax
            self.x_v.set(-500)  # vmin
            self.y_v2.set(2000)  # vmax2
            self.x_v2.set(-500)  # vmin2

        self.draw_plot()

    def go_left_edge(self):
        self.slice_num.set(0)
        self.draw_plot()
        return

    def go_left(self):
        self.wheel_slice(-1)
        return

    def go_right(self):
        self.wheel_slice(1)
        return

    def go_right_edge(self):
        self.slice_num.set(30)
        self.draw_plot()
        return

    def loop_on(self):
        self.show_loop()

    def show_loop(self):
        self.is_in_loop = True
        self.loop_on_button["state"] = "disable"
        self.loop_off_button["state"] = "active"
        global timers
        current_slice = self.slice_num.get()
        if current_slice == 29:
            current_slice = 0
        else:
            current_slice += 1
        self.slice_num.set(current_slice)
        self.draw_plot()
        timers = root.after(1, self.show_loop)

    def loop_off(self):
        self.loop_on_button["state"] = "active"
        self.loop_off_button["state"] = "disable"
        self.is_in_loop = False
        global timers
        root.after_cancel(timers)
        timers = None

    def go_quit(self):
        """
        quit application
        :return:
        """
        root.quit()
        root.destroy()

    def open_patient_folder(self):
        subprocess.run('explorer {}'.format(data_base_folder + self.id), cwd=os.getcwd())

    def save_data(self):
        # ファイル名はUnix time
        ut = time.time()
        directory_path = data_base_folder + self.id
        if not os.path.isdir(directory_path):
            os.mkdir(directory_path)
        np.savez(directory_path + "/" + str(self.SOPUID[0]),
                 time1=self.wave_time1,
                 x1=self.wave1,
                 y1=self.for_wave1,
                 time2=self.wave_time2,
                 x2=self.wave2,
                 y2=self.for_wave2,
                 UID=self.SOPUID,
                 max_x1=self.max_x1,
                 max_x2=self.max_x2,
                 max_y1=self.max_y1,
                 max_y2=self.max_y2,
                 actime=self.actime,
                 study_date=self.study_date)
        fig = plt.figure(figsize=(8, 4))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax1.imshow(self.array1[0, :, :], vmin=self.vmin, vmax=self.vmax, cmap="Greys")
        ax1.axis("off")
        fig.tight_layout()
        if self.roi_center1:
            roinum = 1
            for i, point in enumerate(self.roi_center1):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                ax1.add_patch(rect)
                ax1.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
                         color=wave_colors[i])
                roinum += 1
        ax1.set_xlim([self.xlimlow1, self.xlimhigh1])
        ax1.set_ylim([self.ylimhigh1, self.ylimlow1])

        ax2.imshow(self.array2[0, :, :], vmin=self.vmin2, vmax=self.vmax2, cmap="Greys")
        ax2.axis("off")
        # fig1.tight_layout()
        if self.roi_center2:
            roinum = 1
            for i, point in enumerate(self.roi_center2):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                ax2.add_patch(rect)
                ax2.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
                         color=wave_colors[i])
                roinum += 1
        ax2.set_xlim([self.xlimlow2, self.xlimhigh2])
        ax2.set_ylim([self.ylimhigh2, self.ylimlow2])

        fig.savefig(directory_path + "/" + str(self.SOPUID[0]) + "fig1.png")
        self.update_show_rois("-----------\n")
        self.update_show_rois("Data saved.\n")

    def report(self):
        if self.id:
            patient_save_folder = tk.filedialog.askdirectory(initialdir=data_base_folder + self.id)
        else:
            patient_save_folder = tk.filedialog.askdirectory(initialdir="./data_base/")
        patient = unite_class.build_report(patient_save_folder)
        self.update_show_rois("Report created.\n")

    def calc(self):
        self.show_animation["state"] = "active"
        self.save_data_b["state"] = "active"
        max_x1 = []
        max_x2 = []
        max_y1 = []
        max_y2 = []
        if not self.roi_center2:
            print("Select ROI first")
        else:
            self.show_rois.tag_config('over', foreground='red', font='TkDefaultFont 12 bold')
            self.update_show_rois("--------results-----------" + "\n")
            # self.update_show_rois("#: x, y max shift" + " in mm\n")
            self.marker_chase = []
            self.marker_chase2 = []
            # self.update_show_rois("Vert L-R, S-I\n")

            for i in range(len(self.roi_center1)):
                returned = csrt_tracker.track_MIL(
                    (self.roi_center1[i][0] - roi_size / 2, self.roi_center1[i][1] - roi_size / 2, roi_size, roi_size),
                    self.array1, self.SOPUID[0], self.id, i)
                x_shifts = np.array([it[0] for it in returned])
                y_shifts = np.array([it[1] for it in returned])
                # savgol fileter
                x_shifts = signal.savgol_filter(x_shifts, 7, 5, mode="nearest")
                y_shifts = signal.savgol_filter(y_shifts, 7, 5, mode="nearest")
                self.marker_chase.append([x_shifts, y_shifts])
                max_x_move = (np.max(x_shifts) - np.min(x_shifts)) * self.pixel_spacing1[0]
                max_y_move = (np.max(y_shifts) - np.min(y_shifts)) * self.pixel_spacing1[1]
                max_x1.append(max_x_move)
                max_y1.append(max_y_move)
                calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                calcd_wave_raw = calcd_wave
                calcd_wave = np.delete(calcd_wave, 0)  # ここできってる
                self.for_wave1.append(calcd_wave)
                self.for_wave1_raw.append(calcd_wave_raw)
            # self.update_show_rois("Horiz A-P, S-I\n")
            for i in range(len(self.roi_center2)):
                returned = csrt_tracker.track_MIL(
                    (self.roi_center2[i][0] - roi_size / 2, self.roi_center2[i][1] - roi_size / 2, roi_size, roi_size),
                    self.array2, self.SOPUID[0], self.id, i + 4)
                x_shifts = np.array([it[0] for it in returned])
                y_shifts = np.array([it[1] for it in returned])
                # savgol fileter
                x_shifts = signal.savgol_filter(x_shifts, 7, 5, mode="nearest")
                y_shifts = signal.savgol_filter(y_shifts, 7, 5, mode="nearest")
                self.marker_chase2.append([x_shifts, y_shifts])
                max_x_move = (np.max(x_shifts) - np.min(x_shifts)) * self.pixel_spacing2[0]
                max_y_move = (np.max(y_shifts) - np.min(y_shifts)) * self.pixel_spacing2[1]
                max_x2.append(max_x_move)
                max_y2.append(max_y_move)

                calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                calcd_wave_raw = calcd_wave
                calcd_wave = np.delete(calcd_wave, 0)
                self.for_wave2.append(calcd_wave)
                self.for_wave2_raw.append(calcd_wave_raw)

            LR_from_center = []
            AP_from_center = []
            for x1 in self.marker_chase:
                LR_from_center.append((756.0 - np.average(x1[0])) * self.pixel_spacing1[0])
            for x2 in self.marker_chase2:
                AP_from_center.append((np.average(x2[0]) - 756.0) * self.pixel_spacing2[0])

            c_max_x1, c_max_y1, c_max_x2, c_max_y2 = [], [], [], []
            for n, the_x1 in enumerate(max_x1):
                c_max_x1.append(self.LRAP_correction(the_x1, AP_from_center[n]))
                c_max_y1.append(self.LRAP_correction(max_y1[n], AP_from_center[n]))
                c_max_x2.append(self.LRAP_correction(max_x2[n], LR_from_center[n]))
                c_max_y2.append(self.LRAP_correction(max_y2[n], LR_from_center[n]))

            self.max_x1 = c_max_x1
            self.max_x2 = c_max_x2
            self.max_y1 = c_max_y1
            self.max_y2 = c_max_y2


        self.plot_wave()
        self.calculated = True
        self.change_SI_table()

    def get_maxshift_from_marker_chase(self):
        max_x1, max_y1, max_x2, max_y2 = [], [], [], []
        for j,i in enumerate(self.marker_chase):
            x_shifts, y_shifts = i
            x2_shifts, y2_shifts = self.marker_chase2[j]
            max_x1.append((np.max(x_shifts) - np.min(x_shifts)) * self.pixel_spacing1[0])
            max_y1.append((np.max(y_shifts) - np.min(y_shifts)) * self.pixel_spacing1[1])
            max_x2.append((np.max(x2_shifts) - np.min(x2_shifts)) * self.pixel_spacing1[0])
            max_y2.append((np.max(y2_shifts) - np.min(y2_shifts)) * self.pixel_spacing1[1])

        LR_from_center = []
        AP_from_center = []
        for x1 in self.marker_chase:
            LR_from_center.append((756.0 - np.average(x1[0])) * self.pixel_spacing1[0])
        for x2 in self.marker_chase2:
            AP_from_center.append((np.average(x2[0]) - 756.0) * self.pixel_spacing2[0])

        c_max_x1, c_max_y1, c_max_x2, c_max_y2 = [], [], [], []
        for n, the_x1 in enumerate(max_x1):
            c_max_x1.append(self.LRAP_correction(the_x1, AP_from_center[n]))
            c_max_y1.append(self.LRAP_correction(max_y1[n], AP_from_center[n]))
            c_max_x2.append(self.LRAP_correction(max_x2[n], LR_from_center[n]))
            c_max_y2.append(self.LRAP_correction(max_y2[n], LR_from_center[n]))

        return c_max_x1, c_max_x2, c_max_y1, c_max_y2

    def clear_SI_table(self):
        if self.max_x1:
            for i in range(len(self.max_x1)):
                exec('self.sv{}1.set("")'.format(i + 1))
                exec('self.sv{}2.set("")'.format(i + 1))
                exec('self.sv{}3.set("")'.format(i + 1))
                exec('self.sv{}4.set("")'.format(i + 1))
                exec('self.sv{}5.set("")'.format(i + 1))
                exec('self.label{}2.config(bg="SystemButtonFace", fg="black")'.format(i + 1))
                exec('self.label{}3.config(bg="SystemButtonFace", fg="black")'.format(i + 1))
                exec('self.label{}4.config(bg="SystemButtonFace", fg="black")'.format(i + 1))
                exec('self.label{}5.config(bg="SystemButtonFace", fg="black")'.format(i + 1))

    def change_SI_table(self):
        if self.max_x1:
            max_x1, max_x2, max_y1, max_y2 = self.get_maxshift_from_marker_chase()
            self.max_x1, self.max_x2, self.max_y1, self.max_y2 = max_x1, max_x2, max_y1, max_y2
            for i in range(len(max_x1)):
                exec('self.sv{}1.set(str(i + 1))'.format(i+1))
                exec('self.label{}1.config(fg="{}")'.format(i+1, wave_colors_full[i]))

                exec('self.sv{}2.set({:.2f})'.format(i+1, max_x1[i]))
                if max_x1[i] <= 5:
                    exec('self.label{}2.config(bg="green", fg="white")'.format(i + 1))
                elif max_x1[i] > 5:
                    exec('self.label{}2.config(bg="red", fg="white")'.format(i + 1))

                exec('self.sv{}3.set({:.2f})'.format(i+1, max_x2[i]))
                if max_x2[i] <= 5:
                    exec('self.label{}3.config(bg="green", fg="white")'.format(i + 1))
                elif max_x2[i] > 5:
                    exec('self.label{}3.config(bg="red", fg="white")'.format(i + 1))

                three = max(max_y1[i], max_y2[i])
                exec('self.sv{}4.set({:.2f})'.format(i+1, three))
                if three <= 5:
                    exec('self.label{}4.config(bg="green", fg="white")'.format(i + 1))
                elif three > 5:
                    exec('self.label{}4.config(bg="red", fg="white")'.format(i + 1))

                four = np.sqrt(max_x1[i] ** 2 + max_x2[i] ** 2 + max(max_y1[i], max_y2[i]) ** 2)
                exec('self.sv{}5.set({:.2f})'.format(i+1, four))
                if four <= 5:
                    exec('self.label{}5.config(bg="green", fg="white")'.format(i + 1))
                elif four > 5:
                    exec('self.label{}5.config(bg="red", fg="white")'.format(i + 1))


    def LRAP_correction(self, val_iso, distance):
        """
        return correction factor for ISO center to marker position
        :param val_iso: source to isocenter distance in mm
        :param distance: isocenter to maeker distance in mm
        :return: correction factor
        """
        return (1550.0 + distance) / 1550.0 * val_iso

    def close_dicom(self):
        self.init_draw()
        self.canvas1.draw()
        self.canvas2.draw()
        self.wave_ax.cla()
        self.wave_fig.clf()
        self.wave_canvas.draw()
        self.clear_roi()

    def correlation(self):
        self.cor_fig = Figure(figsize=(12, 6))
        self.cor_ax1 = self.cor_fig.add_subplot(121)
        self.cor_ax2 = self.cor_fig.add_subplot(122)
        self.show_rois.insert("end", '----correlation--------' + "\n")
        self.show_rois.insert("end", '#: fit, lwr, upr, Predict' + "\n")

        for i in range(len(self.roi_center1)):
            # print(np.corrcoef(self.wave1, self.for_wave1)[0,1])
            self.cor_ax1.plot(self.wave1, self.for_wave1[i], "o", mfc="None", c=wave_colors[i],
                              label="#" + str(i + 1))
            X = sm.add_constant(self.wave1)
            re = sm.OLS(self.for_wave1[i], X).fit()
            # print(re.params)
            x_pred_o = np.linspace(-5, 105, 110)
            x_pred = sm.add_constant(x_pred_o)
            y_pred = re.predict(x_pred)
            prstd, iv_l, iv_u = wls_prediction_std(re, exog=x_pred, alpha=0.05)
            self.cor_ax1.plot(x_pred_o, iv_l, "-.", c=wave_colors[i], alpha=0.3)
            self.cor_ax1.plot(x_pred_o, iv_u, "-.", c=wave_colors[i], alpha=0.3)
            # self.cor_ax1.plot(self.wave1, re.fittedvalues,"-", c=wave_colors[i], alpha=0.5)
            self.cor_ax1.plot(x_pred_o, y_pred, "-", c=wave_colors[i], alpha=0.5)
            self.show_rois.insert("end", "L" + str(i + 1) + ": " + f"{y_pred[5]:.2f}" + ", " + f"{iv_l[5]:.2f}" +
                                  ", " + f"{iv_u[5]:.2f}" + ", " + f"{(iv_u[5] - iv_l[5]) / 2:.2f}" + "\n")

        self.cor_ax1.set_title("Left figure")
        self.cor_ax1.legend()
        self.cor_ax1.set_xlabel("Respiratory motion")
        self.cor_ax1.set_ylabel("Internal marker motion")
        self.cor_ax1.set_xlim(-5, 105)
        self.cor_ax1.set_ylim(-5, 105)

        for i in range(len(self.roi_center2)):
            # print(np.corrcoef(self.wave2, self.for_wave2)[0,1])
            self.cor_ax2.plot(self.wave2, self.for_wave2[i], "o", mfc="None", c=wave_colors[i],
                              label="#" + str(i + 1))
            X = sm.add_constant(self.wave2)
            re = sm.OLS(self.for_wave2[i], X).fit()
            # print(re.params)
            x_pred_o = np.linspace(-5, 105, 110)
            x_pred = sm.add_constant(x_pred_o)
            y_pred = re.predict(x_pred)
            prstd, iv_l, iv_u = wls_prediction_std(re, exog=x_pred, alpha=0.05)
            self.cor_ax2.plot(x_pred_o, iv_l, "-.", c=wave_colors[i], alpha=0.3)
            self.cor_ax2.plot(x_pred_o, iv_u, "-.", c=wave_colors[i], alpha=0.3)
            # self.cor_ax1.plot(self.wave1, re.fittedvalues,"-", c=wave_colors[i], alpha=0.5)
            self.cor_ax2.plot(x_pred_o, y_pred, "-", c=wave_colors[i], alpha=0.5)
            self.show_rois.insert("end", "R" + str(
                i + 1) + ": " + f"{y_pred[5]:.2f}" + ", " + f"{iv_l[5]:.2f}" + ", " + f"{iv_u[5]:.2f}"
                                  + ", " + f"{(iv_u[5] - iv_l[5]) / 2:.2f}" + "\n")

        self.cor_ax2.set_title("Right figure")
        self.cor_ax2.legend()
        self.cor_ax2.set_xlabel("Respiratory motion")
        self.cor_ax2.set_ylabel("Internal marker motion")
        self.cor_ax2.set_xlim(-5, 105)
        self.cor_ax2.set_ylim(-5, 105)

        correlation_root = tk.Tk()
        correlation_root.wm_title("Correlation")
        correlation_frame = tk.Frame(master=correlation_root)
        correlation_frame.pack()
        self.correlation_canvas = FigureCanvasTkAgg(self.cor_fig, correlation_frame)
        self.correlation_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.correlation_canvas, correlation_root)
        self.correlation_canvas.get_tk_widget().pack()
        self.correlation_canvas.draw()
        correlation_root.mainloop()

    def animate(self):
        fig_func = plt.figure(figsize=(8, 4))
        ax = fig_func.add_subplot(121)
        ax.set_title("Vertical beam")
        ax.axis("off")
        ax2 = fig_func.add_subplot(122)
        ax2.axis("off")
        ax2.set_title("Horizontal beam")
        fig_func.tight_layout()

        def update(i):
            if i != 0:
                plt.cla()
                ax.cla()
                ax2.cla()
            ax.axis("off")
            ax2.axis("off")
            ax.set_title("Vertical beam")
            ax2.set_title("Horizontal beam")

            # gif_array_1 = self.array1[i, :, :]
            # g1scale = 255.0 / (self.vmax - self.vmin)
            # gif_array_1 = (gif_array_1 - self.vmin) * g1scale
            # ax.imshow(gif_array_1, cmap="Greys", vmin=self.vmin * g1scale, vmax=self.vmax * g1scale)
            ax.imshow(self.array1[i, :, :], cmap="Greys", vmin=self.vmin, vmax=self.vmax)
            # gif_array_2 = self.array2[i, :, :]
            # g2scale = 255.0 / (self.vmax2 - self.vmin2)
            # gif_array_2 = (gif_array_2 - self.vmin2) * g2scale
            # ax2.imshow(gif_array_2, cmap="Greys", vmin=self.vmin2 * g2scale, vmax=self.vmax2 * g2scale)
            ax2.imshow(self.array2[i, :, :], cmap="Greys", vmin=self.vmin2, vmax=self.vmax2)
            for j in range(len(self.marker_chase)):
                x_y = (self.marker_chase[j][0][i], self.marker_chase[j][1][i])
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[j], fill=False,
                                     lw=2)
                ax.add_patch(rect)
                ax.text(self.marker_chase[j][0][i] + roi_size / 2 + 20, self.marker_chase[j][1][i] + roi_size / 2,
                        str(j + 1), size=15,
                        color=wave_colors[j])
            for k in range(len(self.marker_chase2)):
                x_y2 = (self.marker_chase2[k][0][i], self.marker_chase2[k][1][i])
                rect2 = pat.Rectangle(xy=x_y2, width=roi_size, height=roi_size, edgecolor=wave_colors[k], fill=False,
                                      lw=2)
                ax2.add_patch(rect2)
                ax2.text(self.marker_chase2[k][0][i] + roi_size / 2 + 20, self.marker_chase2[k][1][i] + roi_size / 2,
                         str(k + 1), size=15,
                         color=wave_colors[k])
            ax.set_xlim([self.xlimlow1, self.xlimhigh1])
            ax.set_ylim([self.ylimhigh1, self.ylimlow1])
            ax2.set_xlim([self.xlimlow2, self.xlimhigh2])
            ax2.set_ylim([self.ylimhigh2, self.ylimlow2])
            # ax.text(1.0, 1.0, "Vert.", ha="left", va="top", transform=ax.transAxes)

        ami = animation.FuncAnimation(fig_func, update, frames=len(self.array1), interval=80)

        def save_gif():
            ami.save(data_base_folder + self.id + "/" + str(self.SOPUID[0]) + ".gif", writer="pillow", fps=10, dpi=300)
            if save_mpeg4 == "yes":
                ami.save(data_base_folder + self.id + "/" + str(self.SOPUID[0]) + ".mp4", writer="ffmpeg", fps=10,
                         dpi=300)
            self.update_show_rois("animation saved" + "\n")

        save_gif()

    def open_dicom(self):
        self.calc_exe["state"] = "active"
        self.patient_folder["state"] = "active"
        self.show_animation["state"] = "disable"
        self.save_data_b["state"] = "disable"
        self.roi_clear["state"] = "active"
        self.calculated = False
        self.show_rois.delete('1.0', tk.END)
        self.update_show_rois("1. Open DICOM\n")
        self.update_show_rois("2. Click markers\n")
        self.update_show_rois("3. Calculate!\n")
        self.update_show_rois("---------------\n")
        self.clear_SI_table()
        self.marker_chase, self.marker_chase2 = [], []
        try:
            self.clear_roi()
        except:
            pass

        init_dir = config.get("settings", "dicom_folder")
        if self.current_patient_dicom_folder == []:
            dicom_dir = tk.filedialog.askdirectory(initialdir=init_dir)
            self.current_patient_dicom_folder = dicom_dir
            wave_dicoms = folder_viewer.get_wave_dicoms(dicom_dir + "/")
        else:
            dicom_dir = tk.filedialog.askdirectory(initialdir=self.current_patient_dicom_folder)
            self.current_patient_dicom_folder = dicom_dir
            wave_dicoms = folder_viewer.get_wave_dicoms(dicom_dir + "/")
        dicom_window = tk.Toplevel()
        dicom_select_app = folder_viewer.DicomSelectGui(wave_dicoms, master=dicom_window)
        dicom_select_app.mainloop()
        self.open_files = dicom_select_app.get_chosen_dicoms()
        dicom_window.destroy()
        dicom1 = pydicom.dcmread(self.open_files[0])
        dicom2 = pydicom.dcmread(self.open_files[1])
        if dicom1[0x0008, 0x1010].value == "H-SIM1":
            self.direction_1 = "H"
        else:
            self.direction_1 = "V"
        if dicom2[0x0008, 0x1010].value == "H-SIM1":
            self.direction_2 = "H"
        else:
            self.direction_2 = "V"
        if self.direction_1 == "V":
            self.dicom1 = dicom1
            self.dicom2 = dicom2
        else:
            self.dicom1 = dicom2
            self.dicom2 = dicom1

        self.id = self.dicom1.PatientID
        self.update_show_rois("ID: " + str(self.id) + "\n")
        self.SOPUID = []
        self.study_date = [self.dicom1.StudyDate, self.dicom2.StudyDate]
        # print(self.study_date)
        self.SOPUID.append(self.dicom1[0x0008, 0x0018].value)
        self.SOPUID.append(self.dicom2[0x0008, 0x0018].value)
        self.current_slice = 0

        self.wave1, self.wave_time1, self.wave_raw1, self.wave_time_raw1 = wave_analysis.wave_analysis(self.dicom1)
        self.wave2, self.wave_time2, self.wave_raw2, self.wave_time_raw2 = wave_analysis.wave_analysis(self.dicom2)
        self.plot_wave()
        self.array1 = self.dicom1.pixel_array
        self.array2 = self.dicom2.pixel_array
        # print(self.array1.shape) #(30, 1512, 1512)
        self.xlimlow1 = 0
        self.ylimlow1 = 0
        self.xlimhigh1 = self.array1.shape[1]
        self.ylimhigh1 = self.array1.shape[2]
        self.xlimlow2 = 0
        self.ylimlow2 = 0
        self.xlimhigh2 = self.array2.shape[1]
        self.ylimhigh2 = self.array2.shape[2]
        try:
            dummy = self.vmin
        except AttributeError:
            self.vmin = 500
            self.vmax = 7000
            self.vmin2 = 500
            self.vmax2 = 7000

        self.plot_image1(self.array1)
        self.plot_image2(self.array2)
        self.pixel_spacing1 = np.array([float(self.dicom1.PixelSpacing[0]),
                                        float(self.dicom1.PixelSpacing[1])]) * spacing_correction_factor
        self.pixel_spacing2 = np.array([float(self.dicom2.PixelSpacing[0]),
                                        float(self.dicom2.PixelSpacing[1])]) * spacing_correction_factor
        self.actime = [self.dicom1.AcquisitionTime, self.dicom2.AcquisitionTime]

    def plot_image1(self, array):
        if self.imshow1 is None:
            self.ax1.cla()
            self.imshow1 = self.ax1.imshow(array[self.current_slice, :, :], vmin=self.vmin, vmax=self.vmax,
                                           interpolation='none', cmap="Greys")
        else:
            self.imshow1.set_data(array[self.current_slice, :, :])
            self.imshow1.set_clim(vmin=self.vmin, vmax=self.vmax)
            self.ax1.patches = []
            self.ax1.texts = []

        self.ax1.axis("off")
        self.fig1.tight_layout()

        if self.roi_center1 and self.current_slice == 0:
            roinum = 1
            for i, point in enumerate(self.roi_center1):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax1.add_patch(rect)
                self.ax1.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
                              color=wave_colors[i])
                roinum += 1
        elif self.marker_chase:
            roinum = 1
            for i, point in enumerate(self.marker_chase):
                x_y = (point[0][self.current_slice], point[1][self.current_slice])
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax1.add_patch(rect)
                self.ax1.text(point[0][self.current_slice] + 40, point[1][self.current_slice] + 30, str(roinum),
                              size=12,
                              color=wave_colors[i])
                roinum += 1
        self.ax1.set_xlim([self.xlimlow1, self.xlimhigh1])
        self.ax1.set_ylim([self.ylimhigh1, self.ylimlow1])

        self.canvas1.draw()


    def plot_image2(self, array):
        if self.imshow2 is None:
            self.ax2.cla()
            self.imshow2 = self.ax2.imshow(array[self.current_slice, :, :], vmin=self.vmin2, vmax=self.vmax2,
                                           interpolation='none',
                                           cmap="Greys")
        else:
            self.imshow2.set_data(array[self.current_slice, :, :])
            self.imshow2.set_clim(vmin=self.vmin2, vmax=self.vmax2)
            self.ax2.patches = []
            self.ax2.texts = []

        self.ax2.axis("off")
        self.fig2.tight_layout()
        if self.roi_center2 and self.current_slice == 0:
            roinum = 1
            for i, point in enumerate(self.roi_center2):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax2.add_patch(rect)
                self.ax2.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
                              color=wave_colors[i])
                roinum += 1
        elif self.marker_chase2:
            roinum = 1
            for i, point in enumerate(self.marker_chase2):
                x_y = (point[0][self.current_slice], point[1][self.current_slice])
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax2.add_patch(rect)
                self.ax2.text(point[0][self.current_slice] + 40, point[1][self.current_slice] + 30, str(roinum),
                              size=12,
                              color=wave_colors[i])
                roinum += 1
        self.ax2.set_xlim([self.xlimlow2, self.xlimhigh2])
        self.ax2.set_ylim([self.ylimhigh2, self.ylimlow2])
        self.canvas2.draw()

    def plot_wave(self):
        wave_colors = ["g", "b", "c", "m", "y", "k", "w"]
        self.wave_ax.cla()
        self.wave_ax.set_xlim(0, 5)
        self.wave_ax.set_ylim(-5, 105)
        self.wave_ax.set_xlabel("Time (sec)")
        self.wave_ax.set_ylabel("Phase (%)")
        self.wave_fig.tight_layout()

        self.wave_ax2.cla()
        self.wave_ax2.set_xlim(0, 5)
        self.wave_ax2.set_ylim(-5, 105)
        self.wave_ax2.set_xlabel("Time (sec)")
        self.wave_ax2.set_ylabel("Phase (%)")
        self.wave_fig2.tight_layout()

        self.wave_ax.plot(self.wave_time_raw1, self.wave_raw1, "o", label="Resp", c="r", alpha=0.5)
        self.wave_ax2.plot(self.wave_time_raw2, self.wave_raw2, "o", label="Resp", c="r", alpha=0.5)
        self.wave_ax.plot(self.wave_time_raw1[self.current_slice], self.wave_raw1[self.current_slice], "o",
                          c="blue", mfc="None", markeredgewidth=2)
        self.wave_ax2.plot(self.wave_time_raw2[self.current_slice], self.wave_raw2[self.current_slice], "o",
                           c="blue", mfc="None", markeredgewidth=2)

        wave1_max = max(self.wave_time_raw1)
        wave2_max = max(self.wave_time_raw2)
        self.wave_ax.set_xlim(0, wave1_max)
        self.wave_ax2.set_xlim(0, wave2_max)
        if self.for_wave1_raw:
            for i in range(len(self.for_wave1_raw)):
                self.wave_ax.plot(self.wave_time_raw1, self.for_wave1_raw[i], "-.",
                                  label="#" + str(i + 1), c=wave_colors[i], alpha=0.7)
        if self.for_wave2_raw:
            for i in range(len(self.for_wave2_raw)):
                self.wave_ax2.plot(self.wave_time_raw2, self.for_wave2_raw[i], "-.",
                                   label="#" + str(i + 1), c=wave_colors[i], alpha=0.7)
        self.wave_ax.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=10)
        self.wave_canvas.draw()

        self.wave_ax2.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=10)
        self.wave_canvas2.draw()

    def clear_roi(self):
        self.clear_SI_table()
        self.calculated = False
        self.show_animation["state"] = "disable"
        self.save_data_b["state"] = "disable"
        self.show_rois.delete('1.0', tk.END)
        self.roi_center1 = []
        self.roi_center2 = []
        self.plot_image1(self.array1)
        self.plot_image2(self.array2)
        self.wave_ax.cla()
        self.wave_ax.set_xlim(0, 5)
        self.wave_ax.set_ylim(-5, 105)
        self.wave_ax.set_xlabel("Time (sec)")
        self.wave_ax.set_ylabel("Phase (%)")
        self.wave_canvas.draw()
        self.wave_ax2.cla()
        self.wave_ax2.set_xlim(0, 5)
        self.wave_ax2.set_ylim(-5, 105)
        self.wave_ax2.set_xlabel("Time (sec)")
        self.wave_ax2.set_ylabel("Phase (%)")
        self.wave_canvas2.draw()
        self.marker_chase = []
        self.marker_chase2 = []
        self.for_wave1 = []
        self.for_wave2 = []
        self.for_wave1_raw = []
        self.for_wave2_raw = []
        self.plot_wave()
        self.update_show_rois("1. Open DICOM\n")
        self.update_show_rois("2. Click markers\n")
        self.update_show_rois("3. Calculate!\n")
        self.update_show_rois("---------------\n")
        self.update_show_rois("ID: " + str(self.id) + "\n")

    def rotate_preset(self):
        val = self.preset_txt.get()
        if val in preset_list:
            ind = preset_list.index(val)
            if ind < len(preset_list) -1:
                self.preset_txt.set(preset_list[ind + 1])
                self.preset_loop(preset_list[ind + 1])
            elif ind == len(preset_list) - 1:
                self.preset_txt.set(preset_list[0])
                self.preset_loop(preset_list[0])
        else:
            self.preset_txt.set(preset_list[0])
            self.preset_loop(preset_list[0])

    def update_show_rois(self, roitext):
        self.show_rois.insert(tk.END, roitext)

    def on_key_press(self, event):
        if event.keysym == 'Control_L' or event.keysym == 'Control_R':
            self.ctrl_is_held = True
        if event.keysym == "Shift_L" or event.keysym == "Shift_R":
            self.shift_is_held = True

        if self.shift_is_held:
            if event.keysym == "Right" or event.keysym == "6":
                self.go_right_edge()
            elif event.keysym == "Left" or event.keysym == "4":
                self.go_left_edge()

        elif event.keysym == 'Right' or event.keysym == "6":
            self.wheel_slice(1)
        elif event.keysym == 'Left' or event.keysym == "4":
            self.wheel_slice(-1)
        elif event.keysym == "Down" or event.keysym == "minus":
            diff1, diff2 = 100, 100
            zoom_factor = 0.2

            self.xlimhigh1 = self.xlimhigh1 + diff1 * zoom_factor
            self.xlimlow1 = self.xlimlow1 - diff1 * zoom_factor
            self.ylimhigh1 = self.ylimhigh1 + diff2 * zoom_factor
            self.ylimlow1 = self.ylimlow1 - diff2 * zoom_factor

            self.xlimhigh2 = self.xlimhigh2 + diff1 * zoom_factor
            self.xlimlow2 = self.xlimlow2 - diff1 * zoom_factor
            self.ylimhigh2 = self.ylimhigh2 + diff2 * zoom_factor
            self.ylimlow2 = self.ylimlow2 - diff2 * zoom_factor

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "Up" or event.keysym == "plus":
            diff1, diff2 = 100, 100
            zoom_factor = 0.2

            self.xlimhigh1 = self.xlimhigh1 - diff1 * zoom_factor
            self.xlimlow1 = self.xlimlow1 + diff1 * zoom_factor
            self.ylimhigh1 = self.ylimhigh1 - diff2 * zoom_factor
            self.ylimlow1 = self.ylimlow1 + diff2 * zoom_factor

            self.xlimhigh2 = self.xlimhigh2 - diff1 * zoom_factor
            self.xlimlow2 = self.xlimlow2 + diff1 * zoom_factor
            self.ylimhigh2 = self.ylimhigh2 - diff2 * zoom_factor
            self.ylimlow2 = self.ylimlow2 + diff2 * zoom_factor

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "2":
            # pan up
            # press_x, press_y = self.ctrlpress1
            dx = 0
            dy = 10

            self.xlimlow1 = self.xlimlow1 - dx
            self.xlimhigh1 = self.xlimhigh1 - dx
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "8":
            # pan up
            # press_x, press_y = self.ctrlpress1
            dx = 0
            dy = -10

            self.xlimlow1 = self.xlimlow1 - dx
            self.xlimhigh1 = self.xlimhigh1 - dx
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "7":
            # pan up
            # press_x, press_y = self.ctrlpress1
            dx = - 10
            dy = 0

            self.xlimlow1 = self.xlimlow1 - dx
            self.xlimhigh1 = self.xlimhigh1 - dx
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "9":
            # pan up
            # press_x, press_y = self.ctrlpress1
            dx = 10
            dy = 0

            self.xlimlow1 = self.xlimlow1 - dx
            self.xlimhigh1 = self.xlimhigh1 - dx
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy

            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "1":
            dx = - 10
            dy = 0
            self.xlimlow2 = self.xlimlow2 - dx
            self.xlimhigh2 = self.xlimhigh2 - dx
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "3":
            dx = 10
            dy = 0
            self.xlimlow2 = self.xlimlow2 - dx
            self.xlimhigh2 = self.xlimhigh2 - dx
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

        elif event.keysym == "p":
            self.rotate_preset()

        elif event.keysym == "o":
            self.open_dicom()

        elif event.keysym == "space":
            if self.is_in_loop:
                self.loop_off()
            else:
                self.loop_on()

        elif event.keysym == "s":
            self.save_data()

        elif event.keysym == "q":
            self.go_quit()

        elif event.keysym == "c":
            self.clear_roi()

        elif event.keysym == "r":
            self.report()

    def on_key_release(self, event):
        if event.keysym == 'Control_L' or event.keysym == 'Control_R':
            self.ctrl_is_held = False
        if event.keysym == "Shift_L" or event.keysym == "Shift_R":
            self.shift_is_held = False

    def on_key1(self, event):
        print('you pressed', event.key, event.xdata, event.ydata)

    def mouse_scrolled1(self, event):
        """
        by scrolling the mouse, the visualized x-ray image will zoom in/out.
        :param event: tk mouse scroll event
        :return:
        """
        if self.ctrl_is_held:
            zoom_factor = 0.05
            if event.button == "down":
                diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
                diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
                self.xlimhigh1 = self.xlimhigh1 - diff1 * zoom_factor
                self.xlimlow1 = self.xlimlow1 + diff1 * zoom_factor
                self.ylimhigh1 = self.ylimhigh1 - diff2 * zoom_factor
                self.ylimlow1 = self.ylimlow1 + diff2 * zoom_factor
                # here
                diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
                diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
                self.xlimhigh2 = self.xlimhigh2 - diff1 * zoom_factor
                self.xlimlow2 = self.xlimlow2 + diff1 * zoom_factor
                self.ylimhigh2 = self.ylimhigh2 - diff2 * zoom_factor
                self.ylimlow2 = self.ylimlow2 + diff2 * zoom_factor
                # tohere
            else:
                diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
                diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
                self.xlimhigh1 = self.xlimhigh1 + diff1 * zoom_factor
                self.xlimlow1 = self.xlimlow1 - diff1 * zoom_factor
                self.ylimhigh1 = self.ylimhigh1 + diff2 * zoom_factor
                self.ylimlow1 = self.ylimlow1 - diff2 * zoom_factor
                # here
                diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
                diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
                self.xlimhigh2 = self.xlimhigh2 + diff1 * zoom_factor
                self.xlimlow2 = self.xlimlow2 - diff1 * zoom_factor
                self.ylimhigh2 = self.ylimhigh2 + diff2 * zoom_factor
                self.ylimlow2 = self.ylimlow2 - diff2 * zoom_factor
                # tohere
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)
        # こにスライスを
        else:
            if event.button == "down":
                self.wheel_slice(1)
            else:
                self.wheel_slice(-1)

    def mouse_scrolled2(self, event):
        if self.ctrl_is_held:
            zoom_factor = 0.05
            if event.button == "down":
                diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
                diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
                self.xlimhigh1 = self.xlimhigh1 - diff1 * zoom_factor
                self.xlimlow1 = self.xlimlow1 + diff1 * zoom_factor
                self.ylimhigh1 = self.ylimhigh1 - diff2 * zoom_factor
                self.ylimlow1 = self.ylimlow1 + diff2 * zoom_factor

                diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
                diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
                self.xlimhigh2 = self.xlimhigh2 - diff1 * zoom_factor
                self.xlimlow2 = self.xlimlow2 + diff1 * zoom_factor
                self.ylimhigh2 = self.ylimhigh2 - diff2 * zoom_factor
                self.ylimlow2 = self.ylimlow2 + diff2 * zoom_factor
            else:
                diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
                diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
                self.xlimhigh1 = self.xlimhigh1 + diff1 * zoom_factor
                self.xlimlow1 = self.xlimlow1 - diff1 * zoom_factor
                self.ylimhigh1 = self.ylimhigh1 + diff2 * zoom_factor
                self.ylimlow1 = self.ylimlow1 - diff2 * zoom_factor
                diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
                diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
                self.xlimhigh2 = self.xlimhigh2 + diff1 * zoom_factor
                self.xlimlow2 = self.xlimlow2 - diff1 * zoom_factor
                self.ylimhigh2 = self.ylimhigh2 + diff2 * zoom_factor
                self.ylimlow2 = self.ylimlow2 - diff2 * zoom_factor
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)
        else:
            if event.button == "down":
                self.wheel_slice(1)
            else:
                self.wheel_slice(-1)

    def wheel_slice(self, num):
        if self.slice_num == 0:
            pass
        elif self.slice_num == 30:
            pass
        else:
            current = self.slice_num.get()
            self.slice_num.set(current + int(num))
        self.draw_plot()

    def onclick1(self, event):
        if event.button == 3:
            # これは右クリック
            self.rpress1 = event.xdata, event.ydata
        elif self.ctrl_is_held:
            self.ctrlpress1 = event.xdata, event.ydata
        elif self.current_slice == 0 and not self.calculated:
            self.update_show_rois(f"{event.xdata:.2f}" + ", " + f"{event.ydata:.2f}" + "\n")
            self.roi_center1.append([event.xdata, event.ydata])
            self.plot_image1(self.array1)
        elif self.calculated:
            close_points = []
            for i, marker in enumerate(self.marker_chase):
                if np.sqrt(np.power(marker[0][self.current_slice] - (event.xdata - roi_size / 2), 2) + np.power(
                        marker[1][self.current_slice] - (event.ydata - roi_size / 2), 2)) < roi_size:
                    close_points.append([i, marker[0][self.current_slice], marker[1][self.current_slice]])
            if len(close_points) == 1:
                self.closest_point = [close_points[0], event.xdata - roi_size / 2, event.ydata - roi_size / 2]
            elif len(close_points) > 1:
                diams = []
                for i in close_points:
                    diams.append(np.sqrt(
                        np.power(i[1] - (event.xdata - roi_size / 2), 2) + np.power(i[2] - (event.ydata - roi_size / 2),
                                                                                    2)))
                # print(np.argmin(diams))
                # print(close_points[int(np.argmin(diams))][0])
                self.closest_point = [close_points[int(np.argmin(diams))], event.xdata - roi_size / 2,
                                      event.ydata - roi_size / 2]

        else:
            msgrt = tk.Tk()
            msgrt.withdraw()
            res = messagebox.showinfo("info", "choose ROI in the first slice")

    def onclick2(self, event):
        if event.button == 3:
            self.rpress2 = event.xdata, event.ydata
        elif self.ctrl_is_held:
            self.ctrlpress2 = event.xdata, event.ydata
        elif self.current_slice == 0 and not self.calculated:
            self.update_show_rois(f"{event.xdata:.2f}" + ", " + f"{event.ydata:.2f}" + "\n")
            self.roi_center2.append([event.xdata, event.ydata])
            self.plot_image2(self.array2)
        elif self.calculated:
            close_points = []
            for i, marker in enumerate(self.marker_chase2):
                if np.sqrt(np.power(marker[0][self.current_slice] - (event.xdata - roi_size / 2), 2) + np.power(
                        marker[1][self.current_slice] - (event.ydata - roi_size / 2), 2)) < roi_size:
                    close_points.append([i, marker[0][self.current_slice], marker[1][self.current_slice]])
            if len(close_points) == 1:
                self.closest_point2 = [close_points[0], event.xdata - roi_size / 2, event.ydata - roi_size / 2]
            elif len(close_points) > 1:
                diams = []
                for i in close_points:
                    diams.append(np.sqrt(
                        np.power(i[1] - (event.xdata - roi_size / 2), 2) + np.power(i[2] - (event.ydata - roi_size / 2),
                                                                                    2)))
                # print(np.argmin(diams))
                # print(close_points[int(np.argmin(diams))][0])
                self.closest_point2 = [close_points[int(np.argmin(diams))], event.xdata - roi_size / 2,
                                       event.ydata - roi_size / 2]
        else:
            msgrt = tk.Tk()
            msgrt.withdraw()
            res = messagebox.showinfo("info", "choose ROI in the first slice")

    def on_motion1(self, event):
        if self.ctrlpress1 == None:
            if self.rpress1 == None:
                if len(self.closest_point) > 0:
                    marker_info, click_x, click_y = self.closest_point
                    marker_num = marker_info[0]
                    # print(event.xdata)
                    dx = event.xdata - roi_size / 2 - click_x
                    dy = event.ydata - roi_size / 2 - click_y
                    newx = self.marker_chase[marker_num][0][self.current_slice] + dx
                    newy = self.marker_chase[marker_num][1][self.current_slice] + dy
                    self.marker_chase[marker_num][0][self.current_slice] = newx
                    self.marker_chase[marker_num][1][self.current_slice] = newy
                    self.closest_point = [marker_info, newx, newy]
                    for j, wave in enumerate(self.for_wave1):
                        calcd_wave = []
                        x_shifts, y_shifts = self.marker_chase[j]
                        calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                        self.for_wave1_raw[j] = calcd_wave
                        calcd_wave = np.delete(calcd_wave, 0)
                        self.for_wave1[j] = calcd_wave
                    self.draw_plot()
                    self.change_SI_table()
            else:
                rpress_x, rpress_y = self.rpress1
                dx = (event.xdata - rpress_x) * ct_val_scale
                dy = (event.ydata - rpress_y) * ct_val_scale
                current_x, current_y = self.x_v.get(), self.y_v.get()
                # x_v: min, y_v: max
                new_min = current_x - dx - dy
                new_max = current_y + dx - dy
                if new_min < ct_val_min:
                    new_min = ct_val_min
                if new_max > ct_val_max:
                    new_max = ct_val_max
                if new_min > new_max:
                    new_min = new_max - 100
                self.x_v.set(int(new_min))
                self.y_v.set(int(new_max))
                self.draw_plot()
        else:
            press_x, press_y = self.ctrlpress1
            dx = event.xdata - press_x
            dy = event.ydata - press_y
            self.xlimlow1 = self.xlimlow1 - dx
            self.xlimhigh1 = self.xlimhigh1 - dx
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

    def on_motion2(self, event):
        if self.ctrlpress2 == None:
            if self.rpress2 == None:
                if len(self.closest_point2) > 0:
                    marker_info, click_x, click_y = self.closest_point2
                    marker_num = marker_info[0]
                    dx = event.xdata - roi_size / 2 - click_x
                    dy = event.ydata - roi_size / 2 - click_y
                    newx = self.marker_chase2[marker_num][0][self.current_slice] + dx
                    newy = self.marker_chase2[marker_num][1][self.current_slice] + dy
                    self.marker_chase2[marker_num][0][self.current_slice] = newx
                    self.marker_chase2[marker_num][1][self.current_slice] = newy
                    self.closest_point2 = [marker_info, newx, newy]
                    for j, wave in enumerate(self.for_wave2):
                        calcd_wave = []
                        x_shifts, y_shifts = self.marker_chase2[j]
                        calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                        self.for_wave2_raw[j] = calcd_wave
                        calcd_wave = np.delete(calcd_wave, 0)
                        self.for_wave2[j] = calcd_wave
                    self.draw_plot()
                    self.change_SI_table()
            else:
                rpress_x, rpress_y = self.rpress2
                dx = (event.xdata - rpress_x) * ct_val_scale
                dy = (event.ydata - rpress_y) * ct_val_scale
                current_x, current_y = self.x_v2.get(), self.y_v2.get()
                # x_v: min, y_v: max
                new_min = current_x - dx - dy
                new_max = current_y + dx - dy
                if new_min < ct_val_min:
                    new_min = ct_val_min
                if new_max > ct_val_max:
                    new_max = ct_val_max
                if new_min > new_max:
                    new_min = new_max - 100
                self.x_v2.set(int(new_min))
                self.y_v2.set(int(new_max))
                self.draw_plot()
        else:
            press_x, press_y = self.ctrlpress2
            dx = event.xdata - press_x
            dy = event.ydata - press_y
            self.xlimlow2 = self.xlimlow2 - dx
            self.xlimhigh2 = self.xlimhigh2 - dx
            self.ylimlow2 = self.ylimlow2 - dy
            self.ylimhigh2 = self.ylimhigh2 - dy
            self.ylimlow1 = self.ylimlow1 - dy
            self.ylimhigh1 = self.ylimhigh1 - dy
            self.plot_image1(self.array1)
            self.plot_image2(self.array2)

    def on_release1(self, event):
        self.ctrlpress1 = None
        self.rpress1 = None
        self.closest_point = []

    def on_release2(self, event):
        self.ctrlpress2 = None
        self.rpress2 = None
        self.closest_point2 = []

    def start_up(self):
        self.x_v.set(500)
        self.y_v.set(7000)
        self.x_v2.set(500)
        self.y_v2.set(7000)
        self.slice_num.set(0)

    def draw_plot(self, event=None):
        self.vmax = self.y_v.get()
        self.vmin = self.x_v.get()
        if self.vmin > self.vmax:
            self.vmin = self.vmax - 10
            self.x_v.set(self.vmax - 10)
        self.vmax2 = self.y_v2.get()
        self.vmin2 = self.x_v2.get()
        if self.vmin2 > self.vmax2:
            self.vmin2 = self.vmax2 - 10
            self.x_v2.set(self.vmax2 - 10)
        self.current_slice = self.slice_num.get()
        self.plot_image1(self.array1)
        self.plot_image2(self.array2)
        self.plot_wave()
        # self.change_SI_table()

    def init_draw(self):
        global h
        self.fig1 = Figure(figsize=(4, 4), dpi=100, facecolor='black')
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.axis("off")

        self.fig2 = Figure(figsize=(4, 4), dpi=100, facecolor='black')
        self.ax2 = self.fig2.add_subplot(111)
        self.ax2.axis("off")

        self.wave_fig = Figure(figsize=(4, 3), dpi=100)
        self.wave_ax = self.wave_fig.add_subplot(111)
        self.wave_ax.set_xlim(0, 5)
        self.wave_ax.set_ylim(-5, 105)
        self.wave_ax.set_xlabel("Time (sec)")
        self.wave_ax.set_ylabel("Phase (%)")
        self.wave_ax.plot([], [])

        self.wave_fig2 = Figure(figsize=(4, 3), dpi=100)
        self.wave_ax2 = self.wave_fig2.add_subplot(111)
        self.wave_ax2.set_xlim(0, 5)
        self.wave_ax2.set_ylim(-5, 105)
        self.wave_ax2.set_xlabel("Time (sec)")
        self.wave_ax2.set_ylabel("Phase (%)")
        self.wave_ax2.plot([], [])


if __name__ == '__main__':
    root = tk.Tk()
    app = Application(master=root)
    app.mainloop()
