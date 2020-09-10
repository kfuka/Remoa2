import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import matplotlib.patches as pat
import matplotlib.animation as animation
from tkinter import filedialog
from scipy import signal
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import os
import time

import mil_tracker
import wave_analysis
import folder_viewer
import unite_class

roi_size = 30  # pixel
wave_colors = ["g", "b", "c", "m", "y", "k", "w"]
data_base_folder = "./data_base/"


class Application(tk.Frame):
    def __init__(self, master=None):
        self.roi_center1 = []
        self.roi_center2 = []
        self.for_wave1 = []
        self.for_wave2 = []
        super().__init__(master)
        self.master = master
        self.master.title('Respiratory motion analyzer')
        self.grid()
        self.init_draw()
        self.create_widgets()
        self.start_up()
        self.rpress1 = None
        self.rpress2 = None

    def create_widgets(self):
        self.canvas_frame1 = tk.Frame(self.master)
        self.canvas_frame1.grid(row=0, column=0, rowspan=2)
        self.canvas_frame2 = tk.Frame(self.master)
        self.canvas_frame2.grid(row=0, column=2, rowspan=2)
        self.control_frame = tk.Frame(self.master)
        self.control_frame.grid(row=0, column=1, rowspan=2)
        self.control_frame2 = tk.Frame(self.master)
        self.control_frame2.grid(row=0, column=3, rowspan=2)
        self.wave_frame = tk.Frame(self.master)
        self.wave_frame.grid(row=2, column=0, columnspan=1)
        self.wave_frame2 = tk.Frame(self.master)
        self.wave_frame2.grid(row=2, column=2, columnspan=1)
        self.main_frame = tk.Frame(self.master)
        self.main_frame.grid(row=0, column=4, rowspan=3, sticky=tk.N)

        self.canvas1 = FigureCanvasTkAgg(self.fig1, self.canvas_frame1)
        self.canvas1.draw()
        # self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas1.get_tk_widget().grid(row=0, column=0)
        self.canvas1.mpl_connect('button_press_event', self.onclick1)
        self.canvas1.mpl_connect('scroll_event', self.mouse_scrolled1)
        self.canvas1.mpl_connect('motion_notify_event', self.on_motion1)
        self.canvas1.mpl_connect('button_release_event', self.on_release1)

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
        self.x_scale = tk.Scale(self.control_frame,
                                variable=self.x_v,
                                from_=10000,
                                to=0,
                                resolution=100,
                                orient=tk.VERTICAL,
                                length=500,
                                showvalu=0,
                                command=self.draw_plot)
        self.x_scale.grid(row=0, column=0, padx=0, pady=0, sticky=tk.W + tk.N + tk.S)

        self.x_scale2 = tk.Scale(self.control_frame2,
                                 variable=self.x_v2,
                                 from_=10000,
                                 to=0,
                                 resolution=100,
                                 orient=tk.VERTICAL,
                                 length=500,
                                 showvalu=0,
                                 command=self.draw_plot)
        self.x_scale2.grid(row=0, column=0, padx=0, pady=0, sticky=tk.W + tk.N + tk.S)

        self.y_v = tk.IntVar()
        self.y_v2 = tk.IntVar()
        self.y_scale = tk.Scale(self.control_frame,
                                variable=self.y_v,
                                from_=10000,
                                to=0,
                                resolution=100,
                                orient=tk.VERTICAL,
                                showvalue=0,
                                command=self.draw_plot)
        self.y_scale.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S)

        self.y_scale2 = tk.Scale(self.control_frame2,
                                 variable=self.y_v2,
                                 from_=10000,
                                 to=0,
                                 resolution=100,
                                 orient=tk.VERTICAL,
                                 showvalue=0,
                                 command=self.draw_plot)
        self.y_scale2.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S)

        self.open_button = tk.Button(self.main_frame, text="Open DICOM", command=self.open_dicom)
        self.open_button.grid(row=1, column=0, padx=2, pady=2, rowspan=2, sticky=tk.EW + tk.NS)

        self.show_rois = tk.Text(self.main_frame, height=50, width=40, wrap=tk.CHAR)
        self.show_rois.grid(row=0, column=0, padx=0, pady=0, columnspan=2, sticky=tk.N + tk.EW)
        self.update_show_rois("1. Open DICOM\n")
        self.update_show_rois("2. Click markers\n")
        self.update_show_rois("3. Calculate!\n")
        self.update_show_rois("---------------\n")

        self.roi_clear = tk.Button(self.main_frame, text="ROI Clear", command=self.clear_roi)
        self.roi_clear.grid(row=1, column=1, padx=2, pady=2, sticky=tk.EW)

        self.calc_exe = tk.Button(self.main_frame, text="Calculate!", command=self.calc)
        self.calc_exe.grid(row=3, column=0, padx=2, pady=2, rowspan=2, sticky=tk.EW + tk.NS)

        self.print_report = tk.Button(self.main_frame, text="REPORT", command=self.report)
        self.print_report.grid(row=4, column=1, padx=2, pady=2, sticky=tk.EW + tk.NS)

        self.save_data = tk.Button(self.main_frame, text="Save Data", command=self.save_data)
        self.save_data.grid(row=5, column=0, padx=2, pady=2, sticky=tk.EW + tk.NS)
        self.save_data["state"] = "disable"

        self.show_corr = tk.Button(self.main_frame, text="Correlation", command=self.correlation)
        self.show_corr.grid(row=2, column=1, padx=2, pady=2, sticky=tk.EW)

        self.show_animation = tk.Button(self.main_frame, text="Save Animation", command=self.animate)
        self.show_animation.grid(row=3, column=1, padx=2, pady=2, sticky=tk.EW)

        self.quit_button = tk.Button(self.main_frame, text="quit", command=self.go_quit)
        self.quit_button.grid(row=5, column=1, padx=2, pady=2, sticky=tk.EW + tk.NS)

    def go_quit(self):
        root.quit()
        root.destroy()

    def save_data(self):
        # ファイル名はUnix time
        ut = time.time()
        directory_path = data_base_folder + self.id
        if not os.path.isdir(directory_path):
            os.mkdir(directory_path)
        np.savez(directory_path + "/" + str(self.SOPUID[0]),
                 time1 = self.wave_time1,
                 x1=self.wave1,
                 y1=self.for_wave1,
                 time2 = self.wave_time2,
                 x2=self.wave2,
                 y2=self.for_wave2,
                 UID=self.SOPUID,
                 max_x1 = self.max_x1,
                 max_x2 = self.max_x2,
                 max_y1 = self.max_y1,
                 max_y2 = self.max_y2,
                 actime = self.actime,
                 study_date = self.study_date)
        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
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
        patient_save_folder = tk.filedialog.askdirectory(initialdir="./data_base/")
        patient = unite_class.build_report(patient_save_folder)
        self.update_show_rois("Report created.\n")

    def calc(self):
        max_x1 = []
        max_x2 = []
        max_y1 = []
        max_y2 = []
        if not self.roi_center2:
            print("Select ROI first")
        else:
            self.show_rois.tag_config('over', foreground='red', font='TkDefaultFont 12 bold')
            self.update_show_rois("--------results-----------" + "\n")
            self.update_show_rois("#: x, y max shift" + " in mm\n")
            self.marker_chase = []
            self.marker_chase2 = []
            self.update_show_rois("Left x_1, y_1\n")

            for i in range(len(self.roi_center1)):
                returned = mil_tracker.track_MIL(
                    (self.roi_center1[i][0] - roi_size / 2, self.roi_center1[i][1] - roi_size / 2, roi_size, roi_size),
                    self.array1)
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
                self.show_rois.insert("end", "L" + str(i + 1) + ": ")
                if max_x_move > 5.0:
                    self.show_rois.insert("end", f"{max_x_move:.2f}" + ", ", "over")
                else:
                    self.show_rois.insert("end", f"{max_x_move:.2f}" + ", ")
                if max_y_move > 5.0:
                    self.show_rois.insert("end", f"{max_y_move:.2f}" + "\n", "over")
                else:
                    self.show_rois.insert("end", f"{max_y_move:.2f}" + "\n")
                calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                calcd_wave = np.delete(calcd_wave,0)
                self.for_wave1.append(calcd_wave)
            self.update_show_rois("Right x_2, y_2\n")
            for i in range(len(self.roi_center2)):
                returned = mil_tracker.track_MIL(
                    (self.roi_center2[i][0] - roi_size / 2, self.roi_center2[i][1] - roi_size / 2, roi_size, roi_size),
                    self.array2)
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
                self.update_show_rois("R" + str(i + 1) + ": ")
                if max_x_move > 5.0:
                    self.show_rois.insert("end", f"{max_x_move:.2f}" + ", ", "over")
                else:
                    self.show_rois.insert("end", f"{max_x_move:.2f}" + ", ")
                if max_y_move > 5.0:
                    self.show_rois.insert("end", f"{max_y_move:.2f}" + "\n", "over")
                else:
                    self.show_rois.insert("end", f"{max_y_move:.2f}" + "\n")
                calcd_wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
                calcd_wave = np.delete(calcd_wave,0)
                self.for_wave2.append(calcd_wave)
            self.show_rois.insert("end", "-----------------\n")
            self.show_rois.insert("end", "#: X_1, X_2, SI, 3D norm (mm)\n")
            for i in range(len(max_x1)):
                self.show_rois.insert("end", str(i + 1) + ": ")
                if max_x1[i] > 5.0:
                    self.show_rois.insert("end", f"{max_x1[i]:.2f}" + ", ", "over")
                else:
                    self.show_rois.insert("end", f"{max_x1[i]:.2f}" + ", ")
                if max_x2[i] > 5.0:
                    self.show_rois.insert("end", f"{max_x2[i]:.2f}" + ", ", "over")
                else:
                    self.show_rois.insert("end", f"{max_x2[i]:.2f}" + ", ")
                si_large = max(max_y1[i], max_y2[i])
                if si_large > 5.0:
                    self.show_rois.insert("end", f"{si_large:.2f}" + ", ", "over")
                else:
                    self.show_rois.insert("end", f"{si_large:.2f}" + ", ")
                self.show_rois.insert("end", f"{np.sqrt(max_x1[i] ** 2 + max_x2[i] ** 2 + si_large ** 2):.2f}" + "\n")
            self.show_rois.insert("end", "-----------------\n")
            self.max_x1 = max_x1
            self.max_x2 = max_x2
            self.max_y1 = max_y1
            self.max_y2 = max_y2

        self.plot_wave()

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
        self.save_data["state"] = "active"
        for i in range(len(self.roi_center1)):
            self.cor_ax1.plot(self.wave1, self.for_wave1[i], "o", mfc="None", c=wave_colors[i],
                              label=str(i + 1) + "th marker")
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
            self.cor_ax2.plot(self.wave2, self.for_wave2[i], "o", mfc="None", c=wave_colors[i],
                              label=str(i + 1) + "th marker")
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
        ax.axis("off")
        ax2 = fig_func.add_subplot(122)
        ax2.axis("off")
        fig_func.tight_layout()

        def update(i):
            if i != 0:
                plt.cla()
                ax.cla()
                ax2.cla()
            ax.axis("off")
            ax2.axis("off")
            ax.imshow(self.array1[i, :, :], cmap="Greys", vmin=self.vmin, vmax=self.vmax)
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

        ami = animation.FuncAnimation(fig_func, update, frames=len(self.array1), interval=80)

        def save_gif():
            self.filename = "motion"
            filename = filedialog.asksaveasfilename(initialdir="./", initialfile=self.filename, title="save as",
                                                    filetypes=[("gif file", "*.gif")])
            if filename.endswith(".gif"):
                ami.save(filename, writer='pillow', fps=10, dpi=300)
            else:
                ami.save(filename + ".gif", writer='pillow', fps=10, dpi=300)
            self.update_show_rois("animation saved in gif" + "\n")

        save_gif()

    def open_dicom(self):
        dicom_dir = tk.filedialog.askdirectory(initialdir="../")
        wave_dicoms = folder_viewer.get_wave_dicoms(dicom_dir + "/")
        dicom_window = tk.Toplevel()
        dicom_select_app = folder_viewer.DicomSelectGui(wave_dicoms, master=dicom_window)
        dicom_select_app.mainloop()
        self.open_files = dicom_select_app.get_chosen_dicoms()
        dicom_window.destroy()
        self.dicom1 = pydicom.dcmread(self.open_files[0])
        self.dicom2 = pydicom.dcmread(self.open_files[1])
        self.id = self.dicom1.PatientID
        self.SOPUID = []
        self.study_date = [self.dicom1.StudyDate, self.dicom2.StudyDate]
        # print(self.study_date)
        self.SOPUID.append(self.dicom1[0x0008, 0x0018].value)
        self.SOPUID.append(self.dicom2[0x0008, 0x0018].value)
        self.wave1, self.wave_time1 = wave_analysis.wave_analysis(self.dicom1)
        self.wave2, self.wave_time2 = wave_analysis.wave_analysis(self.dicom2)
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
        self.vmin = np.min(500)
        self.vmax = np.max(7000)
        self.vmin2 = np.min(500)
        self.vmax2 = np.max(7000)
        self.x_v.set(self.vmin)
        self.y_v.set(self.vmax)
        self.plot_image1(self.array1)
        self.plot_image2(self.array2)
        self.pixel_spacing1 = self.dicom1.PixelSpacing
        self.pixel_spacing2 = self.dicom2.PixelSpacing
        self.actime = [self.dicom1.AcquisitionTime, self.dicom2.AcquisitionTime]

    def plot_image1(self, array):
        self.ax1.cla()
        self.ax1.imshow(array[0, :, :], vmin=self.vmin, vmax=self.vmax, cmap="Greys")
        self.ax1.axis("off")
        self.fig1.tight_layout()
        if self.roi_center1:
            roinum = 1
            for i, point in enumerate(self.roi_center1):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax1.add_patch(rect)
                self.ax1.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
                              color=wave_colors[i])
                roinum += 1
        self.ax1.set_xlim([self.xlimlow1, self.xlimhigh1])
        self.ax1.set_ylim([self.ylimhigh1, self.ylimlow1])

        self.canvas1.draw()

    def plot_image2(self, array):
        self.ax2.cla()
        self.ax2.imshow(array[0, :, :], vmin=self.vmin2, vmax=self.vmax2, cmap="Greys")
        self.ax2.axis("off")
        self.fig2.tight_layout()
        if self.roi_center2:
            roinum = 1
            for i, point in enumerate(self.roi_center2):
                x_y = (point[0] - roi_size / 2, point[1] - roi_size / 2)
                rect = pat.Rectangle(xy=x_y, width=roi_size, height=roi_size, edgecolor=wave_colors[i], fill=False,
                                     lw=2)
                self.ax2.add_patch(rect)
                self.ax2.text(point[0] + roi_size / 2 + 10, point[1] + roi_size / 2, str(roinum), size=12,
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

        self.wave_ax.plot(self.wave_time1, self.wave1, "o", label="Resp. data", c="r", alpha=0.5)
        self.wave_ax2.plot(self.wave_time2, self.wave2, "o", label="Resp. data", c="r", alpha=0.5)
        wave1_max = max(self.wave_time1)
        wave2_max = max(self.wave_time2)
        self.wave_ax.set_xlim(0, wave1_max)
        self.wave_ax2.set_xlim(0, wave2_max)
        if self.for_wave1:
            for i in range(len(self.for_wave1)):
                self.wave_ax.plot(self.wave_time1, self.for_wave1[i], "-.",
                                  label="marker" + str(i + 1), c=wave_colors[i], alpha=0.7)
        if self.for_wave2:
            for i in range(len(self.for_wave2)):
                self.wave_ax2.plot(self.wave_time1, self.for_wave2[i], "-.",
                                   label="marker" + str(i + 1), c=wave_colors[i], alpha=0.7)
        self.wave_ax.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=10)
        self.wave_canvas.draw()

        self.wave_ax2.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=10)
        self.wave_canvas2.draw()

    def clear_roi(self):
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
        self.for_wave1 = []
        self.for_wave2 = []
        self.plot_wave()
        self.update_show_rois("1. Open DICOM\n")
        self.update_show_rois("2. Click markers\n")
        self.update_show_rois("3. Calculate!\n")
        self.update_show_rois("---------------\n")

    def update_show_rois(self, roitext):
        self.show_rois.insert(tk.END, roitext)

    def mouse_scrolled1(self, event):
        zoom_factor = 0.05
        if event.button == "up":
            diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
            diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
            self.xlimhigh1 = self.xlimhigh1 - diff1 * zoom_factor
            self.xlimlow1 = self.xlimlow1 + diff1 * zoom_factor
            self.ylimhigh1 = self.ylimhigh1 - diff2 * zoom_factor
            self.ylimlow1 = self.ylimlow1 + diff2 * zoom_factor
            #here
            diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
            diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
            self.xlimhigh2 = self.xlimhigh2 - diff1 * zoom_factor
            self.xlimlow2 = self.xlimlow2 + diff1 * zoom_factor
            self.ylimhigh2 = self.ylimhigh2 - diff2 * zoom_factor
            self.ylimlow2 = self.ylimlow2 + diff2 * zoom_factor
            #tohere
        else:
            diff1 = (self.xlimhigh1 - self.xlimlow1) / 2.0
            diff2 = (self.ylimhigh1 - self.ylimlow1) / 2.0
            self.xlimhigh1 = self.xlimhigh1 + diff1 * zoom_factor
            self.xlimlow1 = self.xlimlow1 - diff1 * zoom_factor
            self.ylimhigh1 = self.ylimhigh1 + diff2 * zoom_factor
            self.ylimlow1 = self.ylimlow1 - diff2 * zoom_factor
            #here
            diff1 = (self.xlimhigh2 - self.xlimlow2) / 2.0
            diff2 = (self.ylimhigh2 - self.ylimlow2) / 2.0
            self.xlimhigh2 = self.xlimhigh2 + diff1 * zoom_factor
            self.xlimlow2 = self.xlimlow2 - diff1 * zoom_factor
            self.ylimhigh2 = self.ylimhigh2 + diff2 * zoom_factor
            self.ylimlow2 = self.ylimlow2 - diff2 * zoom_factor
            #tohere
        self.plot_image1(self.array1)
        self.plot_image2(self.array2)

    def mouse_scrolled2(self, event):
        zoom_factor = 0.05
        if event.button == "up":
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

    def onclick1(self, event):
        if event.button == 3:
            self.rpress1 = event.xdata, event.ydata
        else:
            self.update_show_rois(f"{event.xdata:.2f}" + ", " + f"{event.ydata:.2f}" + "\n")
            self.roi_center1.append([event.xdata, event.ydata])
            self.plot_image1(self.array1)

    def onclick2(self, event):
        if event.button == 3:
            self.rpress2 = event.xdata, event.ydata
        else:
            self.update_show_rois(f"{event.xdata:.2f}" + ", " + f"{event.ydata:.2f}" + "\n")
            self.roi_center2.append([event.xdata, event.ydata])
            self.plot_image2(self.array2)

    def on_motion1(self, event):
        if self.rpress1 == None:
            return
        press_x, press_y = self.rpress1
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
        if self.rpress2 == None:
            return
        press_x, press_y = self.rpress2
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
        self.rpress1 = None

    def on_release2(self, event):
        self.rpress2 = None

    def start_up(self):
        self.x_v.set(500)
        self.y_v.set(7000)
        self.x_v2.set(500)
        self.y_v2.set(7000)

    def draw_plot(self, event=None):

        self.vmin = self.x_v.get()
        self.vmax = self.y_v.get()
        self.vmin2 = self.x_v2.get()
        self.vmax2 = self.y_v2.get()
        self.plot_image1(self.array1)
        self.plot_image2(self.array2)

    def init_draw(self):
        global h
        self.fig1 = Figure(figsize=(5, 5), dpi=100, facecolor='black')
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.axis("off")

        self.fig2 = Figure(figsize=(5, 5), dpi=100, facecolor='black')
        self.ax2 = self.fig2.add_subplot(111)
        self.ax2.axis("off")

        self.wave_fig = Figure(figsize=(5, 3), dpi=100)
        self.wave_ax = self.wave_fig.add_subplot(111)
        self.wave_ax.set_xlim(0, 5)
        self.wave_ax.set_ylim(-5, 105)
        self.wave_ax.set_xlabel("Time (sec)")
        self.wave_ax.set_ylabel("Phase (%)")
        self.wave_ax.plot([], [])

        self.wave_fig2 = Figure(figsize=(5, 3), dpi=100)
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
