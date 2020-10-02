"""
Gather all figs, numpy arrays, and create HTML report file.
"""

import numpy as np
import glob
import shutil
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import os

wave_colors = ["g", "b", "c", "m", "y", "k", "w"]
threshold = 0.6


class build_report:
    def __init__(self, patient_folder):
        self.folder = patient_folder
        self.file_list = glob.glob(self.folder + "/*.npz")
        self.gif_list =[]
        self.gif_list = glob.glob(self.folder + "/*.gif")
        self.fig_list = glob.glob(self.folder + "/*fig1.png")
        self.report_file = self.folder + "/report.html"
        shutil.copyfile("./temp/template.html", self.report_file)
        print("number of analyzed files: ", len(self.file_list))
        print("number of figure1 : ", len(self.fig_list))
        print(patient_folder)
        self.table_file = self.folder + "/table1.png"
        self.resp_file = self.folder + "/resp.png"
        self.corr_file = self.folder + "/corr.png"
        self.order_to_order()
        self.table1, self.uids = self.create_table()
        self.colour_table = self.return_colour_table(self.table1)
        # self.plot_table_fig()
        self.plot_table_dyn()
        self.create_fig3_dyn()
        # self.create_fig4()
        self.create_fig4_dyn()
        self.create_html()

    def order_to_order(self):
        """
        change order of dicom files; old -> new
        :return:
        """
        dicom_order_dict = {}
        study_date=[]
        for npfile in self.file_list:
            study_date.append(np.load(npfile)["study_date"][0])
            for fig in self.fig_list:
                if np.load(npfile)["UID"][0] in fig:
                    dicom_order_dict[npfile] = [np.load(npfile)["actime"][0], np.load(npfile)["UID"], fig,
                                                np.load(npfile)["actime"]]
                    break
                elif np.load(npfile)["UID"][1] in fig:
                    dicom_order_dict[npfile] = [np.load(npfile)["actime"][0], np.load(npfile)["UID"], fig,
                                                np.load(npfile)["actime"]]
                    break
        files_sorted = sorted(dicom_order_dict.items(), key=lambda x: x[1], reverse=False)
        new_nps = []
        new_figs = []
        new_times = []
        new_uids = []
        for a_file in files_sorted:
            new_nps.append(a_file[0])
            new_figs.append(a_file[1][2])
            new_times.append(a_file[1][3])
            new_uids.append(a_file[1][1])
        self.study_date = study_date
        self.file_list = new_nps
        self.fig_list = new_figs
        self.new_times = new_times

    def create_table(self):
        """
        Create table using matplotlib.
        table include marker motion for each dimension.
        :return:
        """
        uids = []
        #table1 = np.zeros(4 * len(self.file_list) * 7).reshape((4, len(self.file_list) + 2, 5))
        table1 = np.zeros(4 * 5 * 7).reshape((4, 5 + 2, 5))
        time1s = []
        x1s = []
        y1s = []
        time2s = []
        x2s = []
        y2s = []
        for i in range(len(self.file_list)):
            opened = np.load(self.file_list[i])
            time1s.append(opened["time1"])
            time2s.append(opened["time2"])
            x1s.append(opened["x1"])
            x2s.append(opened["x2"])
            y1s.append(opened["y1"])
            y2s.append(opened["y2"])
            uids = np.append(uids, opened["UID"])

            for j in range(len(opened["max_x1"])):
                table1[j, i, 0] = opened["max_x1"][j]

            for j in range(len(opened["max_y1"])):
                table1[j, i, 1] = opened["max_y1"][j]

            for j in range(len(opened["max_x2"])):
                table1[j, i, 2] = opened["max_x2"][j]

            for j in range(len(opened["max_y2"])):
                table1[j, i, 3] = opened["max_y2"][j]


            for j in range(len(table1)):
                table1[j, i, 4] = np.sqrt(
                    table1[j, i, 0] ** 2 + table1[j, i, 2] ** 2 + max(table1[j, i, 1], table1[j, i, 3]) ** 2)

        for mai in range(table1.shape[0]):
            for yoko in range(table1.shape[2]):
                table1[mai, 5, yoko] = np.mean(table1[mai, :5, yoko])
                table1[mai, 6, yoko] = np.std(table1[mai, :5, yoko], ddof=1)
        table1 = np.round(table1,2)

        self.time1s = time1s
        self.x1s = x1s
        self.y1s = y1s
        self.time2s = time2s
        self.x2s = x2s
        self.y2s = y2s

        print("Number of unique UIDs: ", len(np.unique(uids)))
        if len(np.unique(uids)) == 10:
            print("ok")
        else:
            print("Calculate all data before print report.")
        return table1, uids

    def return_colour_table(self, table1):
        """
        put yellow color if directional motion > 5 mm
        :param table1: marker motion table
        :return:
        """
        colour_table = np.zeros(len(table1[:, 0, 0]) * len(table1[0, :, 0]) * len(table1[0, 0, :])).reshape(
            len(table1[:, 0, 0]), len(table1[0, :, 0]), len(table1[0, 0, :]))
        colour_table = colour_table.astype(str)
        for k in range(len(table1[:, 0, 0])):
            for i in range(len(table1[0, :, 0])):
                for j in range(len(table1[0, 0, :])):
                    if i >= 5 or j >= 4:
                        colour_table[k, i, j] = "white"
                        continue
                    if table1[k, i, j] > 5.0:
                        colour_table[k, i, j] = "yellow"
                    else:
                        colour_table[k, i, j] = "white"
        return colour_table

    def plot_table_dyn(self):
        """
        plot table dynamically.
        :return:
        """
        table_fig=plt.figure(figsize=(8,4),dpi=100)
        for i in range(self.table1.shape[0]):
            v = i+1
            col_labels = ["x_L", "SI_L", "x_R", "SI_R", "3D"]
            row_labels = ["1", "2", "3", "4", "5", "mean", "std"]
            ax1 = table_fig.add_subplot(2,2,v)
            t1 = ax1.table(cellText=self.table1[i, :, :], colLabels=col_labels,
                                 rowLabels=row_labels, loc="center", cellColours=self.colour_table[i])
            t1.auto_set_font_size(False)
            t1.set_fontsize(10)
            t1.scale(1, 1)
            ax1.set_title("Marker motion #" + str(v) + " (mm)", color=wave_colors[i])
            ax1.set_axis_off()
            plt.tight_layout()
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

            for pos in ['right', 'top', 'bottom', 'left']:
                plt.gca().spines[pos].set_visible(False)
        table_fig.savefig(self.table_file)

    def plot_table_fig(self):
        """
        plot all tables. This is not used.
        :return:
        """
        table_fig = plt.figure(figsize=(8, 4), dpi=100)
        plt.subplots_adjust(hspace=0.3)
        table_ax1 = table_fig.add_subplot(221)
        table_ax1.set_title("Marker motion #1 (mm)", color=wave_colors[0])
        table_ax1.set_axis_off()
        table_ax2 = table_fig.add_subplot(222)
        table_ax2.set_title("Marker #2", color=wave_colors[1])
        table_ax2.set_axis_off()
        table_ax3 = table_fig.add_subplot(223)
        table_ax3.set_title("Marker #3", color=wave_colors[2])
        table_ax3.set_axis_off()
        table_ax4 = table_fig.add_subplot(224)
        table_ax4.set_title("Marker #4", color=wave_colors[3])
        table_ax4.set_axis_off()

        col_labels = ["x_L", "SI_L", "x_R", "SI_R", "3D"]
        row_labels = ["1", "2", "3", "4", "5", "mean", "std"]
        colour_table = self.colour_table
        t1 = table_ax1.table(cellText=self.table1[0, :, :], colLabels=col_labels,
                             rowLabels=row_labels, loc="center", cellColours=colour_table[0])
        t2 = table_ax2.table(cellText=self.table1[1, :, :], colLabels=col_labels,
                             rowLabels=row_labels, loc="center", cellColours=colour_table[1])
        t3 = table_ax3.table(cellText=self.table1[2, :, :], colLabels=col_labels,
                             rowLabels=row_labels, loc="center", cellColours=colour_table[2])
        t4 = table_ax4.table(cellText=self.table1[3, :, :], colLabels=col_labels,
                             rowLabels=row_labels, loc="center", cellColours=colour_table[3])
        t1.auto_set_font_size(False)
        t2.auto_set_font_size(False)
        t3.auto_set_font_size(False)
        t4.auto_set_font_size(False)
        t1.set_fontsize(10)
        t2.set_fontsize(10)
        t3.set_fontsize(10)
        t4.set_fontsize(10)
        t1.scale(1, 1)
        t2.scale(1, 1)
        t3.scale(1, 1)
        t4.scale(1, 1)

        plt.tight_layout()
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

        for pos in ['right', 'top', 'bottom', 'left']:
            plt.gca().spines[pos].set_visible(False)

        table_fig.savefig(self.table_file)


    def create_fig3_dyn(self):
        fig3 = plt.figure(figsize=(8, 16))
        new_times = self.new_times

        figure_number = len(self.time1s)
        for i in range(figure_number):
            ax1_3 = fig3.add_subplot(figure_number,2,(i+1)*2-1)
            ax1_3.plot(self.time1s[i], self.x1s[i], "o", label="resp.", c="r", alpha=0.5)
            for j in range(np.array(self.y1s).shape[1]):
                ax1_3.plot(self.time1s[i], self.y1s[i][j], "-.", label="#"+str(j+1), c=wave_colors[j])
            ax1_3.set_title("Acq. time: " + str(new_times[i][0]))
            ax1_3.set_ylabel("Resp. Phase")
            ax1_3.set_xlabel("Time (sec.)")
            ax1_3.legend()

            axy1_3 = fig3.add_subplot(figure_number, 2, (i+1)*2)
            axy1_3.plot(self.time2s[i], self.x2s[i], "o", label="resp.", c="r", alpha=0.5)
            for j in range(np.array(self.y2s).shape[1]):
                axy1_3.plot(self.time2s[i], self.y2s[i][j], "-.", label="#" + str(j + 1), c=wave_colors[j])
            axy1_3.set_title("Acq. time: " + new_times[i][1])
            axy1_3.set_ylabel("Resp. Phase")
            axy1_3.set_xlabel("Time (sec.)")
            axy1_3.legend()

        plt.tight_layout()
        fig3.savefig(self.resp_file)


    def create_fig3(self):
        """
        create respiratory graph but not used now.
        :return:
        """

        fig3 = plt.figure(figsize=(8, 16))

        new_times = self.new_times

        ax1_3 = fig3.add_subplot(5, 2, 1)
        ax1_3.plot(self.time1s[0], self.x1s[0], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y1s).shape[1]):
            ax1_3.plot(self.time1s[0], self.y1s[0][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        ax1_3.set_title("Acq. time: " + str(new_times[0][0]))
        ax1_3.set_ylabel("Resp. Phase")
        ax1_3.legend()

        axy1_3 = fig3.add_subplot(5, 2, 2)
        axy1_3.plot(self.time2s[0], self.x2s[0], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y2s).shape[1]):
            axy1_3.plot(self.time2s[0], self.y2s[0][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        axy1_3.set_title("Acq. time: " + new_times[0][1])
        axy1_3.legend()

        ax2_3 = fig3.add_subplot(5, 2, 3)
        ax2_3.plot(self.time1s[1], self.x1s[1], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y1s).shape[1]):
            ax2_3.plot(self.time1s[1], self.y1s[1][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        ax2_3.set_title("Acq. time: " + str(new_times[1][0]))
        ax2_3.set_ylabel("Resp. Phase")
        ax2_3.legend()

        axy2_3 = fig3.add_subplot(5, 2, 4)
        axy2_3.plot(self.time2s[1], self.x2s[1], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y2s).shape[1]):
            axy2_3.plot(self.time2s[1], self.y2s[1][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        axy2_3.set_title("Acq. time: " + new_times[1][1])
        axy2_3.legend()

        ax3_3 = fig3.add_subplot(5, 2, 5)
        ax3_3.plot(self.time1s[2], self.x1s[2], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y1s).shape[1]):
            ax3_3.plot(self.time1s[2], self.y1s[2][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        ax3_3.set_title("Acq. time: " + str(new_times[2][0]))
        ax3_3.set_ylabel("Resp. Phase")
        ax3_3.legend()

        axy3_3 = fig3.add_subplot(5, 2, 6)
        axy3_3.plot(self.time2s[2], self.x2s[2], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y2s).shape[1]):
            axy3_3.plot(self.time2s[2], self.y2s[2][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        axy3_3.set_title("Acq. time: " + new_times[2][1])
        axy3_3.legend()

        ax4_3 = fig3.add_subplot(5, 2, 7)
        ax4_3.plot(self.time1s[3], self.x1s[3], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y1s).shape[1]):
            ax4_3.plot(self.time1s[3], self.y1s[3][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        ax4_3.set_title("Acq. time: " + str(new_times[3][0]))
        ax4_3.set_ylabel("Resp. Phase")
        ax4_3.legend()

        axy4_3 = fig3.add_subplot(5, 2, 8)
        axy4_3.plot(self.time2s[3], self.x2s[3], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y2s).shape[1]):
            axy4_3.plot(self.time2s[3], self.y2s[3][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        axy4_3.set_title("Acq. time: " + new_times[3][1])
        axy4_3.legend()

        ax5_3 = fig3.add_subplot(5, 2, 9)
        ax5_3.plot(self.time1s[4], self.x1s[4], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y1s).shape[1]):
            ax5_3.plot(self.time1s[4], self.y1s[4][j], "-.", label="#"+str(j+1), c=wave_colors[j])
        ax5_3.set_title("Acq. time: " + str(new_times[4][0]))
        ax5_3.legend()
        ax5_3.set_xlabel("time (sec.)")
        ax5_3.set_ylabel("Resp. Phase")

        axy5_3 = fig3.add_subplot(5, 2, 10)
        axy5_3.plot(self.time2s[4], self.x2s[4], "o", label="resp.", c="r", alpha=0.5)
        for j in range(np.array(self.y2s).shape[1]):
            axy5_3.plot(self.time2s[4], self.y2s[4][j], "-.", label="#"+str(j+1), c=wave_colors[j])

        axy5_3.set_title("Acq. time: " + new_times[4][1])
        axy5_3.legend()
        axy5_3.set_xlabel("time (sec.)")

        plt.tight_layout()
        fig3.savefig(self.resp_file)

    def create_fig4_dyn(self):
        """
        create correlation plot dynamically.
        :return:
        """
        x1_fig = []
        x1s = np.array(self.x1s)
        y1s = np.array(self.y1s)
        x2s = np.array(self.x2s)
        y2s = np.array(self.y2s)

        """
        x1_2 = np.concatenate([x1s, x2s], 1)
        y1_2 = np.concatenate([y1s, y2s], 2)
        x_1_2_flat = np.ravel(x1_2)
        y_1_2_2d = []
        for j in range(y1_2.shape[1]): #マーカーの数だけ
            some_y = y1_2[:, j, :]
            y_1_2_2d.append(np.ravel(some_y))
        """
        newx1s, newy1s = [],[]
        for j in range(y1s.shape[1]): #マーカー数
            x1_for_a_marker = []
            y1_for_a_marker = []
            for i in range(x1s.shape[0]): #測定繰り返し数
                if np.corrcoef(x1s[i,:],y1s[i,j,:])[0,1] > threshold:
                    x1_for_a_marker = np.concatenate([x1_for_a_marker, x1s[i,:]])
                    y1_for_a_marker = np.concatenate([y1_for_a_marker, y1s[i,j,:]])
                if np.corrcoef(x2s[i,:],y2s[i,j,:])[0,1] > threshold:
                    x1_for_a_marker = np.concatenate([x1_for_a_marker, x2s[i, :]])
                    y1_for_a_marker = np.concatenate([y1_for_a_marker, y2s[i, j, :]])
            newx1s.append(x1_for_a_marker)
            newy1s.append(y1_for_a_marker)


        x_1_2_flat = newx1s
        y_1_2_2d = newy1s

        fig4 = plt.figure(figsize=(8,8))
        for i in range(len(y_1_2_2d)): #iはマーカーの数
            v=i+1
            ax1 = fig4.add_subplot(2, 2, v)
            ax1.plot(x_1_2_flat[i], y_1_2_2d[i], "o", c=wave_colors[i], mfc="None", alpha=0.5)
            ax1.set_xlabel("Resp. phase (%)")
            ax1.set_ylabel("Marker phase (%)")
            ax1.set_title("Marker #"+str(i+1),c=wave_colors[i])
            X1 = sm.add_constant(x_1_2_flat[i])
            re1 = sm.OLS(y_1_2_2d[i], X1).fit()
            x1_pred_o = np.linspace(-5, 105, 110)
            x1_pred = sm.add_constant(x1_pred_o)
            y1_pred = re1.predict(x1_pred)
            prstd, iv_l, iv_u = wls_prediction_std(re1, exog=x1_pred, alpha=0.05)
            ax1.plot(x1_pred_o, iv_l, "-.", c=wave_colors[i], alpha=0.3)
            ax1.plot(x1_pred_o, iv_u, "-.", c=wave_colors[i], alpha=0.3)
            ax1.plot(x1_pred_o, y1_pred, "-", c=wave_colors[i], alpha=0.5)
            ax1.set_xlim(-5, 105)
            ax1.set_ylim(-5, 105)
            ax1.text(65, 0, "fit   : " + str(np.round(y1_pred[5], 1)) + "\nlowr: " + str(
                np.round(iv_l[5], 1)) + "\nupr : " + str(np.round(iv_u[5], 1)) + "\npred.: " + str(
                np.round((iv_u[5] - iv_l[5]) / 2, 1)) + "(%)", size=10, color="black")

        plt.tight_layout()
        fig4.savefig(self.corr_file)


    def create_fig4(self):
        """
        Create correlation plot but not used.
        :return:
        """
        x1_fig = []
        x2_fig = []
        y1_fig = []
        y2_fig = []
        x1s = np.array(self.x1s)
        y1s = np.array(self.y1s)
        x2s = np.array(self.x2s)
        y2s = np.array(self.y2s)
        m1 = []
        m2 = []
        m3 = []
        m4 = []

        for i in range(len(x1s)):
            x1_fig = np.concatenate([x1_fig, x1s[i, :]], 0)
            x1_fig = np.concatenate([x1_fig, x2s[i, :]], 0)
            m1 = np.concatenate([m1, y1s[i, 0, :]], 0)
            m1 = np.concatenate([m1, y2s[i, 0, :]], 0)
            m2 = np.concatenate([m2, y1s[i, 1, :]], 0)
            m2 = np.concatenate([m2, y2s[i, 1, :]], 0)
            m3 = np.concatenate([m3, y1s[i, 2, :]], 0)
            m3 = np.concatenate([m3, y2s[i, 2, :]], 0)
            m4 = np.concatenate([m4, y1s[i, 3, :]], 0)
            m4 = np.concatenate([m4, y2s[i, 3, :]], 0)

        fig4 = plt.figure(figsize=(8, 8))
        ax1_4 = fig4.add_subplot(2, 2, 1)
        ax2_4 = fig4.add_subplot(2, 2, 2)
        ax3_4 = fig4.add_subplot(2, 2, 3)
        ax4_4 = fig4.add_subplot(2, 2, 4)

        ax1_4.plot(x1_fig, m1, "o", c=wave_colors[0], mfc="None", alpha=0.5)
        ax2_4.plot(x1_fig, m2, "o", c=wave_colors[1], mfc="None", alpha=0.5)
        ax3_4.plot(x1_fig, m3, "o", c=wave_colors[2], mfc="None", alpha=0.5)
        ax4_4.plot(x1_fig, m4, "o", c=wave_colors[3], mfc="None", alpha=0.5)

        ax3_4.set_xlabel("Resp. phase (%)")
        ax4_4.set_xlabel("Resp. phase (%)")
        ax1_4.set_ylabel("Marker phase (%)")
        ax3_4.set_ylabel("Marker phase (%)")

        X1 = sm.add_constant(x1_fig)
        re1 = sm.OLS(m1, X1).fit()
        x1_pred_o = np.linspace(-5, 105, 110)
        x1_pred = sm.add_constant(x1_pred_o)
        y1_pred = re1.predict(x1_pred)
        prstd, iv_l, iv_u = wls_prediction_std(re1, exog=x1_pred, alpha=0.05)
        ax1_4.plot(x1_pred_o, iv_l, "-.", c=wave_colors[0], alpha=0.3)
        ax1_4.plot(x1_pred_o, iv_u, "-.", c=wave_colors[0], alpha=0.3)
        ax1_4.plot(x1_pred_o, y1_pred, "-", c=wave_colors[0], alpha=0.5)
        ax1_4.set_xlim(-5, 105)
        ax1_4.set_ylim(-5, 105)
        ax1_4.set_title("Marker #1", c=wave_colors[0])
        ax1_4.text(65, 0, "fit   : " + str(np.round(y1_pred[5], 1)) + "\nlowr: " + str(
            np.round(iv_l[5], 1)) + "\nupr : " + str(np.round(iv_u[5], 1)) + "\npred.: " + str(
            np.round((iv_u[5] - iv_l[5]) / 2, 1)) + "(%)", size=10, color="black")

        re2 = sm.OLS(m2, X1).fit()
        y2_pred = re2.predict(x1_pred)
        prstd, iv_l, iv_u = wls_prediction_std(re2, exog=x1_pred, alpha=0.05)
        ax2_4.plot(x1_pred_o, iv_l, "-.", c=wave_colors[1], alpha=0.3)
        ax2_4.plot(x1_pred_o, iv_u, "-.", c=wave_colors[1], alpha=0.3)
        ax2_4.plot(x1_pred_o, y2_pred, "-", c=wave_colors[1], alpha=0.5)
        ax2_4.set_xlim(-5, 105)
        ax2_4.set_ylim(-5, 105)
        ax2_4.set_title("Marker #2", c=wave_colors[1])
        ax2_4.text(65, 0, "fit   : " + str(np.round(y2_pred[5], 1)) + "\nlowr: " + str(
            np.round(iv_l[5], 1)) + "\nupr : " + str(np.round(iv_u[5], 1)) + "\npred.: " + str(
            np.round((iv_u[5] - iv_l[5]) / 2, 1)) + "(%)", size=10, color="black")

        re3 = sm.OLS(m3, X1).fit()
        y3_pred = re3.predict(x1_pred)
        prstd, iv_l, iv_u = wls_prediction_std(re3, exog=x1_pred, alpha=0.05)
        ax3_4.plot(x1_pred_o, iv_l, "-.", c=wave_colors[2], alpha=0.3)
        ax3_4.plot(x1_pred_o, iv_u, "-.", c=wave_colors[2], alpha=0.3)
        ax3_4.plot(x1_pred_o, y3_pred, "-", c=wave_colors[2], alpha=0.5)
        ax3_4.set_xlim(-5, 105)
        ax3_4.set_ylim(-5, 105)
        ax3_4.set_title("Marker #3", c=wave_colors[2])
        ax3_4.text(65, 0, "fit   : " + str(np.round(y3_pred[5], 1)) + "\nlowr: " + str(
            np.round(iv_l[5], 1)) + "\nupr : " + str(np.round(iv_u[5], 1)) + "\npred.: " + str(
            np.round((iv_u[5] - iv_l[5]) / 2, 1)) + "(%)", size=10, color="black")

        re4 = sm.OLS(m4, X1).fit()
        y4_pred = re4.predict(x1_pred)
        prstd, iv_l, iv_u = wls_prediction_std(re4, exog=x1_pred, alpha=0.05)
        ax4_4.plot(x1_pred_o, iv_l, "-.", c=wave_colors[3], alpha=0.3)
        ax4_4.plot(x1_pred_o, iv_u, "-.", c=wave_colors[3], alpha=0.3)
        ax4_4.plot(x1_pred_o, y4_pred, "-", c=wave_colors[3], alpha=0.5)
        ax4_4.set_xlim(-5, 105)
        ax4_4.set_ylim(-5, 105)
        ax4_4.set_title("Marker #4", c=wave_colors[3])
        ax4_4.text(65, 0, "fit   : " + str(np.round(y4_pred[5], 1)) + "\nlowr: " + str(
            np.round(iv_l[5], 1)) + "\nupr : " + str(np.round(iv_u[5], 1)) + "\npred.: " + str(
            np.round((iv_u[5] - iv_l[5]) / 2, 1)) + "(%)", size=10, color="black")


        plt.tight_layout()
        fig4.savefig(self.corr_file)

    def create_html(self):
        """
        Create html report
        :return:
        """
        with open(self.report_file) as inf:
            txt = inf.read()
            soup = BeautifulSoup(txt, features="lxml")

        tag_pid = soup.new_tag("p")
        tag_pid.string = "Study Date: " + str(self.study_date[0]) + "\n"
        soup.body.append(tag_pid)

        tag_pid = soup.new_tag("p")
        tag_pid.string = "Study Date: " + str(self.study_date[0]) + "\n"
        tag_pid.string = "Patient ID: "+os.path.split(self.folder)[1] +"\n"
        soup.body.append(tag_pid)

        if self.gif_list != []:
            tag_fig1 = soup.new_tag('img', src=self.gif_list[0].split("\\")[-1])
            soup.body.append(tag_fig1)
        else:
            tag_fig1 = soup.new_tag('img', src=self.fig_list[0].split("\\")[-1])
            soup.body.append(tag_fig1)

        tag_table1 = soup.new_tag('img', src=self.table_file.rsplit("/")[-1])
        soup.body.append(tag_table1)

        tag_corr = soup.new_tag("img", src=self.corr_file.rsplit("/")[-1])
        soup.body.append(tag_corr)

        tag_resp = soup.new_tag("img", src=self.resp_file.rsplit("/")[-1])
        soup.body.append(tag_resp)



        with open(self.report_file, "w") as outf:
            outf.write(str(soup))


def main():
    patient_list = glob.glob("./data_base/*")
    patient = build_report(patient_list[0])


if __name__ == '__main__':
    main()
