"""
Gather all figs, numpy arrays, and create HTML report file.
"""

import configparser
import glob
import os
import shutil
import datetime

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from bs4 import BeautifulSoup
from statsmodels.sandbox.regression.predstd import wls_prediction_std

wave_colors = ["g", "b", "c", "m", "y", "k", "w"]
config = configparser.ConfigParser()
config.read('./config.ini')
# print(np.float(config.get("settings", "threshold_of_corr_coef")))

threshold = np.float(config.get("settings", "threshold_of_corr_coef"))


class build_report:
    def __init__(self, patient_folder):
        self.folder = patient_folder
        self.file_list = glob.glob(self.folder + "/*.npz")
        self.gif_list = []
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
        study_date = []
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
            active_gyo = 0
            for in_gyo in range(5):
                a_gyo = table1[mai, in_gyo, :]
                if not all(a_gyo == 0):
                    active_gyo += 1

            for yoko in range(table1.shape[2]):
                if active_gyo != 0:
                    table1[mai, 5, yoko] = np.mean(table1[mai, :active_gyo, yoko])
                    table1[mai, 6, yoko] = np.std(table1[mai, :active_gyo, yoko], ddof=1)
                else:
                    table1[mai, 5, yoko] = 0
                    table1[mai, 6, yoko] = 0
        table1 = np.round(table1, 2)

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
            [len(table1[:, 0, 0]), len(table1[0, :, 0]), len(table1[0, 0, :])])
        colour_table = colour_table.astype(str)
        for k in range(len(table1[:, 0, 0])):
            for i in range(len(table1[0, :, 0])):
                for j in range(len(table1[0, 0, :])):
                    if i >= 5:
                        colour_table[k, i, j] = "lightgrey"
                    elif j >= 4:
                        colour_table[k, i, j] = "white"
                    elif table1[k, i, j] > 5.0:
                        if table1[k,i,j] > 10.0:
                            colour_table[k,i,j] = "red"
                        else:
                            colour_table[k, i, j] = "yellow"
                    else:
                        colour_table[k, i, j] = "white"
        return colour_table

    def plot_table_dyn(self):
        """
        plot table dynamically.
        :return:
        """
        if np.array(self.y1s).shape[1] > 2:
            plot_row = 2
        else:
            plot_row = 1

        table_fig = plt.figure(figsize=(8, plot_row * 2), dpi=100)
        # for i in range(self.table1.shape[0]):
        for i in range(np.array(self.y1s).shape[1]):
            v = i + 1
            # col_labels = ["x_L", "SI_L", "x_R", "SI_R", "3D"]
            col_labels = ["LR", "SI_V", "AP", "SI_H", "3D"]
            row_labels = ["1", "2", "3", "4", "5", "mean", "std"]
            ax1 = table_fig.add_subplot(plot_row, 2, v)
            t1 = ax1.table(cellText=self.table1[i, :, :], colLabels=col_labels,
                           rowLabels=row_labels, loc="center", cellColours=self.colour_table[i])
            t1.auto_set_font_size(False)
            t1.set_fontsize(10)
            t1.scale(1, 1)
            ax1.set_title("Marker motion #" + str(v) + " (mm)", color=wave_colors[i])
            ax1.set_axis_off()
            cell_height = 1 / 8.0
            for pos, cell in t1.get_celld().items():
                cell.set_height(cell_height)

            plt.tight_layout()
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
            for pos in ['right', 'top', 'bottom', 'left']:
                plt.gca().spines[pos].set_visible(False)
        table_fig.savefig(self.table_file)

    def create_fig3_dyn(self):
        fig3 = plt.figure(figsize=(8, len(self.time1s) * 2.5))
        new_times = self.new_times

        figure_number = len(self.time1s)
        for i in range(figure_number):
            ax1_3 = fig3.add_subplot(figure_number, 2, (i + 1) * 2 - 1)
            ax1_3.plot(self.time1s[i], self.x1s[i], "o", label="resp.", c="r", alpha=0.5)
            for j in range(np.array(self.y1s).shape[1]):
                ax1_3.plot(self.time1s[i], self.y1s[i][j], "-.", label="#" + str(j + 1), c=wave_colors[j])
            ax1_3.set_title("Acq. time: " + str(new_times[i][0]))
            ax1_3.set_ylabel("Resp. Phase")
            ax1_3.set_xlabel("Time (sec.)")
            if i == 0:
                ax1_3.legend()

            axy1_3 = fig3.add_subplot(figure_number, 2, (i + 1) * 2)
            axy1_3.plot(self.time2s[i], self.x2s[i], "o", label="resp.", c="r", alpha=0.5)
            for j in range(np.array(self.y2s).shape[1]):
                axy1_3.plot(self.time2s[i], self.y2s[i][j], "-.", label="#" + str(j + 1), c=wave_colors[j])
            axy1_3.set_title("Acq. time: " + new_times[i][1])
            axy1_3.set_ylabel("Resp. Phase")
            axy1_3.set_xlabel("Time (sec.)")
            # axy1_3.legend()

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
        newx1s, newy1s = [], []
        for j in range(y1s.shape[1]):  # マーカー数
            x1_for_a_marker = []
            y1_for_a_marker = []
            for i in range(x1s.shape[0]):  # 測定繰り返し数
                if np.corrcoef(x1s[i, :], y1s[i, j, :])[0, 1] > threshold:
                    x1_for_a_marker = np.concatenate([x1_for_a_marker, x1s[i, :]])
                    y1_for_a_marker = np.concatenate([y1_for_a_marker, y1s[i, j, :]])
                if np.corrcoef(x2s[i, :], y2s[i, j, :])[0, 1] > threshold:
                    x1_for_a_marker = np.concatenate([x1_for_a_marker, x2s[i, :]])
                    y1_for_a_marker = np.concatenate([y1_for_a_marker, y2s[i, j, :]])
            newx1s.append(x1_for_a_marker)
            newy1s.append(y1_for_a_marker)

        x_1_2_flat = newx1s
        y_1_2_2d = newy1s
        if y1s.shape[1] > 2:
            plot_row = 2
        else:
            plot_row = 1

        fig4 = plt.figure(figsize=(8, plot_row * 4))
        for i in range(len(y_1_2_2d)):  # iはマーカーの数
            v = i + 1
            ax1 = fig4.add_subplot(plot_row, 2, v)
            ax1.plot(x_1_2_flat[i], y_1_2_2d[i], "o", c=wave_colors[i], mfc="None", alpha=0.5)
            ax1.set_xlabel("Resp. phase (%)")
            ax1.set_ylabel("Marker phase (%)")
            ax1.set_title("Marker #" + str(i + 1), c=wave_colors[i])
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

    def create_html(self):
        """
        Create html report
        :return:
        """
        with open(self.report_file) as inf:
            txt = inf.read()
            soup = BeautifulSoup(txt, features="lxml")

        tag_pid = soup.new_tag("p")
        tag_pid.string = "Study Date: " + str(self.study_date[0]) + ", Report created: " + datetime.datetime.now().strftime("%Y%m%d %H:%M:%S") +"\n"
        soup.body.append(tag_pid)

        tag_pid = soup.new_tag("p")
        tag_pid.string = "Patient ID: " + os.path.split(self.folder)[1] + "\n"
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
