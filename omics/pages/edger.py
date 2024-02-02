from PyQt5.QtWidgets import QWidget, QTableWidgetItem, QLabel
from PyQt5.QtWidgets import *
from ..gui.youtube_ui import Ui_Form
from ..scripts import process
from PyQt5 import QtCore
import os
import pickle
import pathlib
#from ..rnalysis import filtering
import omics.main_module as filtering
#import qpageview
from PyQt5 import QtWebEngineWidgets
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QPixmap

class EdgeR(QWidget):
    submitClicked3 = QtCore.pyqtSignal(list)
    def __init__(self, dat):
        super(EdgeR, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.input_data = dat # a list contain
        self.filtered_file  = ''

        self.start = self.ui.pushButton
        self.start.clicked.connect(self.start_un)

        #self.count_table = self.ui.tableWidget_3
        #self.meta_table = self.ui.tableWidget_2
        self.cond1_combo = self.ui.comboBox_2
        self.cond2_combo = self.ui.comboBox_3

        dfram = self.input_data[3]
        items = dfram['condition']
        self.cond1_combo.addItems(list(set(items)))
        self.cond2_combo.addItems(list(set(items)))


        self.work_flow_combo = self.ui.comboBox
        self.work_flow_combo.addItems(["DESeq2", "EdgeR", "EMOGEA"])
        self.submitClicked3.connect(self.procces_filtering)

        
        ##self.view = QtWebEngineWidgets.QWebEngineView()
        #self.v = qpageview.View()
        #self.view.loadPdf("/Users/ibrahimahmed/projects/GUI2/heatmap_LSK.pdf") 
        ##self.view.load(QtCore.QUrl().fromLocalFile(
            ##"/Users/ibrahimahmed/projects/GUI2/heatmap_LSK.pdf")) 
        
    
        self.label1 = self.ui.label_5
        self.label2 = self.ui.label_6
        self.label3 = self.ui.label_7
        #self.stack1.addWidget(label)
        



        #self.data_view()
        #self.view_meta()
        self.check_recived()

    def check_recived(self):
        
        dframe = self.input_data[0]
        design = self.input_data[1]
        matrix = self.input_data[2]
        
        #df = dframe.to_csv(header=None, index=False)
        print(")))))))))))))))))")
        #print(dframe['samples'])
        print(self.input_data)
        

    def start_un(self):
        text = self.work_flow_combo.currentText()
        index = self.work_flow_combo.currentIndex()
        if text == "EdgeR":
            
            count_data_path = pathlib.Path(__file__).parent.parent.joinpath(self.input_data[0])
            design_matrix_path = str(self.input_data[-1]) #pathlib.Path(__file__).parent.parent.joinpath(self.input_data[0])
           
            filtDesq = filtering.CountFilter(count_data_path)
            filtDesq.differential_expression_edger3(design_matrix_path, "untreated", r_installation_folder="/Library/Frameworks/R.framework/Resources",output_folder="/Users/ibrahimahmed/projects/GUI/result_dir/edger_result")
        
            msg = QMessageBox()
            msg.setText("complete!/n The result directory is: /Users/ibrahimahmed/projects/GUI/result_dir")
            x = msg.exec_()

        if text == "DESeq2":
            count_data_path = pathlib.Path(__file__).parent.parent.joinpath(self.input_data[0])
            design_matrix_path = str(self.input_data[1]) #pathlib.Path(__file__).parent.parent.joinpath(self.input_data[0])
            compare = (("condition", self.cond1_combo.currentText(), self.cond2_combo.currentText()),)
            filtDesq = filtering.CountFilter(count_data_path)
            filtDesq.differential_expression_deseq2(design_matrix_path, compare, r_installation_folder="/Library/Frameworks/R.framework/Resources",output_folder="/Users/ibrahimahmed/projects/GUI/result_dir/DESeq2_result")
            msg = QMessageBox()
            msg.setText("complete!/n The result directory is: /Users/ibrahimahmed/projects/GUI/result_dir")
            x = msg.exec_()

            pixmap1 = QPixmap('/Users/ibrahimahmed/projects/GUI2/heatmap_LSK.jpeg')
            self.label1.setPixmap(pixmap1)
            self.resize(pixmap1.width(), pixmap1.height())

            pixmap2 = QPixmap('/Users/ibrahimahmed/projects/GUI2/clusterProfiler_GSEA_IVANOVA_small.svg')
            self.label2.setPixmap(pixmap2)
       
            pixmap3 = QPixmap("/Users/ibrahimahmed/projects/my_plot.svg")
            self.label3.setPixmap(pixmap3)

        if index == 2:
            print("+++++++++++++++++++++++++++++++")
            print(self.input_data)
            print("+++++++++++++++++++++++++++++++")



    def data_view(self):
        pass
    
    def view_meta(self):
        pass
    
    def procces_filtering(self, fil):
        self.filtered_file = fil



if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    w = EdgeR()
    w.show()
    sys.exit(app.exec_())

