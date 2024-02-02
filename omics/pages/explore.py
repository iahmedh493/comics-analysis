from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QWidget, QFileDialog
from PyQt5.QtWidgets import *
from ..gui.explore_ui import Ui_Form
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5 import QtCore
import pickle
from PyQt5 import QtCore, QtWidgets #, QtWebEngineWidgets
from PyQt5.QtWidgets import (QWidget, QPushButton,
                             QHBoxLayout, QVBoxLayout, QApplication, QLabel, QSpinBox, QComboBox)

from ..scripts.model import ExplorAnalysis
import scanpy as sc
# from pages.toyota import Toyota

import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from PyQt5 import QtWidgets as qtwo
from PyQt5 import QtWebEngineWidgets

#from ..rnalysis import filtering
import omics.main_module as filtering

# Inital setting for plot size
from matplotlib import rcParams
FIGSIZE=(3,3)
rcParams['figure.figsize']=FIGSIZE

import warnings
warnings.filterwarnings('ignore')

#from pages import ComboBoxWidget

import io
from contextlib import redirect_stdout
import sys

class Explore(QWidget):
    #signal1 = QtCore.pyqtSignal(list)
    submitClicked2 = QtCore.pyqtSignal(list)

    def __init__(self, list1):
        super(Explore, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.loaded = list1
        self.view = QtWebEngineWidgets.QWebEngineView()
        self.stack = self.ui.stackedWidget
        self.filtered_file = ''
        self.normalized_file = ''
        
       
        self.apply = self.ui.pushButton
        self.submit = self.ui.pushButton_4
    
        #self.run.clicked.connect(self.show_visualize)

        self.plots = self.ui.comboBox
        self.plots.setPlaceholderText("---Choose function---")
        self.plots.addItems(["PCA", "Volcano"])

        self.apply.clicked.connect(self.apply_plots)
        self.submit.clicked.connect(self.confirm_form)
        self.submitClicked2.connect(self.procces_filtering)

        
        #self.submit.clicked.connect(self.confirm_form)
    

    def confirm_form(self):
        #print(self.loaded)
        temp_list = [self.loaded[-2], self.loaded[-1], self.loaded[-3]]#self.normalized_file]
        self.submitClicked2.emit(self.loaded)
        self.close()

    def apply_plots(self):
        
        
        file_name = 'temp2.csv'
        fun = self.plots.currentText()
        index = self.plots.currentIndex()

        filtDesq = filtering.CountFilter(self.loaded[0])

        filtDesq.differential_expression_edger2("/Users/ibrahimahmed/projects/GUI/Test.csv",
                                        (("condition", "treated", "untreated"),),
                                          r_installation_folder="/Library/Frameworks/R.framework/Resources",
                                          output_folder="/Users/ibrahimahmed/projects/GUI/result_dir/exploration")

        if index == 0:
            self.view.load(QtCore.QUrl().fromLocalFile(
            "/Users/ibrahimahmed/projects/GUI/result_dir/exploration/MDS-Plot.html"))
            self.stack.addWidget(self.view)
            self.view.show()

        #file_name = 'temp2.csv'
        
        if index == 1:
            self.view.load(QtCore.QUrl().fromLocalFile(
            "/Users/ibrahimahmed/projects/GUI/result_dir/exploration/volicano-plot.html"))
            self.stack.addWidget(self.view)
            self.view.show()
        
        
    def procces_filtering(self, fil):
        self.filtered_file = fil

    

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    w = Explore()
    w.show()
    sys.exit(app.exec_())


















