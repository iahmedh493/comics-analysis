
import sys

from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QPixmap

from MainWindow_ui import Ui_MainWindow
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from PyQt5 import uic
from omics.pages.home import Home
from omics.pages.load import LoadData
from omics.pages.explore import Explore
from omics.pages.edger import EdgeR

from omics.pages.enrich import Enrich
from omics.pages.single import SingleCell



class MyMainWindow(QMainWindow):

    def __init__(self):
        super(MyMainWindow, self).__init__()

        # uic.loadUi('ui/home.ui', self)
        # uic.loadUi('toyota.ui', self)
        self.main_ui = Ui_MainWindow()
        self.main_ui.setupUi(self)
        # ******************************************
        # initialize elements
        self.data_loaded = []
        self.explored_dat = []
        
       
        # **********************************************
        # create dict for the function buttons
        self.home_btn = self.main_ui.pushButton
        self.load_data_btn = self.main_ui.pushButton_4
        self.explore_btn = self.main_ui.pushButton_5
        self.edger_btn = self.main_ui.pushButton_7
        self.enrichment_btn = self.main_ui.pushButton_11
        self.single_cell_btn = self.main_ui.pushButton_12

        # **************************************************
        # initialize home page button when start app
        self.home_btn.setDisabled(True)
        self.home_btn.setStyleSheet("background-color: rgb(222,222,222);")
        self.home_btn.setStyleSheet("color: rgb(0,0,0);")
        # *************************************************
        # connect buttons with function
        self.home_btn.clicked.connect(self.show_home)
        self.load_data_btn.clicked.connect(self.show_load_data)
        self.explore_btn.clicked.connect(self.show_explore)
        self.edger_btn.clicked.connect(self.show_edger)
        self.enrichment_btn.clicked.connect(self.show_enrichment)
        self.single_cell_btn.clicked.connect(self.show_single)

       

    def show_home(self):
        self.home = Home()
        self.main_ui.stackedWidget.addWidget(self.home)
        self.main_ui.stackedWidget.setCurrentIndex(1)
        self.home.show()
    
    def show_load_data(self):
        self.load = LoadData()
        self.main_ui.stackedWidget.addWidget(self.load)
        self.main_ui.stackedWidget.setCurrentIndex(2)
        self.load.submitClicked.connect(self.process_loaded)
        self.load.show()
    
    def process_loaded(self, listd):

        self.data_loaded = listd
        print("*****************************")
        print(self.data_loaded)
        print("*****************************")

    def show_explore(self):
        self.explore = Explore(self.data_loaded)
        self.main_ui.stackedWidget.addWidget(self.explore)
        self.explore.submitClicked2.connect(self.process_explore)
        self.main_ui.stackedWidget.setCurrentIndex(3)
        self.explore.show()

    def process_explore(self, list2):
        self.explored_dat = list2


    def show_edger(self):
        self.edger = EdgeR(self.explored_dat)
        self.main_ui.stackedWidget.addWidget(self.edger)
        self.main_ui.stackedWidget.setCurrentIndex(4)
        self.edger.show()

    def show_enrichment(self):
        self.enrich = Enrich()
        self.main_ui.stackedWidget.addWidget(self.enrich)
        self.main_ui.stackedWidget.setCurrentIndex(6)
        self.enrich.show()
    
    def show_single(self):
        self.single1 = SingleCell()
        self.main_ui.stackedWidget.addWidget(self.single1)
        self.main_ui.stackedWidget.setCurrentIndex(7)
        self.single1.show()




if __name__ == "__main__":

    app = QApplication(sys.argv)

    window = MyMainWindow()
    window.show()

    sys.exit(app.exec_())






