from PyQt5.QtWidgets import QWidget

from ..gui.enrich_ui import Ui_Form
from ..scripts import process
import os
import pickle
from PyQt5.QtWidgets import *
from PyQt5.QtCore import QAbstractTableModel, Qt

from PyQt5 import QtCore, QtGui
#from scripts.test2 import PandasModel
from .explore import Explore
import pandas as pd
import pickle
import pathlib
#from rnalysis import filtering
import omics.main_module as filtering
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QPixmap

from omics.utils import io

class PandasModel(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, parnet=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole or role == Qt.EditRole:
                value = self._data.iloc[index.row(), index.column()]
                return str(value)

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            self._data.iloc[index.row(), index.column()] = value
            return True
        return False

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

class Enrich(QWidget):

    def __init__(self):
        super(Enrich, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.file_data = ''

        self.start = self.ui.pushButton
        self.start.clicked.connect(self.start_un)
        self.label1 = self.ui.label_2
        self.label2 = self.ui.label_3

        self.load = self.ui.pushButton_2
        self.load.clicked.connect(self.load_file)

        #self.count_table = self.ui.tableView
    
    def load_file(self):
        f_mata = QFileDialog.getOpenFileNames(self, "Open File", "", "All Files (*);;")
        self.file_data = f_mata[0][0]
        print(self.file_data)
        #self.count_view()
 

    def start_un(self):
        filtDesq = filtering.CountFilter(self.file_data)
        
        filtDesq.differential_expression_pathway(r_installation_folder="/Library/Frameworks/R.framework/Resources",output_folder="/Users/ibrahimahmed/projects/GUI/result_dir/enrichment")
        msg = QMessageBox()
        #msg.setWindowTitle("Filter by row sum")
        msg.setText("complete!")
        x = msg.exec_()

        msg.setText("The result directory is: /Users/ibrahimahmed/projects/GUI/result_dir")
        x = msg.exec_()

        pixmap3 = QPixmap("/Users/ibrahimahmed/projects/GUI2/barplot11.png")
        self.label2.setPixmap(pixmap3)

        pixmap2 = QPixmap("/Users/ibrahimahmed/projects/GUI2/barplot22.png")
        self.label1.setPixmap(pixmap2)

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    w = Enrich()
    w.show()
    sys.exit(app.exec_())

