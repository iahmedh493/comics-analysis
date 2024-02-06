# import pandas as pd
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets
from PyQt5.QtCore import QAbstractTableModel, Qt
from ..gui.toyota_ui import Ui_Form
from PyQt5 import QtCore, QtGui
#from scripts.test2 import PandasModel
from .explore import Explore
import pandas as pd
import pickle
import pathlib
from pathlib import Path
from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtWidgets import QDialog, QComboBox, QApplication, QHeaderView, QPushButton, QTableView

from omics.utils import io
#from ..rnalysis import filtering
import omics.main_module as filtering
from typing import Any, Dict, Union, List, Tuple
from PyQt5.QtCore import QProcess
import os


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
        #return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable
        return QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsSelectable
    
class LoadData(QWidget):
    #signal = QtCore.pyqtSignal(list)
    submitClicked = QtCore.pyqtSignal(list)

    def __init__(self):
        super(LoadData, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.pass_data = [] # hold path of final dataset

        # Input widget_6
        self.result_dir_btn = self.ui.pushButton
        self.load_matrix_btn = self.ui.pushButton_2
        self.load_meta_btn = self.ui.pushButton_3
        self.loaded = self.ui.pushButton_5

        # Summary and validation: groupBox_2
        self.dim_count = self.ui.label_10
        self.dim_count_rowas = self.ui.label_11
        self.count_mata_compar = self.ui.label_12

        # Filtering: Widget
        self.filter_threeshold = self.ui.doubleSpinBox
        self.apply_filter = self.ui.pushButton_4

        self.filter_threeshold.setValue(20.0)
        
        # Normalization: widget_2
        self.yes_normalized = self.ui.radioButton
        self.not_normalized = self.ui.radioButton_2
        self.normalize_label = self.ui.label_8
        self.apply_normalize = self.ui.pushButton_7

        #*******************************************
        # Summary after filtering
        self.summary_after_filter = self.ui.label_12
        #********************************************

        # Summary after filter
        self.group_box2 = self.ui.groupBox_3
        self.dim_count_afte = self.ui.label_12
        self.dim_count_rowas_after = self.ui.label_13
        self.group_box2.setEnabled(False)
        #********************************************
        # Meta edit
        self.meta_table = self.ui.tableView
        self.save_meta_btn = self.ui.pushButton_6
        self.submit = self.ui.pushButton_8
        #********************************************
        self.edit_meta_data = self.ui.pushButton_9
        self.edit_meta_data.clicked.connect(self.edit_meta)
        #********************************************
        # Define variables
        self.count_path = ""
        self.meta_path = ""
        self.result_path = ""
        #*********************************************

        # data
        self.count_df = None
        self.meta_df = None


        self.group_box1 = self.ui.groupBox_2
        #self.group_box1.hide()
        self.group_box1.setEnabled(False)
        ###print(self.group_box1.isEnabled())
        
        # ******************************************************
        # Connect signal and slot
        # *****************************************************
        
        self.load_matrix_btn.clicked.connect(self.load_count_matrix)
        self.load_meta_btn.clicked.connect(self.load_meta_data)
        self.result_dir_btn.clicked.connect(self.get_result_dir)
        self.loaded.clicked.connect(self.validate_input)

        self.apply_filter.clicked.connect(self.filter_data)
        self.apply_normalize.clicked.connect(self.normalize_data)

    
        self.save_meta_btn.clicked.connect(self.IP_File_Import)
        self.submit.clicked.connect(self.accept_data)

        #self.apply_filter.setEnabled(False)
        #self.apply_normalize.setEnabled(False)

        self.filter_widget = self.ui.widget
        self.filter_widget.setEnabled(False)

        self.normalize_widget = self.ui.widget_2
        self.normalize_widget.setEnabled(False)

        self.summary_widget = self.ui.groupBox_3
        self.summary_widget.setEnabled(False)

        self.organism = self.ui.comboBox
        self.gene_id = self.ui.comboBox_2

        self.model_organism = ''
        self.gene_ids = ''

        self.organism.currentIndexChanged.connect(self.get_organism)
        self.gene_id.currentIndexChanged.connect(self.get_gene_id)
        self.frame = self.ui.frame_2
        self.frame.setEnabled(False)
    
    def edit_meta(self):
        qm = QMessageBox()
        ret = qm.question(self,'', "Are you sure to edit meta data?", qm.Yes | qm.No)
        if ret == qm.Yes:
            self.frame.setEnabled(True)
        else:
            self.filter_widget.setEnabled(True)
            self.normalize_widget.setEnabled(True)
            self.group_box1.setEnabled(True)


    def get_organism(self):
        self.model_organism = self.organism.currentText()
        ##self.pass_data.append(self.model_organism)
        #print(self.model_organism)
        

    def get_gene_id(self):
        self.gene_ids = self.gene_id.currentText()
        ##self.pass_data.append(self.gene_ids)
        #print(self.gene_ids)



    def load_count_matrix(self):
        self.count_path = QFileDialog.getOpenFileNames(self, "Open File", "", "All Files (*);;")[0][0]
        print(self.count_path)
        #if not self.count_path[0] == "":
        self.pass_data.append(self.count_path)

    def load_meta_data(self):
        self.meta_path = QFileDialog.getOpenFileNames(self, "Open File", "", "All Files (*);;")[0][0]
        ###print(self.meta_path)
        #if not self.count_path[0] == "":
        self.pass_data.append(self.meta_path)
    
    def get_result_dir(self):
        self.result_path = QFileDialog.getExistingDirectory(self, 'Select Folder')
        self.pass_data.append(self.result_path)
        ###print(self.pass_data)

    def validate_input(self):
        ###print(self.count_path)
        if self.count_path == "":
            self.show_critical_messagebox("Please load count matrix")
        elif self.meta_path == "":
            self.show_critical_messagebox("Please load meta data")
        
        elif self.result_path == "":
            self.show_critical_messagebox("Please select result directory")

        elif self.model_organism == '':
            self.show_critical_messagebox("Please selected model organism")

        elif self.gene_ids == '':
            self.show_critical_messagebox("Please selected id type")
        
        else:
            self.count_df = self.reat_data(self.count_path)
            self.meta_df = self.reat_data(self.meta_path)
            if not self.meta_df.shape[0] == self.count_df.shape[1]:
                self.show_critical_messagebox("The number of items in the design matrix %s does not match the number of \
                                              columns in the count matrix %s." %(self.meta_df.shape[0], self.count_df.shape[1]))
                
            if not sorted(self.meta_df.index) == sorted(self.count_df.columns):
                self.show_warning_messagebox("The sample names in the design matrix do not match the sample names \
                                              in the count matrix:%s != " %(self.meta_df.index, sorted(self.count_df.columns)))

            num_col = "%s columns" % self.count_df.shape[1]
            num_rows = "%s rows" % self.count_df.shape[0]
            messg = num_col + num_rows
            self.dim_count.setText(num_col)
            self.dim_count_rowas.setText(num_rows)
            self.group_box1.setEnabled(True)
            self.load_table2()

            print(self.pass_data)



    def show_critical_messagebox(self, message): 
        msg = QMessageBox() 
        msg.setIcon(QMessageBox.Critical) 
        # setting message for Message Box 
        msg.setText(message) 
        # setting Message box window title 
        msg.setWindowTitle("Critical MessageBox") 
        # declaring buttons on Message Box 
        msg.setStandardButtons(QMessageBox.Ok) 
        # start the app 
        retval = msg.exec_()
    
    def show_warning_messagebox(self, message): 
        msg = QMessageBox() 
        msg.setIcon(QMessageBox.Warning) 
        # setting message for Message Box 
        msg.setText(message) 
        # setting Message box window title 
        msg.setWindowTitle("Warning MessageBox") 
        # declaring buttons on Message Box 
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel) 
        # start the app 
        retval = msg.exec_()

    def reat_data(self, fname: Union[str, Path], drop_columns: Union[str, List[str]] = None):
           
        """
        Load a table.

        :param fname: full path/filename of the .csv file to be loaded into the Filter object
        :type fname: Union[str, Path]
        :param drop_columns: if a string or list of strings are specified, \
        the columns of the same name/s will be dropped from the loaded table.
        :type drop_columns: str, list of str, or None (default=None)

        :Examples:
            >>> from rnalysis import filtering
            >>> d = filtering.Filter("tests/test_files/counted.csv")

        """
        # init from a file (load the csv into a DataFrame/Series)
        assert isinstance(fname, (str, Path))
        self.fname = Path(fname)
        self.df = io.load_table(fname, 0, squeeze=True, drop_columns=drop_columns)
        if isinstance(self.df, pd.Series):
            self.df = self.df.to_frame()
        # check for duplicate indices
        if self.df.index.has_duplicates:
            self.show_warning_messagebox("This loaded data contains multiple rows with the same name/index.")
        
        return self.df
    
    def filter_data(self):
        if self.group_box1.isEnabled():
            #>>> from rnalysis import filtering
            #>>> c = filtering.CountFilter('tests/test_files/counted.csv')
            #>>> c.filter_by_row_sum(5) # remove all rows whose sum is <5
            #Filtered 4 features, leaving 18 of the original 22 features. Filtered inplace.  
            #filter = filtering.CountFilter.from_dataframe(self.count_df, "/data")
            self.filter = filtering.CountFilter(self.count_path, self.result_path)
            self.filter.filter_by_row_sum(5) # remove all rows whose sum is <5
            self.filtered_df = self.filter.df
            threeshold = self.filter_threeshold.value()

            num_col = "%s columns" % self.filtered_df.shape[1]
            num_rows = "%s rows" % self.filtered_df.shape[0]
            
            self.dim_count_afte.setText(num_col)
            self.dim_count_rowas_after.setText(num_rows)
            
        else:
            self.show_warning_messagebox("Make sure all data loaded")
    
    def normalize_data(self):
        print("*********************")
        print("*********************")
        print(self.filter)
        print("*********************")
        print("*********************")
        #self.filter.normalize_tmm()
        ##self.group_box2.setEnabled(True)
        #self.load_table2()
        #print(self.filter.df)

    def load_table2(self):
        '''
        df = self.reat_data(self.meta_path)
        print(df)
        df['sample'] = df.index
        print(df)
        self.model = PandasModel(df)
        self.meta_table.setModel(self.model)
        '''
        file_data = open(self.meta_path, "r").readlines()
        file_list = [item.strip().split(',') for item in file_data]
        
        colnam = file_list[0]
        self.df = pd.DataFrame(file_list, columns=colnam)
        self.model = PandasModel(self.df)
        self.meta_table.setModel(self.model)
       
            

    def accept_data(self): 
       
        self.submitClicked.emit(self.pass_data)

        ###if len(self.pass_data) > 1:
            ###print(self.pass_data)
        ###else:
            ###print("Single files")
        '''
        indexes = self.meta_table.selectionModel().selectedRows()
        for index in reversed(sorted(indexes)):
            #print('Row %d is selected' % index.row())
            print(self.meta_table.model().data(index))
            ##self.meta_table.removeRow(self.meta_table.model().data(index))
            #index = self.currentIndex()
            deleteconfirmation = QtWidgets.QMessageBox.critical(self.parent(), "Delete row", "Really delete the selected sample", QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
            if deleteconfirmation == QtWidgets.QMessageBox.Yes:
                
                self.model.removeRow(index.row())
                print("&*&*&**&&*&**&*&*&*&*&&**&")
                print("&*&*&**&&*&**&*&*&*&*&&**&")
                print(self.model._data)
                self.meta_table.setModel(self.model)
                self.meta_table.selectRow(index.row())
                #self.meta_table.removeRow(index.row())
                #self.model.submit()
                #logbox.warning("Payrate ID %i deleted" %self.model.index(index.row(),0).data()) #prints data from row being deleted and column 0 (ID)
                return
            else:
                return

        self.model.submitAll()
        self.model.select()
        print(self.model)
           
        
        self.submitClicked.emit(self.pass_data)
        if len(self.pass_data) > 1:
            print(self.pass_data)
        else:
            print("Single files")
        '''

    def get_files(self):
        if len(self.files_list) > 0:
            return self.files_list

    def IP_File_Import(self):
        
        newModel = self.meta_table.model()
        data = []
        for row in range(newModel.rowCount(0)):
            rowRes = []
            for column in range(newModel.columnCount()):
                index = newModel.index(row, column)
                item = newModel.data(index)
                if item != '':
                    rowRes.append(item)
            data.append(rowRes)
        dataFrame = pd.DataFrame(data[1:], columns=data[0])
        self.pass_data.append(dataFrame)
        #***************** Revisit******************
        dataFrame.set_index('samples', inplace=True)
        
        dataFrame.to_csv(pathlib.Path(__file__).parent.parent.joinpath("Test.csv"))#, index=False, header=False)
        #io.save_table(dataFrame, "Test.csv")
        self.pass_data.append(pathlib.Path(__file__).parent.parent.joinpath("Test.csv"))
        self.filter_widget.setEnabled(True)
        self.normalize_widget.setEnabled(True)
        self.group_box1.setEnabled(True)
        
        ###print(pathlib.Path(__file__).parent.parent.joinpath("Test.csv"))

 
    
if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    w = LoadData()
    w.show()
    sys.exit(app.exec_())

'''
        self.meta_file = self.f_mata[0][0]
        file_data = open(self.meta_file, "r").readlines()
        file_list = [item.strip().split(',') for item in file_data]
        print("************************")
        print(file_list)
        print("************************")
        colnam = file_list[0]
        self.df = pd.DataFrame(file_list, columns=colnam)
        self.model = PandasModel(self.df)
        self.meta_table.setModel(self.model)
        '''

   