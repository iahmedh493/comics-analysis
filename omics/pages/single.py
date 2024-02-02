
from PyQt5.QtWidgets import QWidget, QFileDialog
from PyQt5.QtWidgets import *
from ..gui.singleCell_ui import Ui_Form



class SingleCell(QWidget):
    
    def __init__(self):
        super(SingleCell, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        
    

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    w = SingleCell()
    w.show()
    sys.exit(app.exec_())


















