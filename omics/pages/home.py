from PyQt5.QtWidgets import *
from ..gui.hom_ui import Ui_Form

# from test2 import PandasModel

import pandas as pd


class Home(QWidget):

    def __init__(self):
        super(Home, self).__init__()

        self.ui = Ui_Form()
        self.ui.setupUi(self)