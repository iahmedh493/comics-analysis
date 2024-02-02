
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class CheckableComboBox(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.combo_box = CheckableComboBoxChild(self)

        self.delete_btn = QPushButton("Delete Selected")
        self.delete_btn.clicked.connect(self.combo_box.deleteSelected)

        layout = QVBoxLayout(self)
        layout.addWidget(self.combo_box)
        layout.addWidget(self.delete_btn)

class CheckableComboBoxChild(QComboBox):
    # Subclass Delegate to increase item height
    class Delegate(QStyledItemDelegate):
        def sizeHint(self, option, index):
            size = super().sizeHint(option, index)
            size.setHeight(20)
            return size

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Make the combo editable to set a custom text
        self.setEditable(True)
        self.setInsertPolicy(QComboBox.NoInsert)
        self.lineEdit().setPlaceholderText("--Select Actions--")
        edit = self.lineEdit()
        self.setLineEdit(edit)
        #self.completer = QCompleter()
        #self.completer.setFilterMode(Qt.MatchContains)
        #self.completer.setCaseSensitivity(Qt.CaseInsensitive)
        #edit.setCompleter(self.completer)
        #self.completer.setModel(self.model())
        edit.returnPressed.connect(self.insertCustomItem)

        # Make the lineedit the same color as QPushButton
        palette = qApp.palette()
        palette.setBrush(QPalette.Base, palette.button())
        self.lineEdit().setPalette(palette)

        # Use custom delegate
        self.setItemDelegate(CheckableComboBoxChild.Delegate())

        # Update the text when an item is toggled
        self.model().dataChanged.connect(self.updateText)

        # Hide and show popup when clicking the line edit
        self.lineEdit().installEventFilter(self)

        # Prevent popup from closing when clicking on an item
        self.view().viewport().installEventFilter(self)

    def resizeEvent(self, event):
        # Recompute text to elide as needed
        self.updateText()
        super().resizeEvent(event)

    def eventFilter(self, object, event):
        if object == self.view().viewport():
            if event.type() == QEvent.MouseButtonRelease:
                if self.lineEdit().hasFocus():
                    return True
                index = self.view().indexAt(event.pos())
                item = self.model().item(index.row())
                if item.checkState() == Qt.Checked:
                    item.setCheckState(Qt.Unchecked)
                else:
                    item.setCheckState(Qt.Checked)
                return False
        return False

    def timerEvent(self, event):
        # After timeout, kill timer, and reenable click on line edit
        self.killTimer(event.timerId())

    def updateText(self):
        texts = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                texts.append(self.model().item(i).text())
        text = ", ".join(texts)

        # Compute elided text (with "...")
        metrics = QFontMetrics(self.lineEdit().font())
        elidedText = metrics.elidedText(text, Qt.ElideRight, self.lineEdit().width())
        self.lineEdit().setText(elidedText)

    def addItem(self, text, data=None, checked=False):
        item = QStandardItem()
        item.setText(text)
        if data is None:
            item.setData(text)
        else:
            item.setData(data)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Unchecked if not checked else Qt.Checked, Qt.CheckStateRole)
        self.model().appendRow(item)

    def addItems(self, texts, datalist=None):
        for i, text in enumerate(texts):
            try:
                data = datalist[i]
            except (TypeError, IndexError):
                data = None
            self.addItem(text, data)

    def insertCustomItem(self):
        text = self.lineEdit().text().strip()
        # Split the text by comma, and lowercase and strip() each piece
        typedItemsOriginal = [item.strip() for item in text.split(",") if item.strip()]
        typedItemsLower = [item.lower() for item in typedItemsOriginal]
        # Uncheck all items
        for i in range(self.model().rowCount()):
            self.model().item(i).setData(Qt.Unchecked, Qt.CheckStateRole)
        # Loop through each item in the text and check it, if it exists in
        # lowercase
        for i in range(len(typedItemsOriginal)):
            for j in range(self.model().rowCount()):
                if self.model().item(j).text().lower() == typedItemsLower[i]:
                    self.model().item(j).setData(Qt.Checked, Qt.CheckStateRole)
                    break
            else:
                # If the item doesn't exist, add it to the list
                self.addItem(typedItemsOriginal[i], checked=True)
        self.updateText()
        self.showPopup()

    def currentData(self):
        # Return the list of selected items data
        res = []
        for i in range(self.model().rowCount()):
            if self.model().item(i).checkState() == Qt.Checked:
                res.append(self.model().item(i).data())
        return res

    def deleteSelected(self):
        text = self.lineEdit().text().strip()
        # Split the text by comma, and lowercase and strip() each piece
        checked_items = [item.strip() for item in text.split(",") if item.strip()]
        if not checked_items:
            QMessageBox.information(self, "No Selection", "No options were selected, to delete!")
            return
        confirmDeletion = QMessageBox.question(self, 'Confirm Deletion', f"Are you sure you want to delete {', '.join(checked_items)}?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        anyRemoved = False
        if confirmDeletion == QMessageBox.Yes:
            for item in checked_items:
                index = self.findText(item)
                if index != -1:
                    self.removeItem(index)
                    anyRemoved = True
        if anyRemoved:
          self.updateText()
          self.showPopup()