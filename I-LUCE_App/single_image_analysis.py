import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QAction, QTextEdit, 
    QFileDialog, QMenuBar, QWidget
)
from PyQt5.QtCore import Qt

class SingleImageAnalysisApp(QApplication):
    def __init__(self):
        super().__init__(sys.argv)
        self.window = SingleImageAnalysisWindow()
        self.window.show()
        sys.exit(self.exec_())

class SingleImageAnalysisWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Single Image Analysis")
        self.setGeometry(600, 100, 1000, 900)
        self.setStyleSheet("background-color: white;")

        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')

        open_action = QAction('Apri', self)
        save_action = QAction('Salva', self)
        exit_action = QAction('Esci', self)

        file_menu.addAction(open_action)
        file_menu.addAction(save_action)
        file_menu.addAction(exit_action)

        open_action.triggered.connect(self.openFile)
        save_action.triggered.connect(self.saveFile)
        exit_action.triggered.connect(self.close)

        layout.setMenuBar(menubar)

        self.text_edit = QTextEdit()
        layout.addWidget(self.text_edit)

        central_widget = QWidget(self)
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def openFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly

        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self, "Open a File", "", "Text Files (*.txt);;All Files (*)", options=options
        )

        if file_path:
            try:
                with open(file_path, 'r') as file:
                    file_content = file.read()
                self.text_edit.setPlainText(file_content)
            except Exception as e:
                print(str(e))

    def saveFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly

        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getSaveFileName(
            self, "Save a File", "", "Text Files (*.txt);;All Files (*)", options=options
        )

        if file_path:
            try:
                with open(file_path, 'w') as file:
                    file_content = self.text_edit.toPlainText()
                    file.write(file_content)
            except Exception as e:
                print(str(e))

if __name__ == '__main__':
    app = SingleImageAnalysisApp()
