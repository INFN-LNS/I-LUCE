import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenuBar, QAction, QWidget
from rcf_analysis import RCFAnalysisWindow
from TOF_analysis import TOFAnalysisWindow



class MyApplication(QApplication):
    def __init__(self):
        super().__init__(sys.argv)
        self.main_window = MyMainWindow()
        self.main_window.show()
        sys.exit(self.exec_())

class MyMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Imposta il titolo della finestra principale
        self.setWindowTitle("App con Barra di Menu")

        # Imposta le dimensioni della finestra principale
        self.setGeometry(100, 100, 800, 600)

        # Crea la barra dei menu
        menubar = self.menuBar()

        # Crea il menu "Analysis" nella barra dei menu
        file_menu = menubar.addMenu("Analysis")

        # Crea il comando "RCF Analysis" nel menu "Analysis"
        rcf_action = QAction("RCF Analysis", self)
        rcf_action.triggered.connect(self.open_rcf_analysis_window)
        file_menu.addAction(rcf_action)
        
        # Crea il comando "TOF Analysis" nel menu "Analysis"
        TOF_action = QAction("TOF Analysis", self)
        TOF_action.triggered.connect(self.open_TOF_analysis_window)
        file_menu.addAction(TOF_action)
        
        # Crea il comando "Chiudi" nel menu "File"
        close_action = QAction("Chiudi", self)
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)

    def open_rcf_analysis_window(self):
        self.rcf_window = RCFAnalysisWindow()  # Usa la classe dal file separato
        self.rcf_window.show()
        
    def open_TOF_analysis_window(self):
        self.TOF_window = TOFAnalysisWindow()  # Usa la classe dal file separato
        self.TOF_window.show()

if __name__ == '__main__':
    app = MyApplication()
