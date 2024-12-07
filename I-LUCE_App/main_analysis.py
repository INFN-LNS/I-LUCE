
# Analisi
import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenuBar, QAction, 
                             QWidget, QVBoxLayout, QLabel, QHBoxLayout,
                             QSpacerItem, QSizePolicy)

from PyQt5.QtGui import QPixmap, QFont

from PyQt5.QtCore import Qt

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
        self.setGeometry(100, 100, 300, 100)
        
        # Imposta il colore di sfondo desiderato
        self.setStyleSheet("background-color: rgb(0, 50, 0);")  

        self.setWindowTitle('I-LUCE analisys platform')
        
        # Crea la bar®ra dei menu
        menubar = self.menuBar()

        # Crea il menu "Ana®lysis" nella barra principale dei menu
        #
        file_menu = menubar.addMenu("Analysis")

        # Crea il comando "RCF Analysis" nel menu "Analysis"
        #
        rcf_action = QAction("RCF Analysis", self)
        rcf_action.triggered.connect(self.open_rcf_analysis_window)
        file_menu.addAction(rcf_action)
        
        # Crea il comando "TOF Analysis" nel menu "Analysis"
        TOF_action = QAction("TOF Analysis", self)
        TOF_action.triggered.connect(self.open_TOF_analysis_window)
        file_menu.addAction(TOF_action)
        
        # Crea il menu "Calibration" nella barra principale dei menu
        #
        file_menu = menubar.addMenu("Calibration")
        
        # Crea il comando "RCF Dose Calibration" nel menu "Calibration"
        #
        rcf_dose_calibration_action = QAction("RCF Dose Calibration", self)
        #rcf_dose_calibration_action.triggered.connect(self.open_rcf_dose_calibration_window)
        file_menu.addAction(rcf_dose_calibration_action)
        
        
        self.displayWidgets()  # Aggiunto per visualizzare i widget

    def open_rcf_analysis_window(self):
        self.rcf_window = RCFAnalysisWindow()
        self.rcf_window.show()
        
    def open_TOF_analysis_window(self):
        self.TOF_window = TOFAnalysisWindow()
        self.TOF_window.show()
        
    def displayWidgets(self):
        """
        Set up widgets using QHBoxLayout and QVBoxLayout.
        """
        # Create label and button widgets
        
        title_1 = QLabel("I-LUCE analysis platform")
        title_1.setFont(QFont('Arial', 24))
        title_1.setAlignment(Qt.AlignCenter)    # Centra il testo
        title_1.setStyleSheet("color: white;")  # Imposta il colore del testo in bianco
        
        title_2 = QLabel("G.A.P. Cirrone")
        title_2.setFont(QFont('Arial', 24))
        title_2.setAlignment(Qt.AlignCenter)    # Centra il testo
        title_2.setStyleSheet("color: white;")  # Imposta il colore del testo in bianco
        
        title_3 = QLabel("INFN-LNS")
        title_3.setFont(QFont('Arial', 24))
        title_3.setAlignment(Qt.AlignCenter)    # Centra il testo
        title_3.setStyleSheet("color: white;")  # Imposta il colore del testo in bianco
        
        logo_label = QLabel()
        logo_pixmap = QPixmap("Logo_INFN.jpeg")  
        logo_pixmap = logo_pixmap.scaledToWidth(300)
        logo_label.setPixmap(logo_pixmap)
        logo_label.setAlignment(Qt.AlignCenter)
        
        
        # Create vertical layout and add widgets 
        v_box = QVBoxLayout()
        v_box.addWidget(title_1)
        v_box.addStretch(1)
        v_box.addWidget(title_2)
        v_box.addStretch(1)
        v_box.addWidget(title_3)
        
        # Create an horizontal layout
        h_box = QHBoxLayout()
        h_box.addLayout(v_box)
        h_box.addStretch(1)
        h_box.addWidget(logo_label)
      
        # Set main layout of the window
        central_widget = QWidget(self)
        central_widget.setLayout(h_box)
        self.setCentralWidget(central_widget)

if __name__ == '__main__':
    app = MyApplication()
