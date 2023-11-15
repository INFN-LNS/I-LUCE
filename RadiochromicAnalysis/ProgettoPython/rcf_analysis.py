#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 22:43:47 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenuBar, QAction, 
                             QWidget, QLabel, QWidgetAction, QFileDialog,
                             QTextEdit, QVBoxLayout, QDialog, QHBoxLayout)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt

from image_viewer import ImageViewer


class RCFAnalysisWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("RCF Analysis")
        self.setGeometry(800, 200, 1000, 800)
        
        # =============================  MENU BAR =============================
        #
        # Crea la barra dei menu
        menubar_RCF_analysis = self.menuBar()
        
        # Crea la voce "Open" nella barra dei menu
        file_menu = menubar_RCF_analysis.addMenu("File")
        
        # Crea il comando "Open calibration file" come sotto-voce di "File"
        open_calibration_act = QAction("Open calibration file", self)
        open_calibration_act.triggered.connect(self.open_calibration_file)
        file_menu.addAction(open_calibration_act)
        
        # Crea il comando "Open image" come sotto-voce di "File"
        open_radiochromic_image_act = QAction("Open image", self)
        open_radiochromic_image_act.triggered.connect(self.open_radiochromic_image)
        file_menu.addAction(open_radiochromic_image_act)
        
        # Crea il comando "Chiudi" nel menu "File"
        close_action = QAction("Close application", self)
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)
        
        # Crea la voce "Tools" nella barra del menù
        Tools_menu = menubar_RCF_analysis.addMenu("Tools")
        
        # Crea il comando "Plot calibration" come sotto_voce di "Tools"
        plot_calibration_act = QAction("Plot calibration curve", self)
        plot_calibration_act.triggered.connect(self.plot_calibration_curve)
        Tools_menu.addAction(plot_calibration_act)
        
        
        
        # ========================  APPLICATION LAYOUT ========================
        
        # Global vertical layout 
        #
        v_box_global = QVBoxLayout()
        
        # Space for the calibration file text
        #
        self.calibration_file_text = QTextEdit()
        # Imposta le dimensioni desiderate
        self.calibration_file_text.setFixedSize(500, 200) 
        self.calibration_file_text.setAlignment(Qt.AlignTop)
        
        # Label of the calibration box text
        #
        label_calibration_file_text = QLabel("CALIBRATION FILE")
        
        # Crea un QLabel vuoto per visualizzare l'immagine
        #
        self.image_container = QLabel()
        # Imposta le dimensioni desiderate
        self.image_container.setFixedSize(400, 400)  
        # Imposta l'allineamento al centro (opzionale)
        self.image_container.setAlignment(Qt.AlignTop)
        self.image_container.setStyleSheet("""background-color: skyblue;
                                           color: white;
                                           border-style: outset;
                                           border-width: 3px;
                                           border-radius: 5px;
                                           font: bold 24px 'Times New Roman';
                                           qproperty-alignment: AlignCenter""")
        
        # Label of the image container
        #
        label_image_container = QLabel("LOADED IMAGE")

        # Logs box space
        #
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFixedSize(1000, 150) 

        # Label of the log box
        label_log_text = QLabel("APPLICATION LOG")
                
        v_box_calibration_text = QVBoxLayout()
        v_box_calibration_text.addWidget(label_calibration_file_text)
        v_box_calibration_text.addStretch(2)
        v_box_calibration_text.addWidget(self.calibration_file_text)
            
        v_box_image = QVBoxLayout()
        v_box_image.addWidget(label_image_container)
        v_box_image.addStretch()
        v_box_image.addWidget(self.image_container)
        
        v_box_log = QVBoxLayout()
        v_box_log.addWidget(label_log_text)
        v_box_log.addWidget(self.log_text)
        
        h_box_1 = QHBoxLayout()
        h_box_1.addLayout(v_box_calibration_text)
        h_box_1.addStretch()
        h_box_1.addLayout(v_box_image)
  
        v_box_global.addLayout(h_box_1)
        v_box_global.addStretch()
        v_box_global.addLayout(v_box_log)
       
        # Set main layout of the window
        central_widget = QWidget(self)
        central_widget.setLayout(v_box_global)
        self.setCentralWidget(central_widget)


        # =======================  OPEN CALIBRATION FILE ======================
        # This function ask to open a calibration file already prepared, 
        # store the 
        # fit parameters in the thre variables P1, P2 and P3 and show the same
        # file in a text box
        #
    def open_calibration_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly

        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self, "Open a Calibration File", "", "Text Files (*.txt);;All Files (*)", options=options
        )

        if file_path:
            try:
                with open(file_path, 'r') as file:
                    # Leggi le linee del file
                    file_content = file.readlines()
                    # Aggiungi le operazioni al log
                    log_message = "Calibration file opened: {}\n".format(file_path)

                if len(file_content) >= 5:  # Verifica che ci siano almeno 5 linee nel file
                    fifth_line = file_content[4].strip()  # Estrai la quinta linea e rimuovi spazi vuoti

                    # Dividi la quinta linea in tre numeri utilizzando lo spazio come separatore
                    numbers = fifth_line.split()
                
                    if len(numbers) == 3:  # Verifica che ci siano esattamente 3 numeri
                        P1, P2, P3 = map(float, numbers)
                        # Ora puoi utilizzare le variabili P1, P2 e P3
                        print("Fitting parameter P1:", P1)
                        print("Fitting parameter P2:", P2)
                        print("Fitting parameter P3:", P3)
                        log_message += "Three parameters have been defined for the 3-order polynomial fitting curve:\n \
                        P3*X^3 + P2*X^2 + P1*X \n"
                            
                        log_message += "The parameters are:\nP1 = " + str(P1) +"\n" \
                            + "P2 = " + str(P2) +"\n" + "P3 = " + str(P3) +"\n"
                        

                    else:
                        self.show_error_dialog("La quinta linea del file non contiene 3 numeri.")
                else:
                    self.show_error_dialog("Il file non ha abbastanza linee.")
            except Exception as e:
                self.show_error_dialog(str(e))
                
        if file_path:
            try:
                with open(file_path, 'r') as file:
                    file_content = file.read()
                self.calibration_file_text.setPlainText(file_content)
       
            except Exception as e:
                self.show_error_dialog(str(e))
                
            # Aggiungi altre operazioni al log
            self.log_text.append(log_message)
    
            
    
    # ======================= OPEN RADIOCHROMIC IMAGE =========================
    # This function: 
    # 1.- Open the TIFF file containing the image of a radiochromic film
    # 2.- Show it in a window inside the application and allows for
    #       performing operations on it
    #
    def open_radiochromic_image(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self, 
            "Open a Radiochromic Image", 
            "", 
            "TIFF Files (*.tiff *.tif);;All Files (*)", 
            options=options
            )
        
        if file_path:
            try:
                # Crea una nuova istanza della finestra ImageViewer
                image_viewer = ImageViewer(file_path)
                
                # Mostra la finestra dell'immagine
                image_viewer.exec_()
                
            except Exception as e:
                    self.show_error_dialog(str(e))
      
      # ======================= PLOT CALIBRTION CURVE =========================
              
    def plot_calibration_curve(self):
        print("Ciao Ciao")

        import matplotlib.pyplot as plt
       # dose = list(0,1,1000)
       # od = 
        
        figure20, ax = plt.subplots(figsize=(6, 6))
        plt.plot()

        
                
    def show_error_dialog(self, message):
        error_dialog = QDialog(self)
        error_dialog.setWindowTitle("Error")
        error_dialog.setGeometry(700, 300, 400, 150)
        error_label = QLabel(message)
        error_layout = QVBoxLayout()
        error_layout.addWidget(error_label)
        error_dialog.setLayout(error_layout)
        error_dialog.exec_()
        