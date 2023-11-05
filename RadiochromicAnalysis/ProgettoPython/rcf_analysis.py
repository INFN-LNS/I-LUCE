#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 22:43:47 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenuBar, QAction, 
                             QWidget, QLabel, QWidgetAction, QFileDialog,
                             QTextEdit, QVBoxLayout, QDialog)


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
        
 
        # ========================  APPLICATION LAYOUT ========================
        
        # Create vertical layout and add widgets 
        v_box = QVBoxLayout()
        
        
        # Here the calibration file is shown
        #
        self.calibration_file_text = QTextEdit()
        # Imposta le dimensioni desiderate
        self.calibration_file_text.setFixedSize(500, 200) 
        
        # LABEL of the calibration box
        title_label_calibration_text = QLabel("CALIBRATION FILE")
       
        # This is a space where the softwre logs shoudl be shown
        #
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFixedSize(1000, 150) 

        # LABEL of the calibration box
        title_label_application_log = QLabel("APPLICATION LOG")

        
        # Add the defined widgets to the vertical layout
        #
        v_box.addWidget(title_label_calibration_text)
        v_box.addWidget(self.calibration_file_text)
        v_box.addStretch()
        v_box.addWidget(title_label_application_log)
        v_box.addWidget(self.log_text)
        
        # Set main layout of the window
        central_widget = QWidget(self)
        central_widget.setLayout(v_box)
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
                        # Ora puoi utilizzare le variabili number1, number2 e number3
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
    # 2.- ...
    #
    def open_radiochromic_image(self):
        
                
    def show_error_dialog(self, message):
        error_dialog = QDialog(self)
        error_dialog.setWindowTitle("Error")
        error_dialog.setGeometry(700, 300, 400, 150)
        error_label = QLabel(message)
        error_layout = QVBoxLayout()
        error_layout.addWidget(error_label)
        error_dialog.setLayout(error_layout)
        error_dialog.exec_()
        