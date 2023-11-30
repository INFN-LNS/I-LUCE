#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 22:43:47 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenuBar, QAction, 
                             QWidget, QLabel, QWidgetAction, QFileDialog,
                             QTextEdit, QVBoxLayout, QDialog, QHBoxLayout,
                             QPushButton)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QGraphicsScene, QGraphicsView


import numpy as np
import matplotlib.pyplot as plt
import io
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

import cv2

from image_processor import ImageProcessor


class PlotDialog(QDialog):
    def __init__(self, x1, y1, x2, y2):
        super(PlotDialog, self).__init__()

        # Imposta il layout della finestra
        layout = QVBoxLayout()

        # Crea un oggetto FigureCanvasQtAgg per ospitare il grafico matplotlib
        self.canvas = FigureCanvasQTAgg(Figure())
        layout.addWidget(self.canvas)

        # Imposta il grafico sulla canvas
        self.ax = self.canvas.figure.add_subplot(111)
        self.ax.plot(x1, y1,'o', x2, y2,'--')
        self.ax.set_title('Calibration points and polynomial fit')
        self.ax.set_xlabel('Net optical density')
        self.ax.set_ylabel('Dose')

        # Aggiungi layout alla finestra
        self.setLayout(layout)
        self.setGeometry(100, 100, 800, 600)
        self.setWindowTitle('Plot Calibration Curve')
        

class RCFAnalysisWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("RCF Analysis")
        self.setGeometry(800, 200, 1000, 800)
        
        # Aggiungi un'istanza di ImageProcessor
        self.image_processor = ImageProcessor()
        
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
          
        # Crea il comando "Chiudi" nel menu "File"
        close_action = QAction("Close application", self)
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)
        
        # Crea il comando "Open image to analyse" come sotto-voce di "File"
        #
        open_radiochromic_image_act = QAction("Open image to analyse", self)
        open_radiochromic_image_act.triggered.connect(self.open_radiochromic_image)
        file_menu.addAction(open_radiochromic_image_act)
        
        # Crea la voce "Tools" nella barra del menù
        Tools_menu = menubar_RCF_analysis.addMenu("Tools")
        
        # Crea il comando "Plot calibration" come sotto_voce di "Tools"
        plot_calibration_act = QAction("Plot calibration curve", self)
        plot_calibration_act.triggered.connect(self.plot_polynomial_curve)
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
        
        # Empty QLabel to show the dose profile 
        #
        self.profile_container = QLabel()
        # Imposta le dimensioni desiderate
        self.profile_container.setFixedSize(600, 400)  
        # Imposta l'allineamento al centro (opzionale)
        self.profile_container.setAlignment(Qt.AlignTop)
        self.profile_container.setStyleSheet("""background-color: white;
                                           color: white;
                                           border-style: outset;
                                           border-width: 3px;
                                           border-radius: 5px;
                                           font: bold 24px 'Times New Roman';
                                           qproperty-alignment: AlignCenter""")
        
        # Label of the figure container
        #
        label_profile_container = QLabel("Dose profile curve")

        # Logs box space
        #
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFixedSize(1000, 150)
        self.log_text.setStyleSheet("""background-color: black;
                                           color: green;
                                           border-style: outset;
                                           border-width: 3px;
                                           border-radius: 5px;
                                           font: bold 20px 'Times New Roman';""")
                                          

        # Label of the log box
        label_log_text = QLabel("APPLICATION LOG")
       
        
        # Crea un layout verticale per la scritta "Calibration file" e il QTextEdit
        v_box_calibration_text = QVBoxLayout()

        # Aggiungi la scritta "Calibration file"
        v_box_calibration_text.addWidget(label_calibration_file_text, alignment=Qt.AlignTop)

        # Aggiungi lo spazio per il QTextEdit
        v_box_calibration_text.addStretch(2)
        v_box_calibration_text.addWidget(self.calibration_file_text)

        # Imposta i margini del layout verticale
        v_box_calibration_text.setContentsMargins(0, 0, 0, 200)
            
        v_box_profile = QVBoxLayout()
        v_box_profile.addWidget(label_profile_container)
        v_box_profile.addStretch()
        v_box_profile.addWidget(self.profile_container)
        
        v_box_log = QVBoxLayout()
        v_box_log.addWidget(label_log_text)
        v_box_log.addWidget(self.log_text)
        
        h_box_1 = QHBoxLayout()
        h_box_1.addLayout(v_box_calibration_text)
        h_box_1.addStretch()
        h_box_1.addLayout(v_box_profile)
  
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
                    
                    # Estrai i dati dalle colonne Dose e Net OD
                    self.dose_values = []
                    self.net_od_values = []
                    for line in file_content:
                        if line.startswith("end"):
                            break
                        columns = line.split()
                        if len(columns) == 2:  # Verifica che ci siano esattamente 2 colonne
                            self.dose_values.append(float(columns[0]))
                            self.net_od_values.append(float(columns[1]))
                            
                
                    if len(numbers) == 3:  # Verifica che ci siano esattamente 3 numeri
                        P1, P2, P3 = map(float, numbers)
                   
                        log_message += "Three parameters have been defined for the 3-order polynomial fitting curve:\n \
                        P3*X^3 + P2*X^2 + P1*X \n"
                            
                        log_message += "The parameters are:\nP1 = " + str(P1) +"\n" \
                            + "P2 = " + str(P2) +"\n" + "P3 = " + str(P3) +"\n\n"
                            
                        log_message += "The dose values are: " + str(self.dose_values) + "\n"
                        log_message += "The net-OD values are: " + str(self.net_od_values) + "\n"
                        

                        if len(numbers) == 3:
                            self.P1, self.P2, self.P3 = map(float, numbers)
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
    #
    
    def open_radiochromic_image(self):
        # Chiedi all'utente di selezionare un'immagine utilizzando QFileDialog
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Image File", "", "Image Files (*.tif *.png *.jpg);;All Files (*)", options=options)

        # Verifica se l'utente ha selezionato un file
        if file_path:
            # Carica l'immagine utilizzando OpenCV
            image = cv2.imread(file_path)

            # Mostra solo il canale rosso
            red_channel = image[:, :, 2]

            # Utilizza ImageProcessor per mostrare l'immagine
            self.image_processor.show_image(red_channel)

            # Aggiungi operazioni al log
            log_message = f"Radiochromic image opened: {file_path}\n"
            self.log_text.append(log_message)
      
      # ======================= PLOT CALIBRTION CURVE =========================
              
    def plot_polynomial_curve(self):
        # Verifica se il file di calibrazione è stato caricato
        if not hasattr(self, 'P1') or not hasattr(self, 'P2') or not hasattr(self, 'P3'):
            self.show_error_dialog('Calibration file not loaded. Please load a calibration file first.')
            return

        # Generazione del vettore x da 1 a 1000
        x = np.arange(0, 1, 0.1)

        # Calcolo dei valori y del polinomio
        #
        y = self.P3 * x**3 + self.P2 * x**2 + self.P1

        # Creazione della finestra di dialogo per il grafico
        plot_dialog = PlotDialog(self.net_od_values, self.dose_values, x, y)
        plot_dialog.exec_()
                
    def show_error_dialog(self, message):
        error_dialog = QDialog(self)
        error_dialog.setWindowTitle("Error")
        error_dialog.setGeometry(700, 300, 400, 150)
        error_label = QLabel(message)
        error_layout = QVBoxLayout()
        error_layout.addWidget(error_label)
        error_dialog.setLayout(error_layout)
        error_dialog.exec_()
        