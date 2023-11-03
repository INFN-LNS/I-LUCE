#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 22:32:49 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QPushButton, QFileDialog, QTextEdit, QMainWindow, QApplication, QWidget, QMainWindow

class SingleImageAnalysisWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Single Image Analysis")
        self.setGeometry(600, 100, 1000, 900)  # Imposta le dimensioni della finestra
        
        # Imposta il colore di sfondo
        self.setStyleSheet("background-color:white;")


        layout = QVBoxLayout()
        
        
        # Aggiungi una barra di menu verticale allineata a sinistra
        menu_layout = QVBoxLayout()
        
        # Regola le dimensioni e l'allineamento dei pulsanti
        button_style = "font-size: 14px; text-align: left; background-color:#003300"

        # Aggiungi il pulsante "Open Calibration File" al menu
        open_calibration_button = QPushButton("Open Calibration File")
        open_calibration_button.setStyleSheet(button_style)
        open_calibration_button.setFixedSize(200, 50)
        open_calibration_button.setStyleSheet("background-color: #003300; color: white;")

        
        menu_layout.addWidget(open_calibration_button)
        layout.addLayout(menu_layout)
        
        
        # Aggiungi un'area di testo per visualizzare il contenuto del file
        self.text_edit = QTextEdit()
        layout.addWidget(self.text_edit)
        
        # Collega il pulsante "Open Calibration File" al metodo di apertura del file
        open_calibration_button.clicked.connect(self.open_calibration_file)

        self.setLayout(layout)


        # Aggiungi un pulsante di chiusura
        close_button = QPushButton("Chiudi")
        layout.addWidget(close_button)
    
        close_button.clicked.connect(self.accept)  # Chiudi la finestra quando il pulsante viene premuto

        self.setLayout(layout)

    def open_calibration_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly

        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Open a Calibration File", "", "Text Files (*.txt);;All Files (*)", options=options)

        if file_path:
            # Apre e legge il contenuto del file
            with open(file_path, 'r') as file:
                file_content = file.read()

            # Visualizza il contenuto del file nell'area di testo
            self.text_edit.setPlainText(file_content)



if __name__ == '__main__':
    app = QApplication([])
    window = SingleImageAnalysisWindow()
    window.show()
    app.exec_()





        
 