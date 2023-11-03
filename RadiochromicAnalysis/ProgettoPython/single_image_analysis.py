#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 22:32:49 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QPushButton

class SingleImageAnalysisWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Single Image Analysis")
        self.setGeometry(600, 100, 1000, 900)  # Imposta le dimensioni della finestra
        
        # Imposta il colore di sfondo
        self.setStyleSheet("background-color:#003300;")


        layout = QVBoxLayout()

        # Aggiungi i tuoi widget e elementi grafici per la finestra di "Single Image Analysis" qui
        label = QLabel("Contenuto della finestra Single Image Analysis")
        layout.addWidget(label)
        
        # Aggiungi un pulsante di chiusura
        close_button = QPushButton("Chiudi")
        layout.addWidget(close_button)

        close_button.clicked.connect(self.accept)  # Chiudi la finestra quando il pulsante viene premuto


        self.setLayout(layout)