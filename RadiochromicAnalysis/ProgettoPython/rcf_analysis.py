#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 22:43:47 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenuBar, QAction, QWidget, QLabel, QWidgetAction


class RCFAnalysisWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("RCF Analysis")
        self.setGeometry(200, 200, 400, 300)
        
        # Crea la barra dei menu
        menubar_RCF_analysis = self.menuBar()
        
        # Crea il menu "File" nella barra dei menu
        file_menu_RCF_analysis = menubar_RCF_analysis.addMenu("File_2")
        
        # Crea il comando "Chiudi" nel menu "File"
        close_action = QAction("Chiudi", self)
        close_action.triggered.connect(self.close)
        file_menu_RCF_analysis.addAction(close_action)
        
 

