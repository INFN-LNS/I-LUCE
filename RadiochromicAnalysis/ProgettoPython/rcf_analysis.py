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
        self.setGeometry(800, 200, 1000, 800)
        
        # Crea la barra dei menu
        menubar_RCF_analysis = self.menuBar()
        
        # Crea il menu "Open calibration file" nella barra dei menu
        open_menu = menubar_RCF_analysis.addMenu("Open")
        
        # Crea il comando "Open image" come sotto-voce di "Open"
        open_image_action = QAction("Open image", self)
        open_menu.addAction(open_image_action)
        
        
        
        # Crea il comando "Chiudi" nel menu "File"
        close_action = QAction("Chiudi", self)
        close_action.triggered.connect(self.close)
        open_menu.addAction(close_action)
        
 

