#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 21:13:52 2023

@author: cirrone
"""
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QStackedWidget, QDialog
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt

class MyGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("La Tua Applicazione")
        self.setGeometry(100, 100, 800, 600)  # Imposta le dimensioni della finestra principale

        # Imposta il colore di sfondo
        self.setStyleSheet("background-color: darkgreen; color: white; font-size: 16px; background-color: darkgreen;")

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        main_layout = QVBoxLayout()

        # Intestazione con logo, nome istituto, nome applicazione e versione
        header_layout = QVBoxLayout()

        logo_label = QLabel()
        logo_pixmap = QPixmap("LOGO_INFN_NEWS_sito.jpeg")  # Inserisci il percorso del tuo logo
        logo_pixmap = logo_pixmap.scaledToWidth(300)
        logo_label.setPixmap(logo_pixmap)
        logo_label.setAlignment(Qt.AlignCenter)

        institute_label = QLabel("Nome del tuo istituto")
        institute_label.setAlignment(Qt.AlignCenter)

        app_label = QLabel("La Tua Applicazione v1.0")  # Sostituisci "v1.0" con la tua versione
        app_label.setAlignment(Qt.AlignCenter)

        name_label = QLabel("Il tuo nome")
        name_label.setAlignment(Qt.AlignCenter)

        header_layout.addWidget(logo_label)
        header_layout.addWidget(institute_label)
        header_layout.addWidget(app_label)
        header_layout.addWidget(name_label)

        main_layout.addLayout(header_layout)

        # Barra laterale con pulsanti per le varie sezioni
        sidebar_layout = QVBoxLayout()

        # Regola le dimensioni e l'allineamento dei pulsanti
        button_style = "font-size: 14px; text-align: left; background-color:#003300"

        calibration_button = QPushButton("Calibration")
        calibration_button.setStyleSheet(button_style)
        calibration_button.setFixedSize(200, 50)
        single_image_analysis_button = QPushButton("Single Image Analysis")
        single_image_analysis_button.setStyleSheet(button_style)
        single_image_analysis_button.setFixedSize(200, 50)
        unfolding_button = QPushButton("Unfolding")
        unfolding_button.setStyleSheet(button_style)
        unfolding_button.setFixedSize(200, 50)
        exit_button = QPushButton("Esci")
        exit_button.setStyleSheet(button_style)
        exit_button.setFixedSize(200, 50)

        # Collega i pulsanti alle rispettive azioni
        calibration_button.clicked.connect(self.open_calibration)
        single_image_analysis_button.clicked.connect(self.open_single_image_analysis)
        unfolding_button.clicked.connect(self.open_unfolding)
        exit_button.clicked.connect(self.close)

        sidebar_layout.addWidget(calibration_button)
        sidebar_layout.addWidget(single_image_analysis_button)
        sidebar_layout.addWidget(unfolding_button)
        sidebar_layout.addWidget(exit_button)

        main_layout.addLayout(sidebar_layout)

        # Stack di widget per le pagine delle diverse sezioni
        self.pages = QStackedWidget()

        calibration_page = QWidget()
        single_image_analysis_page = QWidget()
        unfolding_page = QWidget()

        self.pages.addWidget(calibration_page)
        self.pages.addWidget(single_image_analysis_page)
        self.pages.addWidget(unfolding_page)

        main_layout.addWidget(self.pages)

        central_widget.setLayout(main_layout)

    def open_calibration(self):
        calibration_window = QDialog(self)
        calibration_window.setWindowTitle("Calibration")
        calibration_window.setGeometry(200, 200, 400, 300)  # Imposta le dimensioni della finestra
        calibration_window.exec_()

    def open_single_image_analysis(self):
        single_image_analysis_window = QDialog(self)
        single_image_analysis_window.setWindowTitle("Single Image Analysis")
        single_image_analysis_window.setGeometry(200, 200, 400, 300)  # Imposta le dimensioni della finestra
        single_image_analysis_window.exec_()

    def open_unfolding(self):
        unfolding_window = QDialog(self)
        unfolding_window.setWindowTitle("Unfolding")
        unfolding_window.setGeometry(200, 200, 400, 300)  # Imposta le dimensioni della finestra
        unfolding_window.exec_()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MyGUI()
    window.show()
    sys.exit(app.exec_())
