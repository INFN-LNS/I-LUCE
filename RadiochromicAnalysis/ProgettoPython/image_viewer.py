#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:33:42 2023

@author: cirrone
"""
from PyQt5.QtWidgets import QDialog, QLabel, QVBoxLayout
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt

class ImageViewer(QDialog):
    def __init__(self, image_path):
        super().__init__()

        self.setWindowTitle("Image Viewer")
        self.setGeometry(100, 100, 400, 400)

        # Crea un QLabel per visualizzare l'immagine
        image_label = QLabel()
        image_label.setAlignment(Qt.AlignCenter)

        # Carica l'immagine da file
        image = QPixmap(image_path)
        image_label.setPixmap(image)

        # Crea un layout verticale per la finestra dell'immagine
        layout = QVBoxLayout()
        layout.addWidget(image_label)

        # Imposta il layout per la finestra dell'immagine
        self.setLayout(layout)

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication

    app = QApplication(sys.argv)
    image_path = "path_all_immagine"
    viewer = ImageViewer(image_path)
    viewer.show()
    sys.exit(app.exec_())
