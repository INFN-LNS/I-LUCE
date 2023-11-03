import sys
from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QLabel, QPushButton, QFileDialog, QTextEdit, QSizePolicy
)

class SingleImageAnalysisApp(QApplication):
    def __init__(self):
        super().__init__(sys.argv)
        self.window = SingleImageAnalysisWindow()
        self.window.show()
        sys.exit(self.exec_())

class SingleImageAnalysisWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Single Image Analysis")
        self.setGeometry(600, 100, 1000, 900)
        self.setStyleSheet("background-color: white;")

        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        # Menu layout
        menu_layout = QVBoxLayout()
        button_style = "color: white; font-size: 18px; text-align: left; background-color:#003300"

        open_calibration_button = QPushButton("Open Calibration File")
        open_calibration_button.setStyleSheet(button_style)
        open_calibration_button.setFixedSize(200, 50)

        menu_layout.addWidget(open_calibration_button)
        layout.addLayout(menu_layout)

        # RIQUADRO DI TESTO DOVE VIENE MOSTRATO IL CALIBRATION FILE
        self.text_edit = QTextEdit()
        self.text_edit.setFixedSize(300, 100)  # Imposta le dimensioni desiderate
        layout.addWidget(self.text_edit)

        # LOG TXT WINDOW
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        layout.addWidget(self.log_text)

        # CLOSE BUTTON
        close_button = QPushButton("Close this analysis window")
        close_button.setStyleSheet(button_style)
        
        # Aggiungi il layout vuoto espandibile sopra del pulsante Close
        spacer = QLabel()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(spacer)
        layout.addWidget(close_button)

        open_calibration_button.clicked.connect(self.open_calibration_file)
        close_button.clicked.connect(self.accept)
        self.setLayout(layout)

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
                self.text_edit.setPlainText(file_content)
                
                
                # Aggiungi altre operazioni al log
                self.log_text.append(log_message)
               
            except Exception as e:
                self.show_error_dialog(str(e))
                
                
               
    def show_error_dialog(self, message):
        error_dialog = QDialog(self)
        error_dialog.setWindowTitle("Error")
        error_dialog.setGeometry(700, 300, 400, 150)
        error_label = QLabel(message)
        error_layout = QVBoxLayout()
        error_layout.addWidget(error_label)
        error_dialog.setLayout(error_layout)
        error_dialog.exec_()

if __name__ == '__main__':
    app = SingleImageAnalysisApp()
