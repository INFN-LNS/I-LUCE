import cv2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.patches import Rectangle
from PyQt5.QtWidgets import QSizePolicy, QDialog, QVBoxLayout

class ImageProcessor:
    def __init__(self):
        self.image_canvas = None  # Riferimento alla canvas dell'immagine
        self.image_window = None  # Riferimento alla finestra dell'immagine
        self.image = None
        self.roi_rect = None
        self.roi_mean_value = None

    def show_image(self, image):
        # Memorizza l'immagine per l'accesso successivo
        self.image = image

        # Mostra l'immagine nel canale rosso con matplotlib
        fig, ax = plt.subplots()
        ax.imshow(image, cmap='gray')
        ax.set_title("Red Channel Image")
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")

        # Aggiungi la patch della ROI alla figura
        self.roi_rect = Rectangle((0, 0), 1, 1, edgecolor='r', facecolor='none')
        ax.add_patch(self.roi_rect)

        # Connetti il gestore degli eventi del mouse
        self.image_canvas = FigureCanvas(fig)
        self.image_canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.image_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Crea una finestra per l'immagine e mostra la canvas
        self.image_window = QDialog()
        layout = QVBoxLayout()
        layout.addWidget(self.image_canvas)
        self.image_window.setLayout(layout)
        self.image_window.show()

    def on_mouse_press(self, event):
        if event.button == 1:  # Solo il pulsante sinistro
            # Rimuovi la vecchia patch ROI se esiste
            if hasattr(self, 'roi_rect'):
                self.roi_rect.remove()

            self.start_x, self.start_y = int(event.xdata), int(event.ydata)
            self.roi_rect = Rectangle((self.start_x, self.start_y), 1, 1, edgecolor='r', facecolor='none')
            self.image_canvas.figure.axes[0].add_patch(self.roi_rect)
            self.image_canvas.draw()

            # Connetti l'evento di trascinamento
            self.drag_cid = self.image_canvas.mpl_connect('motion_notify_event', self.on_mouse_drag)
            # Connetti l'evento di rilascio del pulsante
            self.release_cid = self.image_canvas.mpl_connect('button_release_event', self.on_mouse_release)
        
    def on_mouse_drag(self, event):
        if hasattr(self, 'start_x') and hasattr(self, 'start_y'):
            cur_x, cur_y = int(event.xdata), int(event.ydata)
            w = cur_x - self.start_x
            h = cur_y - self.start_y

        if w > 0 and h > 0:
            # Aggiorna la posizione e le dimensioni della patch della ROI durante il trascinamento
            self.update_roi_rect(self.start_x, self.start_y, w, h)

    def on_mouse_release(self, event):
        # Disconnetti gli eventi di trascinamento e rilascio del pulsante
        self.image_canvas.mpl_disconnect(self.drag_cid)
        self.image_canvas.mpl_disconnect(self.release_cid)

        # Calcola e mostra il valore medio della ROI
        w = int(event.xdata) - self.start_x
        h = int(event.ydata) - self.start_y
        self.calculate_roi_mean_value(self.start_x, self.start_y, w, h)
        
        
    def update_roi_rect(self, x, y, w, h):
        # Aggiorna la posizione e le dimensioni della patch della ROI
        self.roi_rect.set_xy((x, y))
        self.roi_rect.set_width(w)
        self.roi_rect.set_height(h)
        self.image_canvas.draw()

    def calculate_roi_mean_value(self, x, y, w, h):
        # Calcola il valore medio dei pixel nella ROI
        roi = self.image[y:y+h, x:x+w]
        mean_value = roi.mean()
        self.roi_mean_value = mean_value
        print(f"ROI Mean Value: {mean_value}")
