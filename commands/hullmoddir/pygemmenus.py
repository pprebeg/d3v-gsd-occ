from ast import literal_eval
from numpy import array
from PySide6.QtCore import Qt


from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDateEdit,
    QDateTimeEdit,
    QDial,
    QDoubleSpinBox,
    QFontComboBox,
    QLabel,
    QLCDNumber,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QTimeEdit,
    QVBoxLayout,
    QWidget,
)

def text_to_nparray(text):
    text = "["+text+"]"
    list = literal_eval(text)
    return array(list)


def text_to_list(text):
    text = "["+text+"]"
    return literal_eval(text)


class CreateFFDBox_menu(QWidget):
    def __init__(self, Hullform):
        super().__init__()
        self.Hullform=Hullform

        self.setWindowTitle("FFD box creation")

        self.layout = QVBoxLayout(self)
        #create widgets:
        self.Guide_label = QLineEdit("Make sure all values are seperated with a comma \",\"")

        #box_center = np.array([50000, 6000, 3000]), box_length = np.array([5000, 5000, 5000]), n_control_points = [3, 3,3]

        self.FFD_Box_center_label = QLabel("Set FFD Box center:")
        self.FFD_Box_center_input = QLineEdit("50000, 6000, 3000")        #self.lineedit.setPlaceholderText()

        self.FFD_Box_dims_label = QLabel("Set FFD Box dimensions:")
        self.FFD_Box_dims_input = QLineEdit("5000, 5000, 5000")

        self.FFD_Box_ncontrol_label = QLabel("Set numbers of FFD box control points:")
        self.FFD_Box_ncontrol_input = QLineEdit("2, 2, 2")

        self.FFD_Box_rotation_label = QLabel("Set FFD box rotation angles:")
        self.FFD_Box_rotation_input = QLineEdit("0, 0, 0")

        self.FFD_Box_knots_label = QLabel("Set number of knots to add (increases accuracy of deformation, but is slower):")
        self.FFD_Box_knots_input = QLineEdit("0, 0, 0")

        self.create_button = QPushButton("Create Box")

        #adding widgets to window
        self.layout.addWidget(self.Guide_label)
        self.layout.addWidget(self.FFD_Box_center_label)
        self.layout.addWidget(self.FFD_Box_center_input)
        self.layout.addWidget(self.FFD_Box_dims_label)
        self.layout.addWidget(self.FFD_Box_dims_input)
        self.layout.addWidget(self.FFD_Box_ncontrol_label)
        self.layout.addWidget(self.FFD_Box_ncontrol_input)
        self.layout.addWidget(self.FFD_Box_rotation_label)
        self.layout.addWidget(self.FFD_Box_rotation_input)
        self.layout.addWidget(self.FFD_Box_knots_label)
        self.layout.addWidget(self.FFD_Box_knots_input)
        self.layout.addWidget(self.create_button)

        self.create_button.setCheckable(True)
        self.create_button.clicked.connect(self.on_create_button_clicked)

    def on_create_button_clicked(self):    #when create button is clicked
        self._box_center = text_to_nparray(self.FFD_Box_center_input.text())
        self._box_length = text_to_nparray(self.FFD_Box_dims_input.text())
        self._n_control_points = text_to_nparray(self.FFD_Box_ncontrol_input.text())
        self._box_rotation_angles = text_to_nparray(self.FFD_Box_rotation_input.text())
        self._knots_to_add = text_to_nparray(self.FFD_Box_knots_input.text())


        self.Hullform.make_ffd_volume(self._box_center, self._box_length, self._n_control_points, self._box_rotation_angles, self._knots_to_add)
        self.Hullform.make_ffd_box_mesh()

        self.delete_window()

    def show_window(self):
        self.show()

    def delete_window(self):
        self.destroy()


class DeformFFDBox_menu(QWidget):
    def __init__(self, Hullform, scale_limit = 3, tick_density = 5):
        super().__init__()
        self.Hullform=Hullform
        self._scale_limit = scale_limit
        self._tick_density = tick_density       #how many intervals are there between integers

        self._active_deformation = True                 #Qt.CheckState.Checked
        self._cpoint_id = [1,1,1]
        self._total_ticks = int(self._tick_density * self._scale_limit)
        self._interval_len = 1/self._tick_density
        start_pos = self._tick_density
        self._x_slider_tick = start_pos
        self._y_slider_tick = start_pos
        self._z_slider_tick = start_pos

        self._x_scaling = start_pos * self._interval_len
        self._y_scaling = start_pos * self._interval_len
        self._z_scaling = start_pos * self._interval_len


        self.setWindowTitle("FFD box deformation")

        self.layout = QVBoxLayout(self)
        #create widgets:
        self.active_deformation_check_label = QLabel("Active Deformation:")
        self.active_deformation_checkbox = QCheckBox()
        self.active_deformation_checkbox.setChecked(self._active_deformation)

        self.Guide_label = QLabel("Make sure all values are seperated with a comma \",\"")

        self.Deformation_cpoint_id_label = QLabel("Select FFD Box Control Point ID:")
        self.Deformation_cpoint_id_input = QLineEdit("1, 1, 1")        #self.lineedit.setPlaceholderText()

        self.Deformation_slider_x_label = QLabel(f"Control Point X scaling: {self._x_scaling:.2f}")
        self.Deformation_slider_x = QSlider(Qt.Horizontal)          #interval im je samo integer........
        self.Deformation_slider_x.setMinimum(0)
        self.Deformation_slider_x.setMaximum(self._total_ticks)
        self.Deformation_slider_x.setSliderPosition(self._x_slider_tick)
        self.Deformation_slider_x.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.Deformation_slider_y_label = QLabel(f"Control Point Y scaling: {self._y_scaling:.2f}")
        self.Deformation_slider_y = QSlider(Qt.Horizontal)
        self.Deformation_slider_y.setMinimum(0)
        self.Deformation_slider_y.setMaximum(self._total_ticks)
        self.Deformation_slider_y.setSliderPosition(self._y_slider_tick)
        self.Deformation_slider_y.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.Deformation_slider_z_label = QLabel(f"Control Point Z scaling: {self._z_scaling:.2f}")
        self.Deformation_slider_z = QSlider(Qt.Horizontal)
        self.Deformation_slider_z.setMinimum(0)
        self.Deformation_slider_z.setMaximum(self._total_ticks)
        self.Deformation_slider_z.setSliderPosition(self._z_slider_tick)
        self.Deformation_slider_z.setTickPosition(QSlider.TickPosition.TicksBothSides)

        self.deform_button = QPushButton("Deform Box")

        #adding widgets to window
        self.layout.addWidget(self.active_deformation_check_label)
        self.layout.addWidget(self.active_deformation_checkbox)


        self.layout.addWidget(self.Guide_label)
        self.layout.addWidget(self.Deformation_cpoint_id_label)
        self.layout.addWidget(self.Deformation_cpoint_id_input)
        self.layout.addWidget(self.Deformation_slider_x_label)
        self.layout.addWidget(self.Deformation_slider_x)
        self.layout.addWidget(self.Deformation_slider_y_label)
        self.layout.addWidget(self.Deformation_slider_y)
        self.layout.addWidget(self.Deformation_slider_z_label)
        self.layout.addWidget(self.Deformation_slider_z)
        self.layout.addWidget(self.deform_button)

        #connect signals:
        self.active_deformation_checkbox.stateChanged.connect(self.on_checkbox_toggle)
        self.Deformation_cpoint_id_input.selectionChanged.connect(self.on_cpoint_change)
        self.Deformation_slider_x.sliderMoved.connect(self.on_x_slider_change)
        self.Deformation_slider_y.sliderMoved.connect(self.on_y_slider_change)
        self.Deformation_slider_z.sliderMoved.connect(self.on_z_slider_change)
        self.deform_button.clicked.connect(self.on_deform_button_clicked)


    def on_checkbox_toggle(self):
        if self.active_deformation_checkbox.checkState() == Qt.CheckState.Checked:
            self._active_deformation = True
        else:
            self._active_deformation = False

    def on_deform_button_clicked(self):    #when deform button is clicked
        move_vector = [self._x_scaling, self._y_scaling, self._z_scaling]
        self.Hullform.move_ffd_pole(ffd_id = 0, pole_id = self._cpoint_id, move_vector = move_vector)       #kasnije promijeni ffd_id, ovo je samo za test
        self.Hullform.make_ffd_box_mesh()
        self.Hullform.ffd_deform_surfaces()

    def on_cpoint_change(self):
        id = text_to_list(self.Deformation_cpoint_id_input.text())
        self._cpoint_id = id

    def on_x_slider_change(self):
        self._x_slider_pos = self.Deformation_slider_x.sliderPosition()
        self._x_scaling = self._x_slider_pos * self._interval_len
        self.Deformation_slider_x_label.setText(f"Control Point X scaling: {self._x_scaling:.2f}")
        if self._active_deformation:
            self.on_deform_button_clicked()

    def on_y_slider_change(self):
        self._y_slider_pos = self.Deformation_slider_y.sliderPosition()
        self._y_scaling = self._y_slider_pos * self._interval_len
        self.Deformation_slider_y_label.setText(f"Control Point Y scaling: {self._y_scaling:.2f}")
        if self._active_deformation:
            self.on_deform_button_clicked()

    def on_z_slider_change(self):
        self._z_slider_pos = self.Deformation_slider_z.sliderPosition()
        self._z_scaling = self._z_slider_pos * self._interval_len
        self.Deformation_slider_z_label.setText(f"Control Point Z scaling: {self._z_scaling:.2f}")
        if self._active_deformation:
            self.on_deform_button_clicked()

    def show_window(self):
        self.show()

    def delete_window(self):
        self.destroy()















