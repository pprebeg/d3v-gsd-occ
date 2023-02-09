import os.path

from PySide6.QtWidgets import QApplication, QMenu, QFormLayout,QWidget,QHeaderView,QSlider,QLineEdit
from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QVBoxLayout,QHBoxLayout,QTableView,QTextEdit,QLabel
from PySide6.QtCore import Slot,Qt,QMetaObject,QCoreApplication
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRect
from PySide6.QtGui import QColor, QPainter
#from PySide6.QtCharts import QtCharts
from PySide6.QtWidgets import QFileDialog,QDialogButtonBox,QProgressDialog,QMenuBar
from PySide6.QtCore import SIGNAL,SLOT
from typing import Dict,List
from uuid import uuid4

import hullformdir.hullgeneratorform

try:
    # d3v imports
    from signals import Signals
    from commands import Command
    from iohandlers import IOHandler
    from core import geometry_manager as manager
    from core.geometry import Geometry
    from mainwin import MainFrame
    from glwin import GlWin
    from application import App
    # d3v-gsd
    from hullformdir.hullform import *
    from hullformdir.hullgeneratorform import HullGeneratorForm
    from hydrostatics_gui import HydrostaticsGUI
    from ship_stability_gui import ShipStabilityGUI
    from hullmoddir.occhullform import OCCHullform
    from hullform_command import HullFormCommand
except BaseException as error:
    print('An exception occurred: {}'.format(error))
except:
    print('Unknown exception occurred during signals connection')

def get_menu_by_name(mb:QMenuBar,name:str)->QMenu:
    actions = mb.actions()
    for action in actions:
        if action.menu():
            if action.text()== name:
                return action
    return None

class HullmodCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        self.app.registerIOHandler(HullmodImporter())
        self.hfcom:HullmodCommand= None
        for com in self.app.commands:
            if isinstance(com,HullFormCommand):
                self.hfcom = com
                break
        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
            manager.geometry_removed.connect(self.onGeometryRemoved)
            manager.visible_geometry_changed.connect(self.onVisibleGeometryChanged)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')

        self.menuOCCForm = QMenu("&OCC Form")
        self.hfcom.menuMain.insertMenu(self.hfcom.menuHullGenForm.menuAction() ,self.menuOCCForm)

        menuHullgenModify = self.menuOCCForm.addAction("&Modify form")
        menuHullgenModify.triggered.connect(self.onModifyHullgenForm)

    def onModifyHullgenForm(self):
        if isinstance(self.hfcom.active_hull_form, OCCHullform):
            self.hfcom.active_hull_form.modify_form()

    @Slot()
    def onVisibleGeometryChanged(self, visible:List[Geometry], loaded:List[Geometry], selected:List[Geometry]):
        for g in visible:
            pass

    @Slot()
    def onSelectedGeometryChanged(self, visible: List[Geometry], loaded: List[Geometry], selected: List[Geometry]):
        pass

    @Slot()
    def onGeometryCreated(self, geometries:List[Geometry]):
        pass


    @Slot()
    def onGeometryRemoved(self, geometries:List[Geometry]):
        for g in geometries:
            if isinstance(g, OCCHullform):
                pass

    @property
    def app(self)->App:
        return self._app

    @property
    def mainwin(self)->MainFrame:
        return self.app.mainFrame

    @property
    def glwin(self)->GlWin:
        return self.mainwin.glWin

class HullmodImporter(IOHandler):
    def __init__(self,force_import=False):
        super().__init__()

    def import_geometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_no_ext, file_extension = os.path.splitext(os.path.basename(fileName))
        hf=None
        if file_extension in self.getImportFormats():
            hf = OCCHullform(fileName, filename_no_ext)
        if hf is not None:
            return hf

    def export_geometry(self, fileName, geometry2export):
        if isinstance(geometry2export,HullForm):
            geometry2export.exportGeometry(fileName)
        om.write_mesh(geometry2export.mesh,fileName)
        pass

    def getExportFormats(self):
       return (".hoc",".igs")

    def getImportFormats(self):
        return (".hoc",".igs")



def createCommand():
    return HullmodCommand()