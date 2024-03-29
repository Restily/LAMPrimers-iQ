# Form implementation generated from reading ui file 'uic.ui'
#
# Created by: PyQt6 UI code generator 6.2.3
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(860, 642)
        MainWindow.setMinimumSize(QtCore.QSize(0, 0))
        MainWindow.setMaximumSize(QtCore.QSize(5000, 2000))
        MainWindow.setStyleSheet("")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.verticalWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalWidget.setMinimumSize(QtCore.QSize(460, 620))
        self.verticalWidget.setMaximumSize(QtCore.QSize(460, 620))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.verticalWidget.setFont(font)
        self.verticalWidget.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.verticalWidget.setObjectName("verticalWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalWidget = QtWidgets.QWidget(self.verticalWidget)
        self.horizontalWidget.setMinimumSize(QtCore.QSize(440, 60))
        self.horizontalWidget.setMaximumSize(QtCore.QSize(440, 60))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.horizontalWidget.setFont(font)
        self.horizontalWidget.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.horizontalWidget.setObjectName("horizontalWidget")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout(self.horizontalWidget)
        self.horizontalLayout_10.setContentsMargins(11, -1, -1, -1)
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.open_file_button = QtWidgets.QPushButton(self.horizontalWidget)
        self.open_file_button.setMinimumSize(QtCore.QSize(150, 45))
        self.open_file_button.setMaximumSize(QtCore.QSize(150, 45))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.open_file_button.setFont(font)
        self.open_file_button.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.open_file_button.setObjectName("open_file_button")
        self.horizontalLayout_10.addWidget(self.open_file_button)
        self.paste_seq_button = QtWidgets.QPushButton(self.horizontalWidget)
        self.paste_seq_button.setMinimumSize(QtCore.QSize(150, 45))
        self.paste_seq_button.setMaximumSize(QtCore.QSize(160, 50))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.paste_seq_button.setFont(font)
        self.paste_seq_button.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.paste_seq_button.setObjectName("paste_seq_button")
        self.horizontalLayout_10.addWidget(self.paste_seq_button)
        self.verticalLayout.addWidget(self.horizontalWidget)
        self.open_widget = QtWidgets.QWidget(self.verticalWidget)
        self.open_widget.setMinimumSize(QtCore.QSize(450, 530))
        self.open_widget.setMaximumSize(QtCore.QSize(450, 530))
        self.open_widget.setObjectName("open_widget")
        self.open_layout = QtWidgets.QVBoxLayout(self.open_widget)
        self.open_layout.setContentsMargins(11, 0, 0, 0)
        self.open_layout.setObjectName("open_layout")
        self.instruction = QtWidgets.QTextBrowser(self.open_widget)
        self.instruction.setMinimumSize(QtCore.QSize(430, 530))
        self.instruction.setMaximumSize(QtCore.QSize(430, 530))
        self.instruction.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.instruction.setReadOnly(True)
        self.instruction.setObjectName("instruction")
        self.open_layout.addWidget(self.instruction)
        self.verticalLayout.addWidget(self.open_widget)
        self.horizontalLayout_3.addWidget(self.verticalWidget)
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setMinimumSize(QtCore.QSize(5, 0))
        self.line.setMaximumSize(QtCore.QSize(5, 610))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.line.setFont(font)
        self.line.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
        self.line.setLineWidth(2)
        self.line.setMidLineWidth(3)
        self.line.setFrameShape(QtWidgets.QFrame.Shape.VLine)
        self.line.setObjectName("line")
        self.horizontalLayout_3.addWidget(self.line)
        self.verticalWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.verticalWidget_2.setMaximumSize(QtCore.QSize(350, 16777215))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.verticalWidget_2.setFont(font)
        self.verticalWidget_2.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.verticalWidget_2.setObjectName("verticalWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalWidget_2)
        self.verticalLayout_2.setContentsMargins(10, 0, 10, 10)
        self.verticalLayout_2.setSpacing(80)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.parameters = QtWidgets.QWidget(self.verticalWidget_2)
        self.parameters.setMaximumSize(QtCore.QSize(400, 500))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.parameters.setFont(font)
        self.parameters.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.parameters.setObjectName("parameters")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.parameters)
        self.verticalLayout_3.setContentsMargins(10, 20, -1, 60)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.parameters_button = QtWidgets.QPushButton(self.parameters)
        self.parameters_button.setMinimumSize(QtCore.QSize(250, 50))
        self.parameters_button.setMaximumSize(QtCore.QSize(250, 50))
        self.parameters_button.setObjectName("parameters_button")
        self.verticalLayout_3.addWidget(self.parameters_button, 0, QtCore.Qt.AlignmentFlag.AlignHCenter)
        self.gridWidget = QtWidgets.QWidget(self.parameters)
        self.gridWidget.setMinimumSize(QtCore.QSize(315, 150))
        self.gridWidget.setMaximumSize(QtCore.QSize(315, 150))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.gridWidget.setFont(font)
        self.gridWidget.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.gridWidget.setObjectName("gridWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridWidget)
        self.gridLayout.setContentsMargins(11, -1, -1, 0)
        self.gridLayout.setVerticalSpacing(0)
        self.gridLayout.setObjectName("gridLayout")
        self.gc_label = QtWidgets.QLabel(self.gridWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.gc_label.setFont(font)
        self.gc_label.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.gc_label.setObjectName("gc_label")
        self.gridLayout.addWidget(self.gc_label, 2, 0, 1, 1)
        self.length_label = QtWidgets.QLabel(self.gridWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.length_label.setFont(font)
        self.length_label.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.length_label.setObjectName("length_label")
        self.gridLayout.addWidget(self.length_label, 1, 0, 1, 1)
        self.tm_label = QtWidgets.QLabel(self.gridWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.tm_label.setFont(font)
        self.tm_label.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.tm_label.setObjectName("tm_label")
        self.gridLayout.addWidget(self.tm_label, 3, 0, 1, 1)
        self.widget = QtWidgets.QWidget(self.gridWidget)
        self.widget.setObjectName("widget")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.widget)
        self.horizontalLayout_7.setSpacing(10)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.length_spin_1 = QtWidgets.QSpinBox(self.widget)
        self.length_spin_1.setMinimumSize(QtCore.QSize(50, 25))
        self.length_spin_1.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.length_spin_1.setFont(font)
        self.length_spin_1.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.length_spin_1.setFrame(False)
        self.length_spin_1.setKeyboardTracking(False)
        self.length_spin_1.setObjectName("length_spin_1")
        self.horizontalLayout_7.addWidget(self.length_spin_1)
        self.length_spin_2 = QtWidgets.QSpinBox(self.widget)
        self.length_spin_2.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.length_spin_2.setFont(font)
        self.length_spin_2.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.length_spin_2.setFrame(False)
        self.length_spin_2.setKeyboardTracking(False)
        self.length_spin_2.setObjectName("length_spin_2")
        self.horizontalLayout_7.addWidget(self.length_spin_2)
        self.gridLayout.addWidget(self.widget, 1, 1, 1, 1)
        self.widget1 = QtWidgets.QWidget(self.gridWidget)
        self.widget1.setObjectName("widget1")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout(self.widget1)
        self.horizontalLayout_9.setSpacing(10)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.tm_spin_1 = QtWidgets.QSpinBox(self.widget1)
        self.tm_spin_1.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.tm_spin_1.setFont(font)
        self.tm_spin_1.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.tm_spin_1.setFrame(False)
        self.tm_spin_1.setKeyboardTracking(False)
        self.tm_spin_1.setObjectName("tm_spin_1")
        self.horizontalLayout_9.addWidget(self.tm_spin_1)
        self.tm_spin_2 = QtWidgets.QSpinBox(self.widget1)
        self.tm_spin_2.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.tm_spin_2.setFont(font)
        self.tm_spin_2.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.tm_spin_2.setFrame(False)
        self.tm_spin_2.setKeyboardTracking(False)
        self.tm_spin_2.setObjectName("tm_spin_2")
        self.horizontalLayout_9.addWidget(self.tm_spin_2)
        self.gridLayout.addWidget(self.widget1, 3, 1, 1, 1)
        self.widget2 = QtWidgets.QWidget(self.gridWidget)
        self.widget2.setObjectName("widget2")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout(self.widget2)
        self.horizontalLayout_8.setSpacing(10)
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.gc_spin_1 = QtWidgets.QSpinBox(self.widget2)
        self.gc_spin_1.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.gc_spin_1.setFont(font)
        self.gc_spin_1.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.gc_spin_1.setFrame(False)
        self.gc_spin_1.setKeyboardTracking(False)
        self.gc_spin_1.setObjectName("gc_spin_1")
        self.horizontalLayout_8.addWidget(self.gc_spin_1)
        self.gc_spin_2 = QtWidgets.QSpinBox(self.widget2)
        self.gc_spin_2.setMaximumSize(QtCore.QSize(50, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.gc_spin_2.setFont(font)
        self.gc_spin_2.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.gc_spin_2.setFrame(False)
        self.gc_spin_2.setKeyboardTracking(False)
        self.gc_spin_2.setObjectName("gc_spin_2")
        self.horizontalLayout_8.addWidget(self.gc_spin_2)
        self.gridLayout.addWidget(self.widget2, 2, 1, 1, 1)
        self.verticalLayout_3.addWidget(self.gridWidget)
        self.widget3 = QtWidgets.QWidget(self.parameters)
        self.widget3.setMaximumSize(QtCore.QSize(330, 150))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.widget3.setFont(font)
        self.widget3.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.widget3.setObjectName("widget3")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.widget3)
        self.verticalLayout_4.setContentsMargins(11, 0, -1, 15)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.loop_primers_flag = QtWidgets.QCheckBox(self.widget3)
        self.loop_primers_flag.setMaximumSize(QtCore.QSize(330, 50))
        self.loop_primers_flag.setDisabled(True)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.loop_primers_flag.setFont(font)
        self.loop_primers_flag.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.loop_primers_flag.setObjectName("loop_primers_flag")
        self.verticalLayout_4.addWidget(self.loop_primers_flag)
        self.probe_flag = QtWidgets.QCheckBox(self.widget3)
        self.probe_flag.setMaximumSize(QtCore.QSize(330, 50))
        self.probe_flag.setDisabled(True)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.probe_flag.setFont(font)
        self.probe_flag.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.probe_flag.setObjectName("probe_flag")
        self.verticalLayout_4.addWidget(self.probe_flag)
        self.verticalLayout_3.addWidget(self.widget3)
        self.verticalLayout_2.addWidget(self.parameters)
        self.verticalWidget1 = QtWidgets.QWidget(self.verticalWidget_2)
        self.verticalWidget1.setMaximumSize(QtCore.QSize(330, 150))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.verticalWidget1.setFont(font)
        self.verticalWidget1.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.verticalWidget1.setObjectName("verticalWidget1")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.verticalWidget1)
        self.verticalLayout_5.setContentsMargins(11, -1, -1, -1)
        self.verticalLayout_5.setSpacing(15)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.loading = QtWidgets.QProgressBar(self.verticalWidget1)
        self.loading.setMinimumSize(QtCore.QSize(325, 20))
        self.loading.setMaximumSize(QtCore.QSize(325, 25))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.loading.setFont(font)
        self.loading.setStyleSheet("m")
        self.loading.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.loading.setProperty("value", 0)
        self.loading.setObjectName("loading")
        self.verticalLayout_5.addWidget(self.loading)
        self.horizontalWidget1 = QtWidgets.QWidget(self.verticalWidget1)
        self.horizontalWidget1.setMaximumSize(QtCore.QSize(325, 50))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.horizontalWidget1.setFont(font)
        self.horizontalWidget1.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.horizontalWidget1.setObjectName("horizontalWidget1")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.horizontalWidget1)
        self.horizontalLayout_4.setContentsMargins(0, 11, 10, -1)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.search_button = QtWidgets.QPushButton(self.horizontalWidget1)
        self.search_button.setMinimumSize(QtCore.QSize(130, 40))
        self.search_button.setMaximumSize(QtCore.QSize(130, 50))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.search_button.setFont(font)
        self.search_button.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.search_button.setObjectName("search_button")
        self.horizontalLayout_4.addWidget(self.search_button)
        self.results_button = QtWidgets.QPushButton(self.horizontalWidget1)
        self.results_button.setMinimumSize(QtCore.QSize(130, 40))
        self.results_button.setMaximumSize(QtCore.QSize(130, 50))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(50)
        self.results_button.setFont(font)
        self.results_button.setLocale(QtCore.QLocale(QtCore.QLocale.Language.English, QtCore.QLocale.Country.UnitedStates))
        self.results_button.setObjectName("results_button")
        self.horizontalLayout_4.addWidget(self.results_button)
        self.verticalLayout_5.addWidget(self.horizontalWidget1)
        self.count_sets = QtWidgets.QLabel(self.verticalWidget1)
        self.count_sets.setText("")
        self.count_sets.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.count_sets.setObjectName("count_sets")
        self.verticalLayout_5.addWidget(self.count_sets)
        self.verticalLayout_2.addWidget(self.verticalWidget1)
        self.horizontalLayout_3.addWidget(self.verticalWidget_2)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.open_file_button.setText(_translate("MainWindow", "Open File"))
        self.paste_seq_button.setText(_translate("MainWindow", "Paste Sequence"))
        self.instruction.setHtml(_translate("MainWindow", 
        """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Document</title>
        </head>
        <body>
            <p>1) You can paste nucleotide sequence from the file (select   &quot;Open File&quot; button) or paste the sequence as a text (select   &quot;Paste Sequence&quot; button).<br> 2) Set required parameters (length of the primers, %GC and   Tm).<br>3) Select &quot;Search&quot; button and wait.<br>4) When the results are ready, you can copy or save them  as *.xls file.</p>
            <p>- Github: https://github.com/Restily/LAMPrimers-iQ</p>
        </body>
        </html>
        """
        ))
        self.parameters_button.setText(_translate("MainWindow", "Primer Design Parameters"))
        self.gc_label.setText(_translate("MainWindow", "GC count (%)"))
        self.length_label.setText(_translate("MainWindow", "Length of primers (nucls)"))
        self.tm_label.setText(_translate("MainWindow", "Melting Temperature (⁰C)"))
        self.loop_primers_flag.setText(_translate("MainWindow", "Search loop primers"))
        self.probe_flag.setText(_translate("MainWindow", "Search hybridization probe"))
        self.search_button.setText(_translate("MainWindow", "Search"))
        self.results_button.setText(_translate("MainWindow", "Open results"))
