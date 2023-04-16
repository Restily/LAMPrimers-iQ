import os

from Bio import SeqRecord, Seq
from PyQt6 import QtWidgets, QtGui, QtCore

from gui.ui_py.params import Ui_Params
from gui.ui_py.ui import Ui_MainWindow
from gui.ui_py.results import Ui_Results

from lamp.config import DesignConfig, LAMPConfig
from lamp.start_lamp import LAMP

from gui.utils import *


# class Loading(QtCore.QObject):



class ParamsWindow(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        self.params_ui = Ui_Params()
        self.params_ui.setupUi(self)

        self.init_UI()


    def init_UI(self):
        self.setFixedSize(350, 500)

        self.setWindowFlags(QtCore.Qt.WindowType.WindowStaysOnTopHint | 
                            QtCore.Qt.WindowType.CustomizeWindowHint)
        
        self.connect_buttons()


    def connect_buttons(self):
        self.params_ui.save_changes.clicked.connect(self.save_changes)
        self.params_ui.cancel_changes.clicked.connect(self.cancel_changes)


    def open_window(self):
        self.set_config_params()

        self.show()


    def set_config_params(self):
        self.params_ui.f3_f2_spin_1.setValue(
            DesignConfig.F3_F2[0]
        )
        self.params_ui.f3_f2_spin_2.setValue(
            DesignConfig.F3_F2[1]
        )

        self.params_ui.f2_f1c_spin_1.setValue(
            DesignConfig.F2_F1C[0]
        )
        self.params_ui.f2_f1c_spin_2.setValue(
            DesignConfig.F2_F1C[1]
        )

        self.params_ui.f1c_b1c_spin_1.setValue(
            DesignConfig.F1C_B1C[0]
        )
        self.params_ui.f1c_b1c_spin_2.setValue(
            DesignConfig.F1C_B1C[1]
        )

        self.params_ui.b2_b1c_spin_1.setValue(
            DesignConfig.B2_B1C[0]
        )
        self.params_ui.b2_b1c_spin_2.setValue(
            DesignConfig.B2_B1C[1]
        )

        self.params_ui.b3_b2_spin_1.setValue(
            DesignConfig.B3_B2[0]
        )
        self.params_ui.b3_b2_spin_2.setValue(
            DesignConfig.B3_B2[1]
        )

        self.params_ui.amplicon_spin_1.setMaximum(500)
        self.params_ui.amplicon_spin_2.setMaximum(500)

        self.params_ui.amplicon_spin_1.setValue(
            DesignConfig.AMPLICON_RANGE[0]
        )
        self.params_ui.amplicon_spin_2.setValue(
            DesignConfig.AMPLICON_RANGE[1]
        )
        self.params_ui.delta_length_spin.setValue(
            DesignConfig.MAX_DIFFERENCE_LENGTHS_PRIMERS
        )

        self.params_ui.na_plus_spin.setValue(
            LAMPConfig.NA_PLUS
        )
        self.params_ui.delta_tm_spin.setValue(
            LAMPConfig.TM_MAX_DIFFERENCE
        )


    def save_changes(self):
        DesignConfig.F3_F2 = [
            self.params_ui.f3_f2_spin_1.value(),
            self.params_ui.f3_f2_spin_2.value()
        ]

        DesignConfig.F2_F1C = [
            self.params_ui.f2_f1c_spin_1.value(),
            self.params_ui.f2_f1c_spin_2.value()
        ]

        DesignConfig.F1C_B1C = [
            self.params_ui.f1c_b1c_spin_1.value(),
            self.params_ui.f1c_b1c_spin_2.value()
        ]

        DesignConfig.B2_B1C = [
            self.params_ui.b2_b1c_spin_1.value(),
            self.params_ui.b2_b1c_spin_2.value()
        ]

        DesignConfig.B3_B2 = [
            self.params_ui.b3_b2_spin_1.value(),
            self.params_ui.b3_b2_spin_2.value()
        ]

        DesignConfig.AMPLICON_RANGE[0] = self.params_ui.amplicon_spin_1.value()
        DesignConfig.AMPLICON_RANGE[1] = self.params_ui.amplicon_spin_2.value()
        DesignConfig.MAX_DIFFERENCE_LENGTHS_PRIMERS = self.params_ui.delta_length_spin.value()

        LAMPConfig.NA_PLUS = self.params_ui.na_plus_spin.value()
        LAMPConfig.TM_MAX_DIFFERENCE = self.params_ui.delta_tm_spin.value()

        self.close()


    def cancel_changes(self):
        self.close()        


class LoadingWindow(QtWidgets.QWidget):
    def __init__(self):
        """
        !!!
        """
        super().__init__()

        self.init_UI()


    def init_UI(self):
        self.setFixedSize(300, 300)
        self.setWindowFlags(QtCore.Qt.WindowType.WindowStaysOnTopHint | 
                            QtCore.Qt.WindowType.CustomizeWindowHint)

        self.label_animation = QtWidgets.QLabel(self)

        self.loading = QtGui.QMovie('./gui/media/loading2.gif')
        self.label_animation.setMovie(self.loading)


    def startAnimation(self):
        self.show()

        self.loading.start()


    def stopAnimation(self):
        self.loading.stop()

        self.close()


class ResultsWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.results_ui = Ui_Results()

        self.results_ui.setupUi(self)
        self.init_UI()

    
    def init_UI(self):
        self.connect_buttons()


    def connect_buttons(self):
        self.results_ui.next_set.clicked.connect(self.next_result)
        self.results_ui.previous_set.clicked.connect(self.previous_result)
        self.results_ui.save_excel.clicked.connect(self.save_file_to_excel)
    

    def open_results(self):
        self.show()

        self.set_results()

    
    def set_results(self):
        if Results.cur_primer_set_ind == 0:
            self.results_ui.previous_set.setDisabled(True)
            
        primer_set = Results.primer_sets[Results.cur_primer_set_ind]

        sequences_str = get_sequences(primer_set)

        self.results_ui.primers_set.setHtml(primer_sets_table(primer_set))
        
        self.results_ui.sequences.setLineWrapMode(QtWidgets.QTextEdit.LineWrapMode.NoWrap)

        self.results_ui.sequences.setHtml(sequences_str)
                    
        if len(Results.primer_sets) == 1:
            self.results_ui.next_set.setDisabled(True)


    def next_result(self):
        self.results_ui.previous_set.setDisabled(False)

        Results.cur_primer_set_ind += 1

        self.set_results()

        if Results.cur_primer_set_ind == len(Results.primer_sets) - 1:
            self.results_ui.next_set.setDisabled(True)
    

    def previous_result(self):
        self.results_ui.next_set.setDisabled(False)

        Results.cur_primer_set_ind -= 1

        self.set_results()


    def save_file_to_excel(self):
        """
        Save primer set to excel
        """
        if Results.current_record.id:
            default_filename = f'{Results.current_record.id}_primers.xlsx'
        else:
            default_filename = 'primers.xlsx'

        save_file = QtWidgets.QFileDialog.getSaveFileName(
            parent=self,
            caption='Save primers to Excel',
            directory=default_filename,
            filter='Excel File (*.xlsx);'
        )

        excel_path = save_file[0]
        
        flag = save_to_excel(
            [Results.primer_sets[Results.cur_primer_set_ind]],
            excel_path
        )


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        self.ui = Ui_MainWindow()
        self.results_window = ResultsWindow()
        self.loading_window = LoadingWindow()
        self.params_window = ParamsWindow()

        self.ui.setupUi(self)
        self.init_UI()


    def init_UI(self):
        self.setWindowTitle('Design primers for LAMP') 
        self.setWindowIcon(QtGui.QIcon('exchanging.png'))

        self.ui.results_button.setDisabled(True)
        self.ui.loading.setVisible(False)

        # Init UI
        self.connect_buttons()
        self.set_flags()
        self.set_values()


    def set_flags(self):
        self.paste_sequence_open = False
        self.open_layout_open = False
        self.selected_record = None
        self.changed = False


    def connect_buttons(self):
        # "Paste Sequence" button
        self.ui.paste_seq_button.clicked.connect(self.paste_sequence)
        
        # "Open File" button
        self.ui.open_file_button.clicked.connect(self.open_file)

        # "Search" button
        self.ui.search_button.clicked.connect(self.start_design_primers)

        # "Open results" button
        self.ui.results_button.clicked.connect(self.results_window.open_results)

        # "Parameters" button
        self.ui.parameters_button.clicked.connect(self.params_window.open_window)


    def set_values(self):
        """
        Set values to spin boxes
        """
        # GC values
        self.ui.gc_spin_1.setValue(40)
        self.ui.gc_spin_2.setValue(60)

        # Tm values
        self.ui.tm_spin_1.setValue(55)
        self.ui.tm_spin_2.setValue(65)

        # Length values
        self.ui.length_spin_1.setValue(20)
        self.ui.length_spin_2.setValue(25)


    def set_config_values(self):
        """
        Set values for design primers
        """
        # GC values
        LAMPConfig.GC_RANGE = [
            self.ui.gc_spin_1.value(),
            self.ui.gc_spin_2.value()
        ]

        # Tm values
        LAMPConfig.TM_RANGE = [
            self.ui.tm_spin_1.value(),
            self.ui.tm_spin_2.value()
        ]

        # Length values
        LAMPConfig.PRIMERS_LENGTH_RANGE = [
            self.ui.length_spin_1.value(),
            self.ui.length_spin_2.value()
        ]


    def start_design_primers(self):
        """
        :button: search_button (Search)

        Start design primers
        """
        self.set_config_values()

        self.ui.results_button.setDisabled(True)
        
        # Check sequence in file or was pasted
        if self.open_layout_open:
            self.search_record = self.selected_record
        elif self.paste_sequence_open:
            self.search_record = SeqRecord.SeqRecord(
                seq=Seq.Seq(self.paste_field.toPlainText()),
                id='LAMPrimers-IQ'
            )

        # self.setup_loading_bar()

        # Create LAMP class
        lamp = LAMP()

        # Start design primers in LAMP class
        self.primer_sets = lamp.start_design_primers(self.search_record)

        # Checking if primers have been found
        if len(self.primer_sets) == 0:
            self.ui.count_sets.setText('No primer sets was found')
        else:
            self.ui.count_sets.setText(f'{len(self.primer_sets)} primer sets was found')
            self.ui.results_button.setDisabled(False)

            # Write in Results class sets
            Results.primer_sets = self.primer_sets
            Results.current_record = self.search_record


    def select_chromosome(self, item):
        """
        Search and write selected chromosome
        """
        chromosome_id = item.text()

        for record in self.records:
            if record.id == chromosome_id:
                self.selected_record = record

                self.chromosome_label.setText(get_record_params(record))


    def open_file(self):
        # self.filename = QtWidgets.QTextBrowser()

        # self.file_layout.addWidget(self.filename)
        # self.file_layout.addWidget(QtWidgets.QPushButton("Open File"))

        # self.ui.open_layout.addLayout(self.file_layout)
        
        # Filter for file extensions
        filter = 'Sequence File (*.fas *.fasta *.fna *.ffn *.faa *.frn *.afa *.mfa);;'

        # Open dialog window
        open_file = QtWidgets.QFileDialog.getOpenFileName(
            parent=self,
            caption='Select a data file',
            directory=os.getcwd(),
            filter=filter
        )

        file_path = open_file[0]

        # Check if file was opened
        if file_path == '':
            return
        
        # Start loading
        self.loading_window.startAnimation()

        # If file was opened delete all widgets
        if self.open_layout_open or self.paste_sequence_open:
            self.selected_record = None

            for i in reversed(range(self.ui.open_layout.count())): 
                self.ui.open_layout.itemAt(i).widget().setParent(None)

        self.open_layout_open = True
        
        # Hide instruction window
        self.ui.instruction.hide()

        # Parse file
        self.records = parse_sequence_file(file_path)

        seq_name = get_filename(file_path)

        self.seq_name = QtWidgets.QLabel()
        self.seq_name.setText(f'Sequence File: {seq_name}')

        self.sequneces = QtWidgets.QListWidget()

        self.sequneces.addItems(get_record_ids(self.records))

        self.sequneces.itemClicked.connect(self.select_chromosome)
        
        self.chromosome_label = QtWidgets.QTextBrowser()
        self.chromosome_label.setMinimumSize(QtCore.QSize(440, 120))
        self.chromosome_label.setMaximumSize(QtCore.QSize(440, 120))

        self.ui.open_layout.addWidget(self.seq_name)
        self.ui.open_layout.addWidget(self.sequneces)
        self.ui.open_layout.addWidget(self.chromosome_label)
        
        self.loading_window.stopAnimation()


    def check_sequence_length(self):
        if not self.changed:
            seq = self.paste_field.toPlainText()
        
            self.changed = True

            seq = seq.replace('\r', '').replace('\n', '').replace(
                ' ', '').replace('\t', '').upper()
            seq = ''.join(filter(lambda x: not x.isdigit(), seq))
            self.paste_field.setPlainText(seq)

            # self.paste_field.moveCursor(QTextCursor.MoveOperation.End)

            self.seq_length.setText(f"Parsed Sequence (5'→3')\t{len(seq)} bp")
        
        self.changed = False


    def paste_sequence(self):
        if self.open_layout_open or self.paste_sequence_open:
            self.selected_record = None

            for i in reversed(range(self.ui.open_layout.count())): 
                self.ui.open_layout.itemAt(i).widget().setParent(None)
        
        self.paste_sequence_open = True

        self.ui.instruction.hide()
        
        font = QtGui.QFont()
        font.setPointSize(10)

        self.paste_field = QtWidgets.QPlainTextEdit()
        self.paste_field.setFont(font)
        self.paste_field.setMinimumSize(QtCore.QSize(440, 500))
        self.paste_field.setMaximumSize(QtCore.QSize(440, 500))

        self.seq_length = QtWidgets.QLabel()
        self.seq_length.setText("Parsed Sequence (5'→3')\t0 bp")
        self.seq_length.setMinimumSize(QtCore.QSize(440, 20))
        self.seq_length.setMaximumSize(QtCore.QSize(440, 20))

        self.ui.open_layout.addWidget(self.seq_length)
        self.ui.open_layout.addWidget(self.paste_field)

        self.paste_field.textChanged.connect(self.check_sequence_length)

