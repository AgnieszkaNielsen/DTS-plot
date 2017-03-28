# -*- coding: utf-8 -*-

import PyQt5.QtWidgets as qt
import sys
import numpy as np

import dts_functions as dfun


language = 'en'
commas = 0

# ## define labels in different languages
langs = dict()
# 0, 1, 2: used in file saving
# 3: label on a colorbar
# 4, 5: x and y axis
# 6, 7: plot title
langs['dk'] = ["gemt_den", "komma", "alle", "temperatur", "kabellængde (m)",
               "tid (hh:mm)", "start", "stop"]
langs['en'] = ["saved_on", "comma", "all", "temperature", "cabel length (m)",
               "time (hh:mm)", "start", "stop"]


class NirasGui(qt.QWidget):

    def __init__(self):

        super().__init__()

        grid = qt.QGridLayout()
        self.setLayout(grid)

        self.button1 = qt.QPushButton('import')
        self.button2 = qt.QPushButton('parameters')
        self.info = qt.QLineEdit()
        self.info.setReadOnly(True)

        grid.addWidget(self.button1, 0, 0, 2, 2)
        grid.addWidget(self.button2, 0, 4, 2, 2)
        grid.addWidget(self.info, 4, 0, 10, 10)

        self.button1.clicked.connect(self.importfiles)
        self.button2.clicked.connect(self.addone)

        self.resize(250, 150)
        self.move(300, 300)
        self.setWindowTitle('Simple')
        self.show()

        self.count = 0

    def printaes(self):
        print(self.count)

    def addone(self):
        self.count += 1

    def importfiles(self):
        filenames = qt.QFileDialog.getOpenFileNames(self, 'Open file',
                                                    '', '(*.tra *.txt)')

        trafiles = [f for f in sorted(filenames[0]) if '.tra' in f]
        txtfiles = [f for f in sorted(filenames[0]) if '.txt' in f]

        if len(trafiles) > 0:
            self.info.setText('import .tra file(s)')

            print("read data from files")
            # initialise the files
            (eltime, timelabels, tempstack, tempstack_comma, lines,
             dataheader, savedate,
             initialtime) = dfun.get_data_initfile(trafiles[0])

            for trafile in trafiles[1:]:
                (eltime, timelabels, tempstack, tempstack_comma,
                 dataheader) = dfun.get_data(trafile, dataheader, eltime,
                                             timelabels, initialtime,
                                             tempstack, tempstack_comma)

            print("save the results (text file)")

            idno = np.arange(0, lines+1)

            saveddata, savenamefig = dfun.savefile(lines, tempstack,
                                                   tempstack_comma,
                                                   commas, langs, language,
                                                   savedate,
                                                   dataheader, idno)

        elif len(txtfiles) > 1:
            self.info.setText('select only ONE .txt file')

        else:
            self.info.setText('import .txt file')


if __name__ == '__main__':
    app = qt.QApplication(sys.argv)
    gui = NirasGui()
    sys.exit(app.exec_())
