# A command-line tool for viewing four-camera movie directories as animation.
# Mostly loads the appropriate QML component and sets its properties from
# command-line arguments.

import sys
from PySide.QtCore import QUrl
from PySide.QtGui import QApplication
from PySide.QtDeclarative import QDeclarativeView

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("template", help="Temlate for image file name, containing "
    "%%1 where camera number should be inserted, and %%2 where frame number "
    "should be inserted. Example: my_scene/cam%%1.%%2")
parser.add_argument("first", type=int, help="Number of first frame")
parser.add_argument("last", type=int, help="Number of last frame")
parser.add_argument("--rate", type=int, default=100,
    help="Frame rate in frames per second")
args = parser.parse_args()

app = QApplication(sys.argv)
view = QDeclarativeView()
view.setSource(QUrl('MovieBrowser.qml'))

view.rootObject().setProperty("scene_template", args.template)
view.rootObject().setProperty("first_frame", args.first)
view.rootObject().setProperty("last_frame", args.last)
view.rootObject().setProperty("frame_rate", args.rate)

view.show()
sys.exit(app.exec_())

