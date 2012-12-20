
import QtQuick 1.0

Rectangle {
    property alias image: img.source
    width: 16
    height: 16
    color: "lightgrey"
    border.width: 1
    
    Image {
        x: 1; y: 1
        id: img
    }
}

