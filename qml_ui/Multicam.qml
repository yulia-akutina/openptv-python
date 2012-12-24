
import QtQuick 1.0

Item {
    id: cams_view
    width: 600; height: 600
    
    property alias images: init_images
    ListModel {
        id: init_images
        ListElement { file: "" }
        ListElement { file: "" }
        ListElement { file: "" }
        ListElement { file: "" }
    }
    
    Component {
        id: cam_image
        Column {
            Text { 
                width: cam_grid.width / 2
                text: "Camera " + (index + 1) + " " + file }
            Image { 
                source: file
                width: cam_grid.width / 2; 
                height: cam_grid.height / 2
            }
        }
    }
    
    Grid {
        id: cam_grid
        rows: 2; columns: 2
        flow: Grid.LeftToRight
        anchors.fill: parent
        
        Repeater {
            model: init_images
            delegate: cam_image
        }
    }
}

