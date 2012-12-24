
import QtQuick 1.0

Item {
    width: cams.width
    height: cams.height + controls.height
    
    property string scene_template: "../pyptv_gui/test/testing_fodder/clean_scene/cam%1_Scene83_%2"

    function pad(n, width, z) {
        z = z || '0';
        n = n + '';
        return n.length >= width ? n : new Array(width - n.length + 1).join(z) + n;
    }

    function load_frame() {
        for (var cam = 0; cam < 4; cam++) {
            cams.images.setProperty(cam, "file",
                scene_template.arg(cam + 1).arg(pad(controls.frame, 4, '0')))
        }
    }
    
    Multicam {
        id: cams
    }
    MovieControl {
        id: controls
        anchors.left: cams.left
        anchors.right: cams.right
        anchors.top: cams.bottom
        
        first_frame: 497
        last_frame: 597
        frame_rate: 100
    }
    Component.onCompleted: {
        controls.frameChanged.connect(load_frame)
    }
}

