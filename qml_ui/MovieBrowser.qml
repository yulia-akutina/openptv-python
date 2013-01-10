
import QtQuick 1.0

Item {
    width: cams.width
    height: cams.height + controls.height
    
    property string scene_template
    property alias first_frame: controls.first_frame
    property alias last_frame: controls.last_frame
    property alias frame_rate: controls.frame_rate

    function pad(n, width, z) {
        z = z || '0';
        n = n + '';
        return n.length >= width ? n : new Array(width - n.length + 1).join(z) + n;
    }

    function load_frame() {
        for (var cam = 0; cam < 4; cam++) {
            cams.images.setProperty(cam, "file",
                scene_template.arg(cam + 1).arg(pad(controls.frame, 4, '0')));
            console.log(image_explorer.image_targets(controls.frame, cam))
            cams.images.setProperty(cam, "target_list",
                image_explorer.image_targets(controls.frame, cam));
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
    }
    Component.onCompleted: {
        controls.frameChanged.connect(load_frame)
    }
}

