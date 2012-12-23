
import QtQuick 1.0

Row {
    id: mc
    property int frame: 0
    property int first_frame: 0
    property int last_frame: 10
    property real frame_rate: 5 
    
    onFrameChanged: frame_indicator.text = "at frame " + frame
    PropertyAnimation {
        id: play_frames
        target: mc
        property: "frame"
        to: last_frame
        duration: (last_frame - from) / frame_rate * 1000 // [ms]
    }

    MiniButton {
        id: back_full
        image: "img/back_full.jpg"
        onClicked: frame = first_frame
    }
    MiniButton {
        id: back_one
        image: "img/back_one.jpg"
        onClicked: if (frame > first_frame) { frame -= 1 }
    }
    MiniButton {
        id: play_pause
        image: "img/play.jpg"
        onClicked: {
            if (!play_frames.running) { 
                play_frames.from = frame
                play_frames.start() 
            } else {
                play_frames.stop()
            }
        }
    }
    MiniButton {
        id: forward_one
        image: "img/forward_one.jpg"
        onClicked: if (frame < last_frame) { frame += 1 }
    }
    MiniButton {
        id: forward_full
        image: "img/forward_full.jpg"
        onClicked: frame = last_frame
    }

    Rectangle {
        id: spacer
        width: 5
        height: 1
        color: "transparent"
    }
    Text {
        id: frame_indicator
        text: "at frame 0   "
    }
}
