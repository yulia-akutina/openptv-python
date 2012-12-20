
import QtQuick 1.0

Row {
    MiniButton {
        id: back_full
        image: "img/back_full.jpg"
    }
    MiniButton {
        id: back_one
        image: "img/back_one.jpg"
    }
    MiniButton {
        id: play_pause
        image: "img/play.jpg"
    }
    MiniButton {
        id: forward_one
        image: "img/forward_one.jpg"
    }
    MiniButton {
        id: forward_full
        image: "img/forward_full.jpg"
    }

    Rectangle {
        id: spacer
        width: 5
        height: 1
        color: "transparent"
    }
    Text {
        id: frame_scrollbar
        text: "at frame 0"
    }
}
