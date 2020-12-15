window.io = require("socket.io-client")
const RTCMultiConnection = require('rtcmulticonnection')

class RemoteStream {

    constructor(audioContext) {
        this.audioContext = audioContext
        this.sink = audioContext.createMediaStreamDestination();
        this.streams = {}
        window.streams = this.streams
    }

    // connecting to destination doesnt work. Currently audio elements are on by default
    connect(sourceNode, destinationNode, roomID="squidback") {
        this.destination = destinationNode
        this.connection = new RTCMultiConnection();
        // by default, socket.io server is assumed to be deployed on your own URL
        // this.connection.socketURL = '/';

        // comment-out below line if you do not have your own socket.io server
        this.connection.socketURL = 'https://rtcmulticonnection.herokuapp.com:443/';
        this.connection.socketMessageEvent = 'audio-conference-demo';
        this.connection.session = { audio: true, video: false };
        this.connection.mediaConstraints = { audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }, video: false };
        this.connection.sdpConstraints.mandatory = { OfferToReceiveAudio: true, OfferToReceiveVideo: false };
        // https://www.rtcmulticonnection.org/docs/iceServers/
        // use your own TURN-server here!
        /*this.connection.iceServers = [{
            'urls': [
                'stun:stun.l.google.com:19302',
                'stun:stun1.l.google.com:19302',
                'stun:stun2.l.google.com:19302',
                'stun:stun.l.google.com:19302?transport=udp',
            ]
        }];*/
        this.connection.dontCaptureUserMedia = true
        this.connection.autoCreateMediaElement = false
        sourceNode.connect(this.sink)
        this.connection.multiPeersHandler.onGettingLocalMedia(this.sink.stream)

        this.connection.onstream = (event) => {
            this.connectStream(event.stream, event.mediaElement)
        };

        this.connection.onstreamended = (event) => {
            this.disconnectStream(event.stream.streamid);
        };

        this.connection.openOrJoin(roomID)

    }

    connectStream(stream, mediaElement) {
        if(stream.type == 'local') return
        const {streamid} = stream;
        mediaElement = document.createElement("audio")
        mediaElement.srcObject = stream;
        mediaElement.muted = true

        this.streamSource = this.audioContext.createMediaStreamSource(stream);
        const gain = this.audioContext.createGain();
        gain.gain.value = 1
        this.streamSource.connect(gain).connect(this.destination);
        this.streams[streamid] = gain;
        console.log("[webrtc] connecting", stream, mediaElement)
        document.body.appendChild(mediaElement)
        this.adjustVolumes();
    }

    disconnectStream(streamid) {
        /*if(stream.type == 'local') return
        console.log("[webrtc] disconnecting", streamid)
        if(this.streams[streamid]) {
            //this.streams[streamid].disconnect();
            delete this.streams[streamid];
        }

        var mediaElement = document.getElementById(streamid);
        if (mediaElement) {
            mediaElement.parentNode.removeChild(mediaElement);
        }*/
    }

    adjustVolumes() {
        const numStreams = Object.values(this.streams).length + 1;
        const newVolume = 0.1 / Math.sqrt(numStreams);

        for(const streamID in this.streams) {
            console.log("[webrtc] mixing", streamID, newVolume)
            this.streams[streamID].gain.value = streamID == this.streamid ? 0 : newVolume
        }
    }
}

module.exports = RemoteStream