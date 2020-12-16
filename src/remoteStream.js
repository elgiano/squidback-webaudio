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
    connect(sourceNode, destinationNode, roomID="squidback", callback) {

        this.connection = new RTCMultiConnection();     
        this.connection.socketURL = 'https://squidback.xyz:443/';
        this.connection.session = { audio: true, video: false };
        this.connection.mediaConstraints = { audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }, video: false };
        this.connection.sdpConstraints.mandatory = { OfferToReceiveAudio: true, OfferToReceiveVideo: false };

        this.connection.dontCaptureUserMedia = true
        this.connection.autoCreateMediaElement = false
        sourceNode.connect(this.sink)
        this.connection.multiPeersHandler.onGettingLocalMedia(this.sink.stream)
        this.destination = destinationNode

        this.connection.onstream = (event) => {
            this.connectStream(event.stream, event.mediaElement)
            callback()
        };

        this.connection.onstreamended = (event) => {
            this.disconnectStream(event.stream.streamid);
            callback()
        };

        this.connection.openOrJoin(roomID)
    }

    connectStream(stream, mediaElement) {
        if(stream.type == 'local') return
        const {streamid} = stream;
        mediaElement = document.createElement("audio")
        mediaElement.id = streamid
        mediaElement.srcObject = stream;
        mediaElement.muted = true

        this.streamSource = this.audioContext.createMediaStreamSource(stream);
        const gain = this.audioContext.createGain();
        gain.gain.value = 1
        this.streamSource.connect(gain).connect(this.destination);
        this.streams[streamid] = gain;
        console.log("[webrtc] connecting", this.connection, stream, mediaElement)
        document.body.appendChild(mediaElement)
        this.adjustVolumes();
    }

    disconnectStream(streamid) {
        if(stream.type == 'local') return
        console.log("[webrtc] disconnecting", streamid)
        if(this.streams[streamid]) {
            this.streams[streamid].disconnect();
            delete this.streams[streamid];
        }

        var mediaElement = document.getElementById(streamid);
        if (mediaElement) {
            mediaElement.parentNode.removeChild(mediaElement);
        }
    }

    adjustVolumes() {
        const numStreams = Object.values(this.streams).length + 1;
        const newVolume = 0.5 / Math.sqrt(numStreams);

        for(const streamID in this.streams) {
            console.log("[webrtc] mixing", streamID, newVolume)
            this.streams[streamID].gain.value = streamID == this.streamid ? 0 : newVolume
        }
    }

    get peerCount() { return Object.values(this.streams).length }
    get isConnected() { return this.peerCount > 0 }
}

module.exports = RemoteStream
