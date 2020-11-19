const SquidbackGraph = require('./graph.js')
const { chromaticFilter } = require('./scales.js')
const Graph = new SquidbackGraph();

class SquidbackWorkletProcess {

    constructor(audioContext, fftSize = 1024) {
        this.audioContext = audioContext;
        this.fftSize = fftSize;
        this.gainIncrement = 0.001;
        this.gainDecrement = 0.001;
        this.minVolume = -10;
        this.maxVolume = -10;

        this.initBuffers();
    }

    async start() {
        await this.audioContext.audioWorklet.addModule('phase-vocoder.js');
        this.initExtremesFilters();
        this.initFFTNode();
        this.initGain();
        //this.initFilter();
        await this.initInput();
        if(this.input) this.connectAll();
    }

    initFilter() {
        /*const scale = chromaticFilter(this.binToFreq, 22000, 12)
        this.filter = new NotchFilterBank(this.audioContext, scale.scale, 100);
        this.freqResp = new Float32Array(this.fftSize/2+1);
        this.freqResp.fill(1);
        this.filter.getFrequencyResponse(freqs, this.freqResp);*/

        const freqs = new Float32Array(this.fftSize/2+1);
        freqs.forEach((f,n)=>freqs[n]=n*this.binToFreq);
        // console.log(scale.scale)
    }

    initGain() {
        this.mainGain = this.audioContext.createGain();
        this.limiter = this.audioContext.createDynamicsCompressor()
        this.limiter.threshold.value = -40
        this.limiter.ratio.value = 20
        this.limiter.attack.value = 0.001
        this.limiter.release.value = 0.01
    }

    initBuffers() {
        this.binToFreq = this.audioContext.sampleRate / this.fftSize;
        const numBins = this.fftSize/2 + 1;
        this.fftBuffer = new Float32Array(numBins);
        this.peaksBuffer = new Int32Array(numBins);
        this.numPeaks = 0;
        this.fbBuffer = new Int32Array(numBins);
        this.numFb = 0;
        this.slopeBuffer = new Float32Array(numBins);
        this.filterBuffer = new Int32Array(numBins);

        this.maxDb = -180;
    }

    now() {
        return this.audioContext.currentTime
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter();
        hipass.type = "highpass"; hipass.frequency.setValueAtTime(this.binToFreq, this.now());
        hipass.Q.setValueAtTime(0.707, this.now());

        const lopass = this.audioContext.createBiquadFilter();
        lopass.type = "lowpass"; lopass.frequency.setValueAtTime(18000, this.now());

        this.hipass = hipass; this.lopass = lopass;
    }

    initFFTNode() {
        this.fftNode = new AudioWorkletNode(this.audioContext, 'phase-vocoder-processor', {processorOptions: {blockSize: this.fftSize}});
        this.initWorkletPort();
    }

    initWorkletPort() {
        this.fftNode.port.onmessage = ({data}) => {
            switch(data[0]) {
                case "spectrum": this.fftBuffer.set(data[1]); break;
                case "slope": this.slopeBuffer.set(data[1]); break;
                case "histo": this.histoBuffer.set(data[1]); break;
                case "peaks":
                    this.numPeaks = data[1]
                    this.peaksBuffer.set(data[2])
                    break;
                case "fb":
                    this.numFb = data[1]
                    this.fbBuffer.set(data[2])
                    this.onFbData();
                    break;
                case 'maxDb': this.maxDb = data[1]; break;
                case 'reductions': this.filterBuffer = data[1]; break;
            }
            //const sum = Math.max(...data);
            // if(data[0] > 0) console.log("Fb candidates:", data[0], data[1])
            // console.log(data)
        };
        requestAnimationFrame(()=>this.showSpectrum());
        requestAnimationFrame(()=>this.showSlope());
    }

    connectAll() {
        this.input
        .connect(this.hipass)
        .connect(this.lopass)
        .connect(this.fftNode)

        //this.filter.connectInput(this.lopass)
        //this.filter
        this.lopass.connect(this.mainGain).connect(this.limiter).connect(this.audioContext.destination)

        //this.lopass.connect(this.mainGain).connect(this.audioContext.destination)
    }

    async initInput() {
        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await navigator.mediaDevices.getUserMedia(constraints);
           this.input = this.audioContext.createMediaStreamSource(stream);
           console.log("GOT INPUT")
         } catch(err) {
           /* handle the error */
            console.error("CANT GET USER INPUT")
         }
    }
    

    onFbData() {
        this.updateGain()
    }

    updateGain() {
        if(this.maxDb <= -180) return
        const gain = this.mainGain.gain;
        gain.cancelAndHoldAtTime(this.now())
        const currentGain = 20 * Math.log10(this.mainGain.gain.value);
        let nextGain = currentGain;
        if(this.maxDb > this.maxVolume) {
            nextGain += (this.minVolume - this.maxDb) * this.gainDecrement
            //console.log("gain--", this.maxDb, nextGain)
        } else {
            nextGain += (this.maxVolume - this.maxDb) * this.gainIncrement
            //console.log("gain++",this.maxDb, nextGain)
        }
        //console.log(this.maxDb,nextGain)
        if(nextGain != currentGain) {
            nextGain = Math.pow(10, nextGain*0.05);
            if(isFinite(nextGain) && nextGain > 0)
                gain.linearRampToValueAtTime(nextGain, this.now())
        }

    }

    // graphs
    showSpectrum(){
        requestAnimationFrame(()=>this.showSpectrum())
        this.fftNode.port.postMessage("getSpectrum");
        this.fftNode.port.postMessage("getPeaks");
        this.fftNode.port.postMessage("getFbCandidates");
        this.fftNode.port.postMessage("getMaxDb");
        this.fftNode.port.postMessage("getMagReductions");

        const bgColor = this.numFb > 0 ? Graph.pitchAmpToHSLA(this.fbBuffer[0] * this.binToFreq, this.fftBuffer[this.fbBuffer[0]]) : 'rgba(100,100,100,1)';


        //Graph.drawBg("canvas#spectrum", 'rgba(0,0,0,0.5)')
        Graph.drawBg("canvas#spectrum", bgColor)
        Graph.drawGain("canvas#spectrum", "red", this.mainGain.gain.value, -100, 80)

        Graph.drawLogLine("canvas#spectrum", this.fftBuffer, -100, 0);
        //Graph.drawLogLine("canvas#spectrum", this.freqResp, 0, 1, "red");

        Graph.drawLogLineInverted("canvas#spectrum", this.filterBuffer, -100, 0, "rgba(0,0,0,0.5)");

        //Graph.drawVerticals("canvas#spectrum", this.numPeaks, this.peaksBuffer, -80, 0, "rgba(255,255,255,1)")
        //Graph.drawVerticals("canvas#spectrum", this.numFb, this.fbBuffer, -80, 0, "rgba(255,0,0,0.5)")
    }

    showSlope(node){
        requestAnimationFrame(()=>this.showSlope())
        this.fftNode.port.postMessage("getSlope");
        //console.log(slopeBuffer)
        Graph.drawBg("canvas#slope", 'white')
        Graph.drawLogLine("canvas#slope", this.slopeBuffer, 0, 500, "black");
    }
}

module.exports = SquidbackWorkletProcess