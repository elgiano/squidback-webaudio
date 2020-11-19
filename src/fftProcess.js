const SquidbackGraph = require('./graph.js')
const Graph = new SquidbackGraph();
const { 
    AutoGain, 
    FFTConvFilter,
    MagnitudesHistory
} = require("./audio-processors/index.js")

class SquidbackFFTProcess {

    constructor(audioContext, fftSize = 1024) {
        this.audioContext = audioContext;
    }

    async start() {
        this.initFFTNode();
        this.initBuffers();
        this.initConv();
        this.initExtremesFilters();
        this.initGain();
        await this.initInput();
        if(this.input) this.connectAll();
        requestAnimationFrame(()=>this.update());
    }

    initGain() {
        this.autoGain = new AutoGain(this.audioContext)
    }

    initBuffers() {
        const numBins = this.fftNode.frequencyBinCount;
        this.numBins = numBins;
        this.binToFreq = this.audioContext.sampleRate / 2 / numBins;
        this.fftBuffer = new Float32Array(numBins);
        this.outFftBuffer = new Float32Array(numBins);
        this.fbBuffer = new Int32Array(numBins);
        this.numFb = 0;

        this.maxDb = -180;

        this.anal = new MagnitudesHistory(numBins, 10, -60);
    }

    initConv() {
        this.conv = new FFTConvFilter(this.audioContext, this.numBins)
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter();
        hipass.type = "highpass"; hipass.frequency.setValueAtTime(this.binToFreq * 2, this.now());
        hipass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        const lopass = this.audioContext.createBiquadFilter();
        lopass.type = "lowpass"; lopass.frequency.setValueAtTime(15000, this.now());
        lopass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        this.hipass = hipass; this.lopass = lopass;
    }

    initFFTNode() {
        this.fftNode = this.audioContext.createAnalyser();
        this.fftNode.fftSize = 512;
        this.fftNode.smoothingTimeConstant = 0.9;

        this.outFftNode = this.audioContext.createAnalyser();
        this.outFftNode.fftSize = this.fftNode.fftSize;
        this.outFftNode.smoothingTimeConstant = this.fftNode.smoothingTimeConstant;
    }

    async initInput() {
        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await navigator.mediaDevices.getUserMedia(constraints);
           this.input = this.audioContext.createMediaStreamSource(stream);
           console.log("[Input] Got access to input device")
         } catch(err) {
            console.error("[Input] Can't access user input device")
         }
    }

    connectAll() {
        const signal = this.input
        .connect(this.hipass)
        .connect(this.lopass)

        signal.connect(this.fftNode)

        signal.connect(this.conv.node).connect(this.autoGain.gain).connect(this.autoGain.limiter).connect(this.audioContext.destination)
        this.conv.connect(this.outFftNode);
    }

    update(){
        requestAnimationFrame(()=>this.update())
        this.updateSpectrum();
        this.autoGain.updateGain(this.anal.maxDb);
        this.conv.updateKernel(this.anal.magReductionsAmp);
        this.drawSpectrum();
    }

    updateSpectrum() {
        this.fftNode.getFloatFrequencyData(this.fftBuffer);
        this.outFftNode.getFloatFrequencyData(this.outFftBuffer);
        this.maxDb = this.fftNode.maxDecibels;
        this.anal.analyseSpectrum(this.fftBuffer, this.maxDb);
    }

    // graphs
    drawSpectrum() {
        //const pitchI = this.anal.fbIndexes[0];
        const pitchI = this.fftBuffer.indexOf(this.anal.maxDb);
        const bgColor = this.anal.nbPeaks > 0 ? Graph.pitchAmpToHSLA(pitchI * this.binToFreq, this.fftBuffer[pitchI]) : 'rgba(100,100,100,0.3)';
        //Graph.drawBg("canvas#spectrum", 'rgba(0,0,0,0.5)')
        Graph.drawBg("canvas#spectrum", bgColor)
        Graph.drawGain("canvas#spectrum", "red", this.autoGain.getCurrentValue(), -100, 80)

        //Graph.drawSmoothLogLine("canvas#spectrum", this.fftBuffer, -100, 0);
        Graph.drawLogLine("canvas#spectrum", this.fftBuffer, -100, 0);

        // Graph.drawSmoothLogLine("canvas#spectrum", this.outFftBuffer, -100, 0, "rgba(0,0,0,0.1)");
        
        //Graph.drawSmoothLogLineInverted("canvas#spectrum", this.anal.magReductions, -100, 0, "rgba(0,0,0,0.5)");
        Graph.drawLogLineInverted("canvas#spectrum", this.anal.magReductions, -100, 0, "rgba(0,0,0,0.5)");

        /*
        Graph.drawHLine("canvas#spectrum",  this.anal.averageDb, -100, 0, "green");
        //Graph.drawHLine("canvas#spectrum", this.anal.maxDb, -100, 0, "red");
        Graph.drawHLine("canvas#spectrum", this.anal.peakThr, -100, 0, "blue");
        */

        Graph.drawVerticals("canvas#spectrum", this.anal.nbPeaks, this.anal.peakIndexes, -80, 0, "rgba(255,255,255,0.01)")
        Graph.drawVerticals("canvas#spectrum", this.numFb, this.fbBuffer, -80, 0, "rgba(255,0,0,0.5)")

        //Graph.drawBg("canvas#slope", "white")
        //Graph.drawLine("canvas#slope", this.conv.irBuffer.getChannelData(0), -1, 1, "rgba(255,0,0,1)");
    }
}

module.exports = SquidbackFFTProcess