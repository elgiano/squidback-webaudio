const SquidbackGraph = require('./graph.js')
const adapter = require('webrtc-adapter')
const { AutoGain, MagnitudesHistory } = require("../audio-processors/index.js")

// this class is here because I've been experimenting with different processes
// (notably one where I attempted to filter by iffting a correction spectrum to time domain and convolving it to the signal
// which was computationally cheap but didn't work well sound-wise)
class SquidbackCommonProcess {

    // subclass responsibilites
    initSubclassNodes() {}
    connectAll() {}
    update() {}

    // parent

    constructor(audioContext, canvas) {
        this.audioContext = audioContext;
        this.graph = new SquidbackGraph(canvas);
        this.drawFunctions = [
            'drawColorBg',
            'drawGain',
            'drawInputSpectrum',
        ]

        this.lopassFreq = 10000
        this.minDb = -180;
        this.maxDb = 10;

    }

    async start(fftSize = 1024) {
        this.initFFTNode(fftSize);
        this.initExtremesFilters();
        this.initBuffers();
        this.initAnal();
        this.initGain();
        this.initPanner();
        this.initSubclassNodes();
        await this.initIO();
        this.connectAll();
        requestAnimationFrame(()=>this.update());
    }

    initGain() {
        this.autoGain = new AutoGain(this.audioContext)
    }

    initPanner() {
        this.panner = this.audioContext.createStereoPanner()
        this.panner.pan.value = 0;
    }

    initBuffers(minFreq, maxFreq) {
        const numBins = this.numBins;
        //this.fftBuffer = new Uint8Array(numBins);
        this.fftBuffer = new Float32Array(numBins);
        //this.outFftBuffer = new Float32Array(numBins);
        //this.fbBuffer = new Int32Array(numBins);

        this.maxVisualDb = 0;
        this.minVisualDb = -180;
    }

    initAnal(numFilters = 30, minFreq, maxFreq) {
        minFreq = minFreq || this.hipassFreq
        maxFreq = maxFreq || this.lopassFreq * 1.2
        this.anal = new MagnitudesHistory(this.audioContext.sampleRate, this.fftSize, numFilters, minFreq, maxFreq, -60);

        this.graph.setupFreqAxis(minFreq, maxFreq, this.binToFreq)
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter(); hipass.type = "highpass"; 
        this.hipassFreq = this.binToFreq * 2
        hipass.frequency.setValueAtTime(this.hipassFreq,  this.audioContext.currentTime);
        hipass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        const lopass = this.audioContext.createBiquadFilter(); lopass.type = "lowpass"; 
        lopass.frequency.setValueAtTime(this.lopassFreq,  this.audioContext.currentTime);
        lopass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        this.hipass = hipass; this.lopass = lopass;
    }

    initFFTNode(fftSize) {
        this.fftNode = this.audioContext.createAnalyser();
        this.fftNode.fftSize = fftSize;
        this.fftNode.smoothingTimeConstant = 0.9;
        this.fftNode.minDecibels = -180
        this.fftNode.maxDecibels = 0

        /*this.outFftNode = this.audioContext.createAnalyser();
        this.outFftNode.fftSize = this.fftNode.fftSize;
        this.outFftNode.smoothingTimeConstant = this.fftNode.smoothingTimeConstant;*/

        this.fftSize = fftSize;
        this.numBins = this.fftNode.frequencyBinCount;
        this.numBinsDb = Math.log10(this.fftNode.frequencyBinCount) * 20;
        this.binToFreq = this.audioContext.sampleRate / 2 / this.numBins;
    }

    async getUserMedia(constraints) {
        let stream;
        if(navigator.mediaDevices) {
            stream = await navigator.mediaDevices.getUserMedia(constraints)
        } else {
            stream = await navigator.getUserMedia(constraints)
        }
        return stream

    }

    async initIO() {
        this.input = this.audioContext.createGain();
        this.output = this.audioContext.createGain();

        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await this.getUserMedia(constraints);
           this.inputDevice = this.audioContext.createMediaStreamSource(stream);
           console.log("[Input] Got access to input device")
         } catch(err) {
            console.error("[Input] Can't access user input device", err)
         }
    }

    // draw graphs

    // draw modularly, call every function in this.drawFunction
    // intended to implement gui options to enable/disable graphs
    draw() {
        const smooth = 0.5;
        this.minVisualDb = -120;//smooth * this.minVisualDb + (1-smooth)*(this.anal.minDb * 1.01);
        this.maxVisualDb = -10;//smooth * this.maxVisualDb + (1-smooth)*(this.anal.maxDb + 10);
        const drawOptions = {
            minDb: this.minVisualDb, maxDb: this.maxVisualDb
        };
        this.drawFunctions.forEach(fn=>this[fn](drawOptions))
    }

    drawBg() {
        this.graph.drawBg('rgba(0,0,0,0.5)')
    }

    findPitch() {
        let max = -Infinity; let maxI = -1;
        for(const bin in this.fftBuffer) {
            const val = this.fftBuffer[bin]
            if(val < max) continue
            max = val; maxI = bin; 
        }
        return [maxI * this.binToFreq, max / 255 * this.fftNode.maxDecibels + this.fftNode.minDecibels]
    }

    drawColorBg() {
        const [pitch, amp] = this.findPitch();
        const bgColor = this.graph.pitchAmpToHSLA(pitch, amp);
        this.graph.drawBg(bgColor)
    }

    drawGain() {
        this.graph.drawGain("red", this.autoGain.getCurrentValueDb(), this.autoGain.minGain, this.autoGain.maxGain);
    }

    drawInputSpectrum(opts) {
        this.graph.drawLogLine(this.fftBuffer, opts.minDb, opts.maxDb,"rgba(255,255,255,0.5)", 1);
    }
    drawInputSpectrumSmooth(opts) {
        this.graph.drawAvgLogLine(this.fftBuffer, opts.minDb, opts.maxDb,"rgba(255,255,255,0.5)", 3, 2);
    }
    /*drawOutputSpectrum(opts) {
        this.graph.drawLogLine(this.outFftBuffer, opts.minDb, opts.maxDb, "rgba(0,0,0,0.05)", 1);
    }
    drawOutputSpectrumSmooth(opts) {
        this.graph.drawAvgLogLine(this.outFftBuffer, opts.minDb, opts.maxDb, "rgba(0,0,0,0.05)", 3, 2);
    }*/

    drawStatLines(opts) {
        this.graph.drawHLine(this.anal.averageDb, opts.minDb, opts.maxDb, "green");
        this.graph.drawHLine(this.anal.maxDb, opts.minDb, opts.maxDb, "red");
        this.graph.drawHLine(this.anal.peakThr, opts.minDb, opts.maxDb, "blue");
    }

    drawPeaks() {
        this.graph.drawLinVerticals(this.anal.peakFinder.nbPeaks, this.anal.peakFinder.peakIndexes)
    }

    clearGraphCache() { this.graph.clearCache() }
}

module.exports = SquidbackCommonProcess