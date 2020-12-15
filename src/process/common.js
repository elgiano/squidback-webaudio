const SquidbackGraph = require('./graph.js')
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
        this.fftBuffer = new Uint8Array(numBins);
        //this.outFftBuffer = new Float32Array(numBins);
        //this.fbBuffer = new Int32Array(numBins);

        this.maxVisualDb = 0;
        this.minVisualDb = -180;
    }

    initAnal(octaveDivisions = 5, minFreq, maxFreq) {
        minFreq = minFreq || this.hipassFreq / 1.2
        maxFreq = maxFreq || this.lopassFreq * 1.2
        this.octaveDivisions = octaveDivisions;
        this.anal = new MagnitudesHistory(this.audioContext.sampleRate, this.fftSize, 5, minFreq, maxFreq, -60);

        const minBin = minFreq / this.binToFreq;
        const maxBin = maxFreq / this.binToFreq;
        this.minFilterBin = Math.round(minBin); this.maxFilterBin = Math.round(maxBin);
        this.graph.setupFreqAxis(minFreq, maxFreq, this.binToFreq)
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter(); hipass.type = "highpass"; 
        this.hipassFreq = Math.min(this.binToFreq * 2, 80)
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
        this.fftNode.minDb = -180
        this.fftNode.maxDb = 0

        /*this.outFftNode = this.audioContext.createAnalyser();
        this.outFftNode.fftSize = this.fftNode.fftSize;
        this.outFftNode.smoothingTimeConstant = this.fftNode.smoothingTimeConstant;*/

        this.fftSize = fftSize;
        this.numBins = this.fftNode.frequencyBinCount;
        this.numBinsDb = Math.log10(this.fftNode.frequencyBinCount) * 20;
        this.binToFreq = this.audioContext.sampleRate / 2 / this.numBins;
    }

    async initIO() {
        this.input = this.audioContext.createGain();
        this.output = this.audioContext.createGain();
        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await navigator.mediaDevices.getUserMedia(constraints);
           this.inputDevice = this.audioContext.createMediaStreamSource(stream);
           console.log("[Input] Got access to input device")
         } catch(err) {
            console.error("[Input] Can't access user input device")
         }
    }

    // draw graphs

    // draw modularly, call every function in this.drawFunction
    // intended to implement gui options to enable/disable graphs
    draw() {
        const smooth = 0.5;
        this.minVisualDb = smooth * this.minVisualDb + (1-smooth)*(this.anal.minDb - 20);
        this.maxVisualDb = smooth * this.maxVisualDb + (1-smooth)*(this.anal.maxDb + 20);
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