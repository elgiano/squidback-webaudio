const SquidbackGraph = require('../graph.js')
const Graph = new SquidbackGraph();
const { FFTConvFilter } = require("../audio-processors/index.js")
const SquidbackCommonProcess = require ('./common.js');

class SquidbackFFTProcess extends SquidbackCommonProcess{

    initSubclassNodes() {
        this.initConv();
        this.maxDb = 180;
        this.drawFunctions = [
            'drawColorBg',
            'drawGain',
            'drawInputSpectrumSmooth',
            'drawOutputSpectrumSmooth',
            'drawPeaks',
            'drawFbs',
            'drawMagReductions',
            'drawIR'
        ]
    }

    initConv() {
        this.conv = new FFTConvFilter(this.audioContext, this.numBins)
    }

    connectAll() {
        const signal = this.input
        .connect(this.autoGain.gain)
        .connect(this.hipass)
        .connect(this.lopass)

        signal.connect(this.fftNode)

        this.conv.connectInput(signal)
        .connect(this.autoGain.limiter)
        .connect(this.panner)
        .connect(this.output)
        .connect(this.audioContext.destination)
        this.conv.connect(this.outFftNode);
    }

    update(){
        requestAnimationFrame(()=>this.update())
        this.updateSpectrum();
        this.autoGain.updateGain(this.anal.maxDb);
        //if(this.count == 0)
        for(const bin in this.anal.magReductionsAmp){
            this.anal.magReductionsAmp[bin] = dbamp(this.anal.magReductions[bin])
        }
        this.conv.updateKernel(this.anal.magReductionsAmp);
        //this.count = ((this.count || 0) + 1) % 1;
        this.draw();
    }

    updateSpectrum() {
        this.fftNode.getFloatFrequencyData(this.fftBuffer);
        this.outFftNode.getFloatFrequencyData(this.outFftBuffer);
        this.maxDb = this.fftNode.maxDecibels;
        this.anal.analyseSpectrum(this.fftBuffer, this.maxDb);

        // const mel = this.mel.filterSpectrum(this.anal.magReductions)
        // console.log("mel", mel)
    }

    // graphs

    drawIR() {
        // TODO: draw on aux graph
        Graph.drawBg("canvas#slope", "white")
        Graph.drawLine("canvas#slope", this.conv.irBuffer.getChannelData(0), -1, 1, "rgba(255,0,0,1)");
    }
}

module.exports = SquidbackFFTProcess