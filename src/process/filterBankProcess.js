const SquidbackCommonProcess = require ('./common.js');
const { 
    MelRebin,
    NotchFilterBank
} = require("../audio-processors/index.js")

class SquidbackFilterBankProcess extends SquidbackCommonProcess {

    initSubclassNodes() {
        this.initFilterBank();
        this.drawFunctions = [
            'drawColorBg',
            'drawGain',
            'drawSlopes',

            'drawMelInput',
            //'drawOutputSpectrumForFilterFreqs',
            'drawFR'
        ]
    }

    initFilterBank(div=5, minFreq, maxFreq) {
        const interval = Math.pow(2, 1/this.octaveDivisions);
        const q = Math.sqrt(interval) / (interval - 1);
        this.filterBank = new NotchFilterBank(this.audioContext, this.anal.mel.filter.freqs, q)
        this.freqResponse = new Float32Array(this.anal.numFilters);
        this.filterBank.getFrequencyResponse(this.anal.mel.filter.freqs, this.freqResponse);
    }

    connectAll() {
        if(this.inputDevice) this.inputDevice.connect(this.input)

        const signal = this.input
        //.connect(this.autoGain.gain)
        .connect(this.hipass)
        .connect(this.lopass)

        signal.connect(this.fftNode)

        this.filterBank.connectInput(signal)
        .connect(this.autoGain.limiter)
        .connect(this.panner)
        .connect(this.output)
        .connect(this.audioContext.destination)

        this.filterBank.connect(this.outFftNode)
        //this.input.connect(this.output).connect(this.audioContext.destination)
    }

    update(){
        requestAnimationFrame(()=>this.update())
        this.updateSpectrum();
        this.autoGain.updateGain(this.anal.maxDb);
        this.draw();
    }

    updateSpectrum() {
        //this.fftNode.getFloatFrequencyData(this.fftBuffer);
        this.fftNode.getByteFrequencyData(this.fftBuffer);
        this.outFftNode.getFloatFrequencyData(this.outFftBuffer);
        this.maxDb = this.fftNode.maxDecibels;
        this.minDb = this.fftNode.minDecibels;
        this.anal.analyseSpectrum(this.fftBuffer, this.minDb, this.maxDb);

        // this.filterGains = this.mel.filterSpectrum(this.anal.magReductions);
        this.filterBank.setGains(this.anal.magReductions);
        this.filterBank.getFrequencyResponse(this.anal.mel.filter.freqs, this.freqResponse);
        const correctionBaselineAmp = Math.max(...this.freqResponse);
        this.anal.normalizeReductions(correctionBaselineAmp);
        this.filterGains = this.anal.magReductions;
        this.filterBank.setGains(this.filterGains);
        this.filterBank.getFrequencyResponse(this.anal.mel.filter.freqs, this.freqResponse);
    }

    // graphs

    drawMelInput(opts) {
        //this.graph.drawLogLine("canvas#spectrum", this.fftBuffer.slice(this.minFilterBin, this.maxFilterBin), opts.minDb, opts.maxDb,"rgba(255,255,255,0.5)", 1);
        this.graph.fillLine(this.anal.melMagDb, opts.minDb, opts.maxDb,"rgba(255,255,255,0.5)");
    }
    /*drawOutputSpectrumForFilterFreqs(opts) {
        //this.graph.drawLogLine("canvas#spectrum", this.outFftBuffer.slice(this.minFilterBin, this.maxFilterBin), opts.minDb, opts.maxDb, "rgba(0,0,0,0.05)", 1);
        this.graph.drawFreqFFT("canvas#spectrum", this.outFftBuffer.slice(this.minFilterBin, this.maxFilterBin), opts.minDb, opts.maxDb, "rgba(0,0,0,0.05)", 1, false);
    }*/

    drawSlopes(opts) {
        //console.log(this.anal.msd.slopes)
        this.graph.drawLine(this.anal.msd.slopes,  -100, 100, "rgba(0,0,0,0.5)")
        this.graph.drawGrayscaleVerticals(this.anal.msd.slopes,  0.01)
    }

    drawFR(opts) {
        // freqResponse is updated in updateSpectrum
        this.graph.drawSmoothCQDB(this.freqResponse , -80, 0, "rgba(30,30,30,0.4)", this.anal.mel.filter.freqs);
        this.filterBank.forEachFreqResponse(this.anal.mel.filter.freqs, spectrum=>{
            this.graph.drawSmoothCQDB(spectrum, -80, 0, "rgba(0,0,0,0.3)", this.anal.mel.filter.freqs);
        })
    }
}

module.exports = SquidbackFilterBankProcess