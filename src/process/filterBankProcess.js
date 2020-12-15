const SquidbackCommonProcess = require ('./common.js');
const { NotchFilterBank } = require("../audio-processors/notchFilterBank.js")

class SquidbackFilterBankProcess extends SquidbackCommonProcess {

    initSubclassNodes() {
        this.initFilterBank();
        this.drawFunctions = [
            'drawColorBg',
            'drawGain',
            //'drawSlopes',

            'drawMelInput',
            //'drawOutputSpectrumForFilterFreqs',
            'drawFR'
        ]
    }

    initFilterBank() {
        //const interval = Math.pow(2, 1/this.octaveDivisions);
        //const q = Math.sqrt(interval) / (interval - 1);

        const intervals = this.anal.mel.freqs.map((v,i,a)=>
            i==0 ? (a[i+1] / v) : (v / a[i-1])
        );
        const  q = intervals.map(v=>Math.sqrt(v) / (v - 1))
        this.filterBank = new NotchFilterBank(this.audioContext, this.anal.mel.freqs, q)
        this.freqResponse = new Float32Array(this.anal.numFilters);
        this.filterBank.getFrequencyResponse(this.anal.mel.freqs, this.freqResponse);
    }

    connectAll() {
        if(this.inputDevice) this.inputDevice.connect(this.input)

        const signal = this.input
        .connect(this.hipass)
        .connect(this.lopass)

        signal.connect(this.fftNode)

        this.filterBank.connectInput(signal)
        .connect(this.autoGain.gain)
        .connect(this.autoGain.limiter)
        .connect(this.panner)
        .connect(this.output)
        .connect(this.audioContext.destination)

        //this.filterBank.connect(this.outFftNode)
        //this.input.connect(this.output).connect(this.audioContext.destination)
    }

    update(){
        requestAnimationFrame(()=>this.update())
        this.updateSpectrum();
        this.autoGain.updateGain(this.currentVolume);
        this.draw();
    }

    updateSpectrum() {
        this.fftNode.getFloatFrequencyData(this.fftBuffer);
        //this.fftNode.getByteFrequencyData(this.fftBuffer);
        //this.outFftNode.getFloatFrequencyData(this.outFftBuffer);
        this.anal.analyseSpectrum(this.fftBuffer, this.fftNode.minDecibels, this.fftNode.maxDecibels);
        this.currentVolume = this.anal.maxDb;// - this.numBinsDb;
        const correctionBaselineAmp = this.findCorrectionBaseline();
        this.anal.normalizeReductions(correctionBaselineAmp);
        this.filterGains = this.anal.magReductions;
        this.filterBank.setGains(this.filterGains);
        this.filterBank.getFrequencyResponse(this.anal.mel.freqs, this.freqResponse);
    }

    findCorrectionBaseline() {
        this.filterBank.setGains(this.anal.magReductions);
        this.filterBank.getFrequencyResponse(this.anal.mel.freqs, this.freqResponse);
        return Math.max(...this.freqResponse);
    }

    // graphs

    drawMelInput(opts) {
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
        this.graph.drawSmoothCQDB(this.freqResponse , -40, 0, "rgba(30,30,30,0.4)", this.anal.mel.freqs);
        this.filterBank.forEachFreqResponse(this.anal.mel.freqs, spectrum=>{
            this.graph.drawSmoothCQDB(spectrum, -40, 0, "rgba(0,0,0,0.3)", this.anal.mel.freqs);
        })
    }
}

module.exports = SquidbackFilterBankProcess