const {MelRebin} = require('./mel-rebin.js')

function ampdb(amp) { return 20 * Math.log10(amp) }

/*
const d3_peaks = require("d3-peaks")
class PeakFinder {
    constructor(maxPeaks = 10, numBins) {
        this.numBins = numBins;
        this.maxPeaks = maxPeaks;
        this.nbPeaks = 0;

        this.peakIndexes = new Int32Array(numBins);
        this.peakHistorySize = 1024 * 10;
        this.rPeakHistorySize = 1 / this.peakHistorySize;
        this.peakHistory = new Float32Array(numBins * this.peakHistorySize);
        this.peakPersistence = new Float32Array(numBins);
        this.minPeakPersistence = this.peakHistorySize * 0.75;

        const ricker = d3_peaks.ricker;
        this.findPeaksFn = d3_peaks.findPeaks()
            //.kernel(ricker)
            //.gapThreshold(10)
            //.minSNR(1)
            //.widths([1,2,3]);
    }

    prepareNewBlock() {
        const lastPeakMemBlock = this.numBins * (this.peakHistorySize - 1);
        for(let bin = 0; bin < this.numBins; ++bin)
            this.peakPersistence[bin] -= this.peakHistory[lastPeakMemBlock + bin];  
        this.nbPeaks = 0;
        this.peakHistory.copyWithin(this.numBins,0);
        this.peakHistory.fill(0, 0, this.numBins);
    }

    findPeaks(data, peakCallback) {
        var peaks = this.findPeaksFn(data);
        peaks.forEach(p=>this.addPeak(p.index, peakCallback));
        // console.log(peaks);
    }

    addPeak(binIndex, callback=()=>{}){
        this.peakIndexes[this.nbPeaks] = binIndex;
        this.peakHistory[binIndex] = this.rPeakHistorySize;
        this.peakPersistence[binIndex] += this.rPeakHistorySize;
        this.nbPeaks++;
        callback(binIndex)
    }
}*/


class MSD {

    constructor(numBins, historySize){
        this.numBins = numBins;
        this.historySize = historySize;
        this.msd = new Float32Array(numBins);
        this.magHistory = new Float32Array(numBins * historySize);
        this.magHistory.fill(0);
        //this.magDiffDiffHistory = new Float32Array(numBins * historySize);
        //this.lastMagnitudes = new Float32Array(numBins);
        this.lastMagDiff = new Float32Array(numBins);
        //this.magDiffDiffNormalize = 1 / historySize;
        this.slopes = new Float32Array(numBins);
    }

    prepareNewBlock() {
        /*const lastMemBlock =  this.numBins * (this.historySize - 1);
        for(let bin = 0; bin < this.numBins; ++bin) {
            this.msd[bin] -= this.magDiffDiffHistory[lastMemBlock + bin];
        }*/
        this.magHistory.copyWithin(this.numBins,0);
        //this.magDiffDiffHistory.copyWithin(this.numBins,0);
    }

    analyzeSpectrum(spectrum, smoothing=0.99) {
        for(let bin = 0; bin < spectrum.length; ++bin)
            this.addForBin(bin, spectrum[bin], smoothing)
    }

    addForBin(binIndex, magDb, smoothing) {
        this.magHistory[binIndex] = (1-smoothing) * magDb + smoothing * this.magHistory[binIndex];
        //const magDiff =  magDb - this.magHistory[binIndex + this.numBins];
        //const magDiffDiff = Math.pow(magDiff - this.lastMagDiff[binIndex], 2);
        //this.magDiffSum[binIndex] += (1-smoothing) * magDiff + smoothing * this.magDiffHistory[binIndex];
        //this.msd[binIndex] += (1-smoothing) * magDiffDiff + smoothing * this.magDiffDiffHistory[binIndex];
        //this.magDiffHistory[binIndex] = magDiffDiff;
        //this.magDiffDiffHistory[binIndex] = magDiffDiff;
        //this.lastMagDiff[binIndex] = magDiff;
        
        this.slopes[binIndex] = this.getSlope(binIndex);
    }

    getSlope(binIndex, framesAgo=(this.historySize-1)) { 
        return this.magHistory[binIndex] - this.magHistory[framesAgo*this.numBins + binIndex]
    }
}

class MagnitudesHistory {

    constructor(sampleRate, fftSize, /*historySize,*/ octaveDivisions = 5, minFreq, maxFreq, minPeakThr=-40) {
        this.mel = new MelRebin(sampleRate, fftSize, octaveDivisions, minFreq, maxFreq);
        this.numFilters = this.mel.filter.freqs.length;
        this.melMagDb = new Float32Array(this.numFilters)
        this.avgThr = 0.5;

        this.historySize = 50;
        this.msd = new MSD(this.numFilters, this.historySize);

        /*this.maxPeaks = 10;
        this.minPeakThr = minPeakThr;
        this.peakFinder = new PeakFinder(this.maxPeaks, this.numFilters)*/

        this.magReductions = new Float32Array(this.numFilters);
        this.correctionSmoothing = 1-1e-3;//-3;
        this.maxCorrection = -300;
        this.correctionMemory = new Float32Array(this.numFilters);
        this.correctionMemoryFactor = 1e-7;
    }

    analyseSpectrum(buffer, minDb, maxDb) {
        this.msd.prepareNewBlock();
        //this.peakFinder.prepareNewBlock();

        const [min, max, nMax, avg] = this.mel.filterUint8SpectrumInPlace(buffer, this.melMagDb, minDb, maxDb); // in place to avoid garbage collection
        this.averageDb = avg; this.minDb = min; this.maxDb = max;
        this.maxBin = nMax;

        this.peakThr = this.averageDb + (max - this.averageDb) * this.avgThr;
        // console.log(this.peakThr, this.averageDb, this.avgThr, max)
        if(this.peakThr < this.minPeakThr) this.peakThr = this.minPeakThr;

        /*this.peakFinder.findPeaks(this.melMagDb, (peakIndex)=>
            this.correctMagnitude(peakIndex)
        );*/

        this.msd.analyzeSpectrum(this.melMagDb);

        for(const bin in this.melMagDb) {
            this.correctMagnitude(bin)
        }
    }

    correctMagnitude(bin) {
        //this.magReductions[bin] += (this.peakThr - this.magDb[bin]) / 2
        let correction = (this.peakThr - this.melMagDb[bin]);// * this.peakFinder.peakPersistence[bin];// * 2;
        // console.log(correction, this.peakThr, this.melMagDb[bin])

        const slope = this.msd.slopes[bin]
        if(correction < 0){
            if(slope < -5)
                correction *= -2
            else if(slope > 5)
                correction *= 0.01
        }

        correction = this.correctionSmoothing * this.magReductions[bin] + (1-this.correctionSmoothing) * correction

        if(correction < this.maxCorrection ) {
            this.magReductions[bin] = correction + this.correctionMemory[bin];
            this.correctionMemory[bin] += correction * this.correctionMemoryFactor;
        }
        else if(correction < 0 ){
            this.magReductions[bin] = correction + this.correctionMemory[bin];
            this.correctionMemory[bin] += correction * this.correctionMemoryFactor;
        } else {
            this.magReductions[bin] = 0;
        }
    }

    normalizeReductions(ampBaseline) {
        const minCorrectionAmp = 1/ampBaseline;
        const minCorrection = ampdb(minCorrectionAmp);
        //console.log(ampBaseline, minCorrectionAmp, minCorrection)
        for(const bin in this.magReductions) {
            this.magReductions[bin] += minCorrection;
            if(this.magReductions[bin] > 0) this.magReductions[bin] = 0
        }
    }

}

module.exports = { MagnitudesHistory }