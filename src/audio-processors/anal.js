const d3_peaks = require("d3-peaks")
const {MelRebin} = require('./mel-rebin.js')

function dbamp(db) {
    return Math.pow(10, db * 0.05)
}
function ampdb(amp) {
    return 20 * Math.log10(amp)
}

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

class PeakFinder {
    constructor(maxPeaks = 10, numBins) {
        this.numBins = numBins;
        this.maxPeaks = maxPeaks;
        this.nbPeaks = 0;
        //this.peakThr = 0;
        this.peakIndexes = new Int32Array(numBins);
        this.peakHistorySize = 1024 * 10;
        this.rPeakHistorySize = 1 / this.peakHistorySize;
        this.peakHistory = new Float32Array(numBins * this.peakHistorySize);
        this.peakPersistence = new Float32Array(numBins);
        this.minPeakPersistence = this.peakHistorySize * 0.75;
        // this.peakWidth = 11;

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

    /*findPeaks_custom(){
        var mags = this.magDb;
        const thr = this.peakThr;
        let blockMax = -Infinity;
        let blockMaxBin = -1;
        let l_bin = 0; let r_bin = 0; 

        const checkBlock = ()=>{
            r_bin = l_bin;
            blockMax = -Infinity;
            blockMaxBin = -1;
            while(r_bin < l_bin + this.peakWidth * 2) {
                if(mags[r_bin] >= thr && mags[r_bin] > blockMax) { blockMax = mags[r_bin]; blockMaxBin = r_bin }
                r_bin++;
            }
        }

        checkBlock();
        let bin = r_bin - this.peakWidth;
        if(blockMaxBin == bin) this.addPeak(bin);
        r_bin++; l_bin++; bin++;

        while(r_bin < this.numBins){
            if(mags[r_bin] >= thr && mags[r_bin] > blockMax) { 
                blockMax = mags[r_bin]; blockMaxBin = r_bin
            } else {
                if(blockMaxBin > 0 && blockMaxBin < l_bin) checkBlock();
                if(blockMaxBin == bin) this.addPeak(bin);
            }
            r_bin++; l_bin++; bin++;
        }
    }*/
}

class MagnitudesHistory {

    constructor(sampleRate, fftSize, /*historySize,*/ octaveDivisions = 5, minFreq, maxFreq, minPeakThr=-40) {
        // this.numBins = numBins;
        // this.rNumBins = 1.0 / numBins;
        // this.historySize = historySize;
        // this.magDb = new Float32Array(numBins);
        // this.magScale = 2 / this.numBins;

        this.mel = new MelRebin(sampleRate, fftSize, octaveDivisions, minFreq, maxFreq);
        this.numFilters = this.mel.filter.freqs.length;
        this.melMagDb = new Float32Array(this.numFilters)
        this.avgThr = 0.5;
        // this.resetDbStats();

        this.historySize = 50;
        this.msd = new MSD(this.numFilters, this.historySize);

        this.maxPeaks = 10;
        this.minPeakThr = minPeakThr;
        this.peakFinder = new PeakFinder(this.maxPeaks, this.numFilters)

        this.fbIndexes = new Int32Array(this.numFilters);
        this.nbFb = 0;
        this.fbHistory = new Int32Array(25);
        this.fbHistoryPos = 0;

        this.magReductions = new Float32Array(this.numFilters);
        this.magReductionsAmp = new Float32Array(this.numFilters);
        this.correctionSmoothing = 1-1e-3;//-3;
        this.maxCorrection = -300;
        this.correctionMemory = new Float32Array(this.numFilters);
        this.correctionMemoryFactor = 1e-7;
    }

    /*resetDbStats() {
        this.minDb = -180;
        this.maxDb = -180;
        this.averageDb = 0;
    }*/

    shiftHistory() {
        // this.resetDbStats();
        this.msd.prepareNewBlock();
        this.peakFinder.prepareNewBlock();
    }

    analyseSpectrum(buffer, minDb, maxDb) {
        this.shiftHistory();
        //buffer.forEach((mag, binIndex)=>this.addForBin(binIndex, mag));

        const [min, max, nMax, avg] = this.mel.filterUint8SpectrumInPlace(buffer, this.melMagDb, minDb, maxDb); // in place to avoid garbage collection
        this.averageDb = avg; this.minDb = min; this.maxDb = max;
        this.maxBin = nMax;
        // this.averageDb = this.averageDb * this.rNumBins;
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
        //this.findFeedbackCandidates();
        //this.normalizeCorrections();
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
            //this.magReductionsAmp[bin] = 0;
            this.correctionMemory[bin] += correction * this.correctionMemoryFactor;
        }
        else if(correction < 0 ){
            this.magReductions[bin] = correction + this.correctionMemory[bin];
            //this.magReductionsAmp[bin] = dbamp(this.magReductions[bin]);
            this.correctionMemory[bin] += correction * this.correctionMemoryFactor;
        } else {
            this.magReductions[bin] = 0;
            //this.magReductionsAmp[bin] = 1;
        }
    }

    normalizeReductions(ampBaseline) {
        const minCorrectionAmp = 1/ampBaseline;
        const minCorrection = ampdb(minCorrectionAmp);
        //console.log(ampBaseline, minCorrectionAmp, minCorrection)
        for(const bin in this.magReductions) {
            this.magReductions[bin] += minCorrection;
            //this.magReductionsAmp[bin] /= minCorrectionAmp;
            if(this.magReductions[bin] > 0) this.magReductions[bin] = 0
            //if(this.magReductionsAmp[bin] > 1) this.magReductions[bin] = 1
        }
    }

}

module.exports = { MSD, MagnitudesHistory }