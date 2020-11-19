const d3_peaks = require("d3-peaks")

function dbamp(db) {
    return Math.pow(10, db * 0.05)
}

class MSD {

    constructor(numBins, historySize){
        this.msd = new Float32Array(numBins);
        this.magDiff = new Float32Array(numBins);
        this.magDiffDiffHistory = new Float32Array(numBins * historySize);
        this.lastMagnitudes = new Float32Array(numBins);
        this.lastMagDiff = new Float32Array(numBins);
        this.magDiffDiffNormalize = 1 / historySize;
    }

    addForBin(binIndex, magDb, smoothing) {
        this.magDiff[binIndex] = magDb - this.lastMagnitudes[binIndex];
        const magDiffDiff = Math.pow(this.magDiff[binIndex] - this.lastMagDiff[binIndex],2);
        this.msd[binIndex] += (1-smoothing) * magDiffDiff + smoothing * this.magDiffDiffHistory[binIndex] ;
        this.magDiffDiffHistory[binIndex] = magDiffDiff;

        this.lastMagDiff[binIndex] = this.magDiff[binIndex];
        this.lastMagnitudes[binIndex] = this.magDb[binIndex];
    }

}

class MagnitudesHistory {

    constructor(numBins, historySize, minPeakThr=-40) {
        this.numBins = numBins;
        this.rNumBins = 1.0 / numBins;
        this.historySize = historySize;
        this.magDb = new Float32Array(numBins);
        this.magScale = 2 / this.numBins;

        this.msd = new MSD(this.numBins, historySize);

        this.maxDb = -180;
        this.averageDb = 0;
        this.avgThr = 0.5;

        this.maxPeaks = 10;
        this.nbPeaks = 0;
        this.minPeakThr = minPeakThr;
        this.peakThr = 0;
        this.peakIndexes = new Int32Array(numBins);
        this.peakHistorySize = 100;
        this.peakHistory = new Float32Array(numBins * this.peakHistorySize);
        this.peakPersistence = new Float32Array(numBins);
        this.minPeakPersistence = this.peakHistorySize * 0.75;
        this.peakWidth = 11;

        this.fbIndexes = new Int32Array(this.numBins);
        this.nbFb = 0;
        this.fbHistory = new Int32Array(25);
        this.fbHistoryPos = 0;

        this.magReductions = new Float32Array(this.numBins);
        this.magReductionsAmp = new Float32Array(this.numBins);
        this.magReductionsUnitDb = 0.1;
        this.magReductionsUnitAmp = dbamp(this.magReductionsUnitDb);

        const ricker = d3_peaks.ricker;
          this.findPeaksFn = d3_peaks.findPeaks()
            //.kernel(ricker)
            .gapThreshold(10)
            .minSNR(1)
            .widths([1,2,3]);

        console.log("FFT constructed")
    }

    shiftHistory() {
        //const lastMemBlock =  this.numBins * (this.historySize - 1);
        const lastPeakMemBlock = this.numBins * (this.peakHistorySize - 1);
        for(let bin = 0; bin < this.numBins; ++bin) {
            // this.msd.msd[bin] -= this.msd.magDiffDiffHistory[lastMemBlock + bin];
            this.peakPersistence[bin] -= this.peakHistory[lastPeakMemBlock + bin];  
            //if(this.magReductions[bin]<0) this.magReductions[bin] += this.magReductionsUnitDb
            //if(this.magReductionsAmp[bin]<1) this.magReductions[bin] *= this.magReductionsUnitAmp
        }

        this.maxDb = -180;
        this.averageDb = 0;
        this.nbPeaks = 0;
        // this.msd.magDiffDiffHistory.copyWithin(this.numBins,0);
        this.peakHistory.copyWithin(this.numBins,0);
    }

    analyseSpectrum(buffer, max) {
        this.shiftHistory();
        buffer.forEach((mag, binIndex)=>this.addForBin(binIndex, mag));
        this.averageDb = this.averageDb * this.rNumBins;
        this.peakThr = this.averageDb + (max - this.averageDb) * this.avgThr;
        if(this.peakThr < this.minPeakThr) this.peakThr = this.minPeakThr;

        this.findPeaks();
        this.findFeedbackCandidates();
    }

    addForBin(binIndex, magDb) {
        this.magDb[binIndex] = magDb;
        // this.msd.addForBin(binIndex, magDb, this.smoothing)
        if(magDb > this.maxDb) this.maxDb = magDb;
        this.averageDb += magDb;
        this.peakHistory[binIndex] = 0;
        //this.checkPeak(binIndex);
        this.correctMagnitude(binIndex)
    }

    findPeaks() {
        var peaks = this.findPeaksFn(this.magDb);
        peaks.forEach(p=>this.addPeak(p.index))
        // console.log(peaks);
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

    addPeak(binIndex){
        this.peakIndexes[this.nbPeaks] = binIndex;
        this.peakHistory[binIndex] = 1;
        this.peakPersistence[binIndex]++;
        this.nbPeaks++;
    }
    
    findFeedbackCandidates() {
        this.nbFb = 0;
        let minMsdBin = -1;
        let minMsd = Infinity;

        //console.log(this.peakHistory)
        for (let i = 0; i < this.nbPeaks; ++i) {
            const bin = this.peakIndexes[i];
            if(this.peakPersistence[bin] >= this.minPeakPersistence) {
            /*if(this.msd[bin] < minMsd && this.smoothedMagDb[bin-1] < this.smoothedMagDb[bin] && this.smoothedMagDb[bin+1] < this.smoothedMagDb[bin]) {
                minMsd = this.msd[bin];
                minMsdBin = bin;
            }*/
            this.fbIndexes[this.nbFb++] = bin;
            }
            //this.correctMagnitude(bin)

        }

        
        /*if(minMsdBin >= 0 && minMsd <= 0.01) {
            this.fbHistory[this.fbHistoryPos++] = minMsdBin;
            if(this.fbHistory.reduce((t,v)=>v==minMsdBin?t+1:t,0) >= this.fbHistory.length * 0.5) {
                this.fbIndexes[this.nbFb++] = minMsdBin;
            }
        } else {
            this.fbHistory[this.fbHistoryPos++] = -1;
        }
        this.fbHistory.filter((v,i,a)=>v >= 0 && a.indexOf(v)===i).forEach(bin=>{
            this.fbIndexes[this.nbFb] = bin;
            this.nbFb++;
        })
        if(this.fbHistoryPos > this.fbHistory.length) this.fbHistoryPos = 0;
        */
    }

    correctMagnitude(bin) {
        //this.magReductions[bin] += (this.peakThr - this.magDb[bin]) / 2
        const correction = (this.peakThr - this.magDb[bin]);
        if(correction <= 0 ){
            this.magReductions[bin] = correction;
            this.magReductionsAmp[bin] = dbamp(this.magReductions[bin])
        } else {
            this.magReductions[bin] = 0;
            this.magReductionsAmp[bin] = 1;
        }
    }

}

module.exports = { MSD, MagnitudesHistory }