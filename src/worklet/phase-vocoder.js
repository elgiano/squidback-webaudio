"use strict";

const FFTProcessor = require('./fft-processor.js');

const BUFFERED_BLOCK_SIZE = 512;
const DB_MIN = -180;

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

    constructor(numBins, historySize) {
        this.numBins = numBins;
        this.rNumBins = 1.0 / numBins;
        this.historySize = historySize;
        this.magDb = new Float32Array(numBins);
        this.magScale = 2 / this.numBins;

        this.msd = new MSD(this.numBins, historySize);

        this.smoothing = 0.99;
        this.smoothedMagDb = new Float32Array(numBins);
        this.averageDb = 0;
        this.avgThr = 0.5;
        this.smoothedMax = DB_MIN;

        this.maxPeaks = 10;
        this.nbPeaks = 0;
        this.minPeakThr = -40;
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

        console.log("FFT constructed")
    }

    shiftHistory() {
        //const lastMemBlock =  this.numBins * (this.historySize - 1);
        const lastPeakMemBlock = this.numBins * (this.peakHistorySize - 1);
        for(let bin = 0; bin < this.numBins; ++bin) {
            // this.msd.msd[bin] -= this.msd.magDiffDiffHistory[lastMemBlock + bin];
            this.peakPersistence[bin] -= this.peakHistory[lastPeakMemBlock + bin];  
            if(this.magReductions[bin]<0) this.magReductions[bin] += 0.1
        }

        this.peakThr = this.averageDb * this.rNumBins;
        this.peakThr += (this.smoothedMax - this.peakThr) * this.avgThr;
        if(this.peakThr < this.minPeakThr) this.peakThr = this.minPeakThr
        this.averageDb = 0;
        this.smoothedMax = DB_MIN;
        this.nbPeaks = 0;
        // this.msd.magDiffDiffHistory.copyWithin(this.numBins,0);
        this.peakHistory.copyWithin(this.numBins,0);

    }


    addForBin(binIndex, mag) {
        let magDb =  mag == 0 ? DB_MIN : 20 * Math.log10(mag * this.magScale);
        this.smoothedMagDb[binIndex] = (1-this.smoothing) * magDb + this.smoothing * this.smoothedMagDb[binIndex];
        this.magDb[binIndex] = magDb;
        // this.msd.addForBin(binIndex, magDb, this.smoothing)

        this.averageDb += this.smoothedMagDb[binIndex];
        if(this.smoothedMagDb[binIndex] > this.smoothedMax) this.smoothedMax = this.smoothedMagDb[binIndex]
        this.peakHistory[binIndex] = 0;
        //this.checkPeak(binIndex);
    }

    findPeaks(){
        const mags = this.smoothedMagDb;
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
    }

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
            this.correctMagnitude(bin)

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
        this.magReductions[bin] += (this.peakThr - this.magDb[bin]) / 2
    }

}

class PhaseVocoderProcessor extends FFTProcessor {

    constructor(options) {
        options.processorOptions = {
            blockSize: options.processorOptions.blockSize || BUFFERED_BLOCK_SIZE,
        };
        super(options);
        this.magHistory = new MagnitudesHistory(this.numBins, 11);
        this.phases = new Float32Array(this.numBins)
        this.magnitudes = new Float32Array(this.numBins)

        this.listenToNodeMessages();
    }

    listenToNodeMessages() {
        this.responders = {
            getSpectrum: ()=>this.sendMagnitudes(),
            getSlope: ()=>this.sendSlope(),
            getPeaks: ()=>this.sendPeaks(),
            getHisto: ()=>this.sendHisto(),
            getFbCandidates: ()=>this.sendFb(),
            getMaxDb: ()=>this.sendMaxDb(),
            getMagReductions: ()=>this.sendMagReductions()
        };
        this.port.onmessage = event => {
            const method = event.data;
            // console.log("requested", method)
            const responder = this.responders[method];
            if(responder) responder();
        };
        console.log("Processor listening")
    };
    sendMagnitudes(data) { this.port.postMessage(["spectrum", this.magHistory.smoothedMagDb]) }
    sendSlope() { this.port.postMessage(["slope", this.magHistory.msd.msd]) }
    sendHisto() { this.port.postMessage(["histo", this.magHistory.msd.magDiffDiffHistory]) }
    sendPeaks() {  this.port.postMessage(["peaks", this.magHistory.nbPeaks, this.magHistory.peakIndexes]) }
    sendFb() {  this.port.postMessage(["fb",this.magHistory.nbFb, this.magHistory.fbIndexes]) }
    sendMaxDb() {  this.port.postMessage(["maxDb",this.magHistory.smoothedMax]) }
    sendMagReductions() { this.port.postMessage(["reductions", this.magHistory.magReductions]) }

    processFFT(complexSpectrum, output, parameters) {
        this.magHistory.shiftHistory();
        this.computeMagnitudes(complexSpectrum);
        this.magHistory.findPeaks();
        this.magHistory.findFeedbackCandidates();
        let j = 0;
        for(let i=0; i<this.numBins; i++, j+=2) {
            const correction = 10 ** (this.magHistory.magReductions[i] * 0.05);
            this.freqComplexBuffer[j] *= correction;
            this.freqComplexBuffer[j+1] *= correction;
        }

        
    }

    /** Compute squared magnitudes for peak finding **/
    computeMagnitudes(fft) {
        var i = 0, j = 0;
        while (i < this.numBins) {
            
            let real = fft[j];
            let imag = fft[j + 1];
            this.phases[i] = imag/real;
            this.magHistory.addForBin(i, Math.sqrt(real ** 2 + imag ** 2))
            i+=1;
            j+=2;
        }
    }
}

registerProcessor("phase-vocoder-processor", PhaseVocoderProcessor);

