const FFT = require('fft.js');
const KissFFT = require('kissfft-js')

class FFTConvFilter {
    constructor(audioContext, numBins) {
        this.audioContext = audioContext;
        this.numBins = numBins;
        this.irBuffer = audioContext.createBuffer(2, numBins*2, audioContext.sampleRate)
        this.irBuffer.getChannelData(0)[0] = 1
        this.irBuffer.getChannelData(0)[numBins*2-1] = 1


        /*this.fft = new FFT(numBins * 2);
        this.freqComplexBuffer = this.fft.createComplexArray();
        this.timeComplexBuffer = this.fft.createComplexArray();
        this.realScale = numBins * 2;*/
        this.fftSize = numBins * 2;
        this.fft = new KissFFT.FFTR(numBins * 2);
        this.freqComplexBuffer = new Float32Array(numBins*2);
        this.realScale = 1 / (numBins * 2);


        this.numNodes = 1;
        this.nodes = new Array( this.numNodes);
        this.crossfades = new Array( this.numNodes);
        for(let i = 0; i < this.nodes.length; i++) {
            this.nodes[i] = audioContext.createConvolver();
            this.nodes[i].normalize = false
            this.nodes[i].buffer = this.irBuffer;

            this.crossfades[i] = audioContext.createGain();
            this.crossfades[i].gain.value = 0;
            this.nodes[i].connect(this.crossfades[i]);
        };

        this.activeNode = 0;
        this.crossfades[0].gain.value = 1;
    }

    updateKernel(reductionSpectrum, fade=0.3) {
        this.freqComplexBuffer.fill(0);
        //let avg = 0;
        for(let i = 0, j = 0; i < this.numBins; ++i, j+=2) {
            this.freqComplexBuffer[j] = reductionSpectrum[i] * this.realScale;
            //avg += this.freqComplexBuffer[j]
            this.freqComplexBuffer[j+1] = 0;  
        };
        //avg /= this.numBins;
        this.freqComplexBuffer[0] = 0;//avg;
        //console.log("spec",avg)
        /*this.fft.completeSpectrum(this.freqComplexBuffer)
        this.fft.inverseTransform(this.timeComplexBuffer, this.freqComplexBuffer);
        this.fft.fromComplexArray(this.timeComplexBuffer, this.irBuffer.getChannelData(0));
        this.setBufferWithCrossfade(fade);*/
        const ifft = this.fft.inverse(this.freqComplexBuffer);
        console.log(ifft.length)
        this.irBuffer.getChannelData(0).set(ifft/*.subarray(0,this.numBins)*/);
        this.setBufferWithCrossfade(fade);
    }

    setBufferWithCrossfade(crossFadeTime) {
        const nextNode = (this.activeNode + 1) % this.nodes.length;
        this.nodes[nextNode].buffer = this.irBuffer;
        this.crossfades[this.activeNode].gain.exponentialRampToValueAtTime(1e-45, this.audioContext.currentTime+crossFadeTime);
        this.crossfades[nextNode].gain.exponentialRampToValueAtTime(1, this.audioContext.currentTime+crossFadeTime);
        this.activeNode = nextNode;
    }

    connect(destination) { 
        this.crossfades.forEach(node=>node.connect(destination));
        return destination
    }
    connectInput(input) { 
        this.nodes.forEach(node=>input.connect(node));
        return this
    }
}

module.exports = { FFTConvFilter }