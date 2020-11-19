const FFT = require('fft.js');

class FFTConvFilter {
    constructor(audioContext, numBins) {
        this.numBins = numBins;
        this.irBuffer = audioContext.createBuffer(2, numBins, audioContext.sampleRate)
        this.irBuffer.getChannelData(0)[0] = (1)

        this.fft = new FFT(numBins * 2);
        this.freqComplexBuffer = this.fft.createComplexArray();
        this.timeComplexBuffer = this.fft.createComplexArray();

        this.node = audioContext.createConvolver();
        this.node.normalize = false
        this.node.buffer = this.irBuffer;
    }

    updateKernel(reductionSpectrum) {
        this.freqComplexBuffer.fill(0);
        for(let i = 0, j = 0; i < this.numBins; ++i, j+=2) {
            this.freqComplexBuffer[j] = reductionSpectrum[i];
            this.freqComplexBuffer[j+1] = reductionSpectrum[i];
        };
        this.fft.completeSpectrum(this.freqComplexBuffer)
        this.fft.inverseTransform(this.timeComplexBuffer, this.freqComplexBuffer);
        this.fft.fromComplexArray(this.timeComplexBuffer, this.irBuffer.getChannelData(0));
        this.node.buffer = this.irBuffer
    }

    connect(destination) { return this.node.connect(destination) }
}

module.exports = { FFTConvFilter }