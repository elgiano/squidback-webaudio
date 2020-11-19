"use strict";
const FFT = require('fft.js');
const WEBAUDIO_BLOCK_SIZE = 128;

function genHannWindow(length) {
    let win = new Float32Array(length);
    for (var i = 0; i < length; i++) {
        win[i] = 0.5 * (1 - Math.cos(2 * Math.PI * i / length));
    }
    return win;
}

/** FFT Node */
class FFTProcessor extends AudioWorkletProcessor {
    constructor(options) {
        super(options);

        this.nbInputs = options.numberOfInputs;
        this.nbOutputs = options.numberOfOutputs;

        this.blockSize = options.processorOptions.blockSize;
         // TODO for now, the only support hop size is the size of a web audio block
        this.hopSize = WEBAUDIO_BLOCK_SIZE;

        this.nbOverlaps = this.blockSize / this.hopSize;

        // pre-allocate input buffers (will be reallocated if needed)
        this.inputBuffers = new Array(this.nbInputs);
        this.inputBuffersHead = new Array(this.nbInputs);
        this.inputBuffersToSend = new Array(this.nbInputs);
        // default to 1 channel per input until we know more
        for (var i = 0; i < this.nbInputs; i++) {
            this.allocateInputChannels(i, 1);
        }
        // pre-allocate input buffers (will be reallocated if needed)
        this.outputBuffers = new Array(this.nbOutputs);
        this.outputBuffersToRetrieve = new Array(this.nbOutputs);
        // default to 1 channel per output until we know more
        for (var i = 0; i < this.nbOutputs; i++) {
            this.allocateOutputChannels(i, 1);
        }

        // prepare FFT and pre-allocate buffers
        this.fftSize = this.blockSize;
        this.numBins = this.fftSize / 2 + 1;
        this.fft = new FFT(this.fftSize);
        this.freqComplexBuffer = this.fft.createComplexArray();
        this.timeComplexBuffer = this.fft.createComplexArray();

        this.hannWindow = genHannWindow(this.blockSize);

    }

    /** Handles dynamic reallocation of input/output channels buffer
     (channel numbers may vary during lifecycle) **/
    reallocateChannelsIfNeeded(inputs, outputs) {
        for (var i = 0; i < this.nbInputs; i++) {
            let nbChannels = inputs[i].length;
            if (nbChannels != this.inputBuffers[i].length) {
                this.allocateInputChannels(i, nbChannels);
            }
        }

        for (var i = 0; i < this.nbOutputs; i++) {
            let nbChannels = outputs[i].length;
            if (nbChannels != this.outputBuffers[i].length) {
                this.allocateOutputChannels(i, nbChannels);
            }
        }
    }

    allocateInputChannels(inputIndex, nbChannels) {
        // allocate input buffers

        this.inputBuffers[inputIndex] = new Array(nbChannels);
        for (var i = 0; i < nbChannels; i++) {
            this.inputBuffers[inputIndex][i] = new Float32Array(this.blockSize + WEBAUDIO_BLOCK_SIZE);
            this.inputBuffers[inputIndex][i].fill(0);
        }

        // allocate input buffers to send and head pointers to copy from
        // (cannot directly send a pointer/subarray because input may be modified)
        this.inputBuffersHead[inputIndex] = new Array(nbChannels);
        this.inputBuffersToSend[inputIndex] = new Array(nbChannels);
        for (var i = 0; i < nbChannels; i++) {
            this.inputBuffersHead[inputIndex][i] = this.inputBuffers[inputIndex][i] .subarray(0, this.blockSize);
            this.inputBuffersToSend[inputIndex][i] = new Float32Array(this.blockSize);
        }
    }

    allocateOutputChannels(outputIndex, nbChannels) {
        // allocate output buffers
        this.outputBuffers[outputIndex] = new Array(nbChannels);
        for (var i = 0; i < nbChannels; i++) {
            this.outputBuffers[outputIndex][i] = new Float32Array(this.blockSize);
            this.outputBuffers[outputIndex][i].fill(0);
        }

        // allocate output buffers to retrieve
        // (cannot send a pointer/subarray because new output has to be add to exising output)
        this.outputBuffersToRetrieve[outputIndex] = new Array(nbChannels);
        for (var i = 0; i < nbChannels; i++) {
            this.outputBuffersToRetrieve[outputIndex][i] = new Float32Array(this.blockSize);
            this.outputBuffersToRetrieve[outputIndex][i].fill(0);
        }
    }

    /** Read next web audio block to input buffers **/
    readInputs(inputs) {
        // when playback is paused, we may stop receiving new samples
        if (inputs[0].length == 0 || inputs[0][0].length == 0) {
            for (var i = 0; i < this.nbInputs; i++) {
                for (var j = 0; j < this.inputBuffers[i].length; j++) {
                    this.inputBuffers[i][j].fill(0, this.blockSize);
                }
            }
            return;
        }

        for (var i = 0; i < this.nbInputs; i++) {
            for (var j = 0; j < this.inputBuffers[i].length; j++) {
                let webAudioBlock = inputs[i][j];
                this.inputBuffers[i][j].set(webAudioBlock, this.blockSize);
            }
        }
    }

    /** Write next web audio block from output buffers **/
    writeOutputs(outputs) {
        for (var i = 0; i < this.nbInputs; i++) {
            for (var j = 0; j < this.inputBuffers[i].length; j++) {
                let webAudioBlock = this.outputBuffers[i][j].subarray(0, WEBAUDIO_BLOCK_SIZE);
                outputs[i][j].set(webAudioBlock);
            }
        }
    }

    /** Shift left content of input buffers to receive new web audio block **/
    shiftInputBuffers() {
        for (var i = 0; i < this.nbInputs; i++) {
            for (var j = 0; j < this.inputBuffers[i].length; j++) {
                this.inputBuffers[i][j].copyWithin(0, WEBAUDIO_BLOCK_SIZE);
            }
        }
    }

    /** Shift left content of output buffers to receive new web audio block **/
    shiftOutputBuffers() {
        for (var i = 0; i < this.nbOutputs; i++) {
            for (var j = 0; j < this.outputBuffers[i].length; j++) {
                this.outputBuffers[i][j].copyWithin(0, WEBAUDIO_BLOCK_SIZE);
                this.outputBuffers[i][j].subarray(this.blockSize - WEBAUDIO_BLOCK_SIZE).fill(0);
            }
        }
    }

    /** Copy contents of input buffers to buffer actually sent to process **/
    prepareInputBuffersToSend() {
        for (var i = 0; i < this.nbInputs; i++) {
            for (var j = 0; j < this.inputBuffers[i].length; j++) {
                this.inputBuffersToSend[i][j].set(this.inputBuffersHead[i][j]);
            }
        }
    }

    /** Add contents of output buffers just processed to output buffers **/
    handleOutputBuffersToRetrieve() {
        for (var i = 0; i < this.nbOutputs; i++) {
            for (var j = 0; j < this.outputBuffers[i].length; j++) {
                for (var k = 0; k < this.blockSize; k++) {
                    this.outputBuffers[i][j][k] += this.outputBuffersToRetrieve[i][j][k] / this.nbOverlaps;
                }
            }
        }
    }

    process(inputs, outputs, params) {
        this.reallocateChannelsIfNeeded(inputs, outputs);

        this.readInputs(inputs);
        this.shiftInputBuffers();
        this.prepareInputBuffersToSend()

        for (var i = 0; i < this.nbInputs; i++) {
            for (var j = 0; j < this.inputBuffersToSend[i].length; j++) {
                // big assumption here: output is symetric to input
                var input = this.inputBuffersToSend[i][j];
                var output = this.outputBuffersToRetrieve[i][j];

                this.applyHannWindow(input);
                this.fft.realTransform(this.freqComplexBuffer, input);
                this.processFFT(this.freqComplexBuffer, output, params);

                this.fft.completeSpectrum(this.freqComplexBuffer);
                this.fft.inverseTransform(this.timeComplexBuffer, this.freqComplexBuffer);
                this.fft.fromComplexArray(this.timeComplexBuffer, output);

                this.applyHannWindow(output);
            }
        }


        this.handleOutputBuffersToRetrieve();
        this.writeOutputs(outputs);
        this.shiftOutputBuffers();

        return true;
    }

    processFFT(complexSpectrum, outputs, params) {
        console.assert(false, "Not overriden");
    }

    /** Apply Hann window in-place */
    applyHannWindow(input) {
        for (var i = 0; i < this.blockSize; i++) {
            input[i] = input[i] * this.hannWindow[i];
        }
    }
}

module.exports = FFTProcessor;
