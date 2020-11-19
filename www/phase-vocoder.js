(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
'use strict';

function FFT(size) {
  this.size = size | 0;
  if (this.size <= 1 || (this.size & (this.size - 1)) !== 0)
    throw new Error('FFT size must be a power of two and bigger than 1');

  this._csize = size << 1;

  // NOTE: Use of `var` is intentional for old V8 versions
  var table = new Array(this.size * 2);
  for (var i = 0; i < table.length; i += 2) {
    const angle = Math.PI * i / this.size;
    table[i] = Math.cos(angle);
    table[i + 1] = -Math.sin(angle);
  }
  this.table = table;

  // Find size's power of two
  var power = 0;
  for (var t = 1; this.size > t; t <<= 1)
    power++;

  // Calculate initial step's width:
  //   * If we are full radix-4 - it is 2x smaller to give inital len=8
  //   * Otherwise it is the same as `power` to give len=4
  this._width = power % 2 === 0 ? power - 1 : power;

  // Pre-compute bit-reversal patterns
  this._bitrev = new Array(1 << this._width);
  for (var j = 0; j < this._bitrev.length; j++) {
    this._bitrev[j] = 0;
    for (var shift = 0; shift < this._width; shift += 2) {
      var revShift = this._width - shift - 2;
      this._bitrev[j] |= ((j >>> shift) & 3) << revShift;
    }
  }

  this._out = null;
  this._data = null;
  this._inv = 0;
}
module.exports = FFT;

FFT.prototype.fromComplexArray = function fromComplexArray(complex, storage) {
  var res = storage || new Array(complex.length >>> 1);
  for (var i = 0; i < complex.length; i += 2)
    res[i >>> 1] = complex[i];
  return res;
};

FFT.prototype.createComplexArray = function createComplexArray() {
  const res = new Array(this._csize);
  for (var i = 0; i < res.length; i++)
    res[i] = 0;
  return res;
};

FFT.prototype.toComplexArray = function toComplexArray(input, storage) {
  var res = storage || this.createComplexArray();
  for (var i = 0; i < res.length; i += 2) {
    res[i] = input[i >>> 1];
    res[i + 1] = 0;
  }
  return res;
};

FFT.prototype.completeSpectrum = function completeSpectrum(spectrum) {
  var size = this._csize;
  var half = size >>> 1;
  for (var i = 2; i < half; i += 2) {
    spectrum[size - i] = spectrum[i];
    spectrum[size - i + 1] = -spectrum[i + 1];
  }
};

FFT.prototype.transform = function transform(out, data) {
  if (out === data)
    throw new Error('Input and output buffers must be different');

  this._out = out;
  this._data = data;
  this._inv = 0;
  this._transform4();
  this._out = null;
  this._data = null;
};

FFT.prototype.realTransform = function realTransform(out, data) {
  if (out === data)
    throw new Error('Input and output buffers must be different');

  this._out = out;
  this._data = data;
  this._inv = 0;
  this._realTransform4();
  this._out = null;
  this._data = null;
};

FFT.prototype.inverseTransform = function inverseTransform(out, data) {
  if (out === data)
    throw new Error('Input and output buffers must be different');

  this._out = out;
  this._data = data;
  this._inv = 1;
  this._transform4();
  for (var i = 0; i < out.length; i++)
    out[i] /= this.size;
  this._out = null;
  this._data = null;
};

// radix-4 implementation
//
// NOTE: Uses of `var` are intentional for older V8 version that do not
// support both `let compound assignments` and `const phi`
FFT.prototype._transform4 = function _transform4() {
  var out = this._out;
  var size = this._csize;

  // Initial step (permute and transform)
  var width = this._width;
  var step = 1 << width;
  var len = (size / step) << 1;

  var outOff;
  var t;
  var bitrev = this._bitrev;
  if (len === 4) {
    for (outOff = 0, t = 0; outOff < size; outOff += len, t++) {
      const off = bitrev[t];
      this._singleTransform2(outOff, off, step);
    }
  } else {
    // len === 8
    for (outOff = 0, t = 0; outOff < size; outOff += len, t++) {
      const off = bitrev[t];
      this._singleTransform4(outOff, off, step);
    }
  }

  // Loop through steps in decreasing order
  var inv = this._inv ? -1 : 1;
  var table = this.table;
  for (step >>= 2; step >= 2; step >>= 2) {
    len = (size / step) << 1;
    var quarterLen = len >>> 2;

    // Loop through offsets in the data
    for (outOff = 0; outOff < size; outOff += len) {
      // Full case
      var limit = outOff + quarterLen;
      for (var i = outOff, k = 0; i < limit; i += 2, k += step) {
        const A = i;
        const B = A + quarterLen;
        const C = B + quarterLen;
        const D = C + quarterLen;

        // Original values
        const Ar = out[A];
        const Ai = out[A + 1];
        const Br = out[B];
        const Bi = out[B + 1];
        const Cr = out[C];
        const Ci = out[C + 1];
        const Dr = out[D];
        const Di = out[D + 1];

        // Middle values
        const MAr = Ar;
        const MAi = Ai;

        const tableBr = table[k];
        const tableBi = inv * table[k + 1];
        const MBr = Br * tableBr - Bi * tableBi;
        const MBi = Br * tableBi + Bi * tableBr;

        const tableCr = table[2 * k];
        const tableCi = inv * table[2 * k + 1];
        const MCr = Cr * tableCr - Ci * tableCi;
        const MCi = Cr * tableCi + Ci * tableCr;

        const tableDr = table[3 * k];
        const tableDi = inv * table[3 * k + 1];
        const MDr = Dr * tableDr - Di * tableDi;
        const MDi = Dr * tableDi + Di * tableDr;

        // Pre-Final values
        const T0r = MAr + MCr;
        const T0i = MAi + MCi;
        const T1r = MAr - MCr;
        const T1i = MAi - MCi;
        const T2r = MBr + MDr;
        const T2i = MBi + MDi;
        const T3r = inv * (MBr - MDr);
        const T3i = inv * (MBi - MDi);

        // Final values
        const FAr = T0r + T2r;
        const FAi = T0i + T2i;

        const FCr = T0r - T2r;
        const FCi = T0i - T2i;

        const FBr = T1r + T3i;
        const FBi = T1i - T3r;

        const FDr = T1r - T3i;
        const FDi = T1i + T3r;

        out[A] = FAr;
        out[A + 1] = FAi;
        out[B] = FBr;
        out[B + 1] = FBi;
        out[C] = FCr;
        out[C + 1] = FCi;
        out[D] = FDr;
        out[D + 1] = FDi;
      }
    }
  }
};

// radix-2 implementation
//
// NOTE: Only called for len=4
FFT.prototype._singleTransform2 = function _singleTransform2(outOff, off,
                                                             step) {
  const out = this._out;
  const data = this._data;

  const evenR = data[off];
  const evenI = data[off + 1];
  const oddR = data[off + step];
  const oddI = data[off + step + 1];

  const leftR = evenR + oddR;
  const leftI = evenI + oddI;
  const rightR = evenR - oddR;
  const rightI = evenI - oddI;

  out[outOff] = leftR;
  out[outOff + 1] = leftI;
  out[outOff + 2] = rightR;
  out[outOff + 3] = rightI;
};

// radix-4
//
// NOTE: Only called for len=8
FFT.prototype._singleTransform4 = function _singleTransform4(outOff, off,
                                                             step) {
  const out = this._out;
  const data = this._data;
  const inv = this._inv ? -1 : 1;
  const step2 = step * 2;
  const step3 = step * 3;

  // Original values
  const Ar = data[off];
  const Ai = data[off + 1];
  const Br = data[off + step];
  const Bi = data[off + step + 1];
  const Cr = data[off + step2];
  const Ci = data[off + step2 + 1];
  const Dr = data[off + step3];
  const Di = data[off + step3 + 1];

  // Pre-Final values
  const T0r = Ar + Cr;
  const T0i = Ai + Ci;
  const T1r = Ar - Cr;
  const T1i = Ai - Ci;
  const T2r = Br + Dr;
  const T2i = Bi + Di;
  const T3r = inv * (Br - Dr);
  const T3i = inv * (Bi - Di);

  // Final values
  const FAr = T0r + T2r;
  const FAi = T0i + T2i;

  const FBr = T1r + T3i;
  const FBi = T1i - T3r;

  const FCr = T0r - T2r;
  const FCi = T0i - T2i;

  const FDr = T1r - T3i;
  const FDi = T1i + T3r;

  out[outOff] = FAr;
  out[outOff + 1] = FAi;
  out[outOff + 2] = FBr;
  out[outOff + 3] = FBi;
  out[outOff + 4] = FCr;
  out[outOff + 5] = FCi;
  out[outOff + 6] = FDr;
  out[outOff + 7] = FDi;
};

// Real input radix-4 implementation
FFT.prototype._realTransform4 = function _realTransform4() {
  var out = this._out;
  var size = this._csize;

  // Initial step (permute and transform)
  var width = this._width;
  var step = 1 << width;
  var len = (size / step) << 1;

  var outOff;
  var t;
  var bitrev = this._bitrev;
  if (len === 4) {
    for (outOff = 0, t = 0; outOff < size; outOff += len, t++) {
      const off = bitrev[t];
      this._singleRealTransform2(outOff, off >>> 1, step >>> 1);
    }
  } else {
    // len === 8
    for (outOff = 0, t = 0; outOff < size; outOff += len, t++) {
      const off = bitrev[t];
      this._singleRealTransform4(outOff, off >>> 1, step >>> 1);
    }
  }

  // Loop through steps in decreasing order
  var inv = this._inv ? -1 : 1;
  var table = this.table;
  for (step >>= 2; step >= 2; step >>= 2) {
    len = (size / step) << 1;
    var halfLen = len >>> 1;
    var quarterLen = halfLen >>> 1;
    var hquarterLen = quarterLen >>> 1;

    // Loop through offsets in the data
    for (outOff = 0; outOff < size; outOff += len) {
      for (var i = 0, k = 0; i <= hquarterLen; i += 2, k += step) {
        var A = outOff + i;
        var B = A + quarterLen;
        var C = B + quarterLen;
        var D = C + quarterLen;

        // Original values
        var Ar = out[A];
        var Ai = out[A + 1];
        var Br = out[B];
        var Bi = out[B + 1];
        var Cr = out[C];
        var Ci = out[C + 1];
        var Dr = out[D];
        var Di = out[D + 1];

        // Middle values
        var MAr = Ar;
        var MAi = Ai;

        var tableBr = table[k];
        var tableBi = inv * table[k + 1];
        var MBr = Br * tableBr - Bi * tableBi;
        var MBi = Br * tableBi + Bi * tableBr;

        var tableCr = table[2 * k];
        var tableCi = inv * table[2 * k + 1];
        var MCr = Cr * tableCr - Ci * tableCi;
        var MCi = Cr * tableCi + Ci * tableCr;

        var tableDr = table[3 * k];
        var tableDi = inv * table[3 * k + 1];
        var MDr = Dr * tableDr - Di * tableDi;
        var MDi = Dr * tableDi + Di * tableDr;

        // Pre-Final values
        var T0r = MAr + MCr;
        var T0i = MAi + MCi;
        var T1r = MAr - MCr;
        var T1i = MAi - MCi;
        var T2r = MBr + MDr;
        var T2i = MBi + MDi;
        var T3r = inv * (MBr - MDr);
        var T3i = inv * (MBi - MDi);

        // Final values
        var FAr = T0r + T2r;
        var FAi = T0i + T2i;

        var FBr = T1r + T3i;
        var FBi = T1i - T3r;

        out[A] = FAr;
        out[A + 1] = FAi;
        out[B] = FBr;
        out[B + 1] = FBi;

        // Output final middle point
        if (i === 0) {
          var FCr = T0r - T2r;
          var FCi = T0i - T2i;
          out[C] = FCr;
          out[C + 1] = FCi;
          continue;
        }

        // Do not overwrite ourselves
        if (i === hquarterLen)
          continue;

        // In the flipped case:
        // MAi = -MAi
        // MBr=-MBi, MBi=-MBr
        // MCr=-MCr
        // MDr=MDi, MDi=MDr
        var ST0r = T1r;
        var ST0i = -T1i;
        var ST1r = T0r;
        var ST1i = -T0i;
        var ST2r = -inv * T3i;
        var ST2i = -inv * T3r;
        var ST3r = -inv * T2i;
        var ST3i = -inv * T2r;

        var SFAr = ST0r + ST2r;
        var SFAi = ST0i + ST2i;

        var SFBr = ST1r + ST3i;
        var SFBi = ST1i - ST3r;

        var SA = outOff + quarterLen - i;
        var SB = outOff + halfLen - i;

        out[SA] = SFAr;
        out[SA + 1] = SFAi;
        out[SB] = SFBr;
        out[SB + 1] = SFBi;
      }
    }
  }
};

// radix-2 implementation
//
// NOTE: Only called for len=4
FFT.prototype._singleRealTransform2 = function _singleRealTransform2(outOff,
                                                                     off,
                                                                     step) {
  const out = this._out;
  const data = this._data;

  const evenR = data[off];
  const oddR = data[off + step];

  const leftR = evenR + oddR;
  const rightR = evenR - oddR;

  out[outOff] = leftR;
  out[outOff + 1] = 0;
  out[outOff + 2] = rightR;
  out[outOff + 3] = 0;
};

// radix-4
//
// NOTE: Only called for len=8
FFT.prototype._singleRealTransform4 = function _singleRealTransform4(outOff,
                                                                     off,
                                                                     step) {
  const out = this._out;
  const data = this._data;
  const inv = this._inv ? -1 : 1;
  const step2 = step * 2;
  const step3 = step * 3;

  // Original values
  const Ar = data[off];
  const Br = data[off + step];
  const Cr = data[off + step2];
  const Dr = data[off + step3];

  // Pre-Final values
  const T0r = Ar + Cr;
  const T1r = Ar - Cr;
  const T2r = Br + Dr;
  const T3r = inv * (Br - Dr);

  // Final values
  const FAr = T0r + T2r;

  const FBr = T1r;
  const FBi = -T3r;

  const FCr = T0r - T2r;

  const FDr = T1r;
  const FDi = T3r;

  out[outOff] = FAr;
  out[outOff + 1] = 0;
  out[outOff + 2] = FBr;
  out[outOff + 3] = FBi;
  out[outOff + 4] = FCr;
  out[outOff + 5] = 0;
  out[outOff + 6] = FDr;
  out[outOff + 7] = FDi;
};

},{}],2:[function(require,module,exports){
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

},{"fft.js":1}],3:[function(require,module,exports){
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


},{"./fft-processor.js":2}]},{},[3])
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJub2RlX21vZHVsZXMvZmZ0LmpzL2xpYi9mZnQuanMiLCJzcmMvd29ya2xldC9mZnQtcHJvY2Vzc29yLmpzIiwic3JjL3dvcmtsZXQvcGhhc2Utdm9jb2Rlci5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNmQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5TkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbigpe2Z1bmN0aW9uIHIoZSxuLHQpe2Z1bmN0aW9uIG8oaSxmKXtpZighbltpXSl7aWYoIWVbaV0pe3ZhciBjPVwiZnVuY3Rpb25cIj09dHlwZW9mIHJlcXVpcmUmJnJlcXVpcmU7aWYoIWYmJmMpcmV0dXJuIGMoaSwhMCk7aWYodSlyZXR1cm4gdShpLCEwKTt2YXIgYT1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK2krXCInXCIpO3Rocm93IGEuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixhfXZhciBwPW5baV09e2V4cG9ydHM6e319O2VbaV1bMF0uY2FsbChwLmV4cG9ydHMsZnVuY3Rpb24ocil7dmFyIG49ZVtpXVsxXVtyXTtyZXR1cm4gbyhufHxyKX0scCxwLmV4cG9ydHMscixlLG4sdCl9cmV0dXJuIG5baV0uZXhwb3J0c31mb3IodmFyIHU9XCJmdW5jdGlvblwiPT10eXBlb2YgcmVxdWlyZSYmcmVxdWlyZSxpPTA7aTx0Lmxlbmd0aDtpKyspbyh0W2ldKTtyZXR1cm4gb31yZXR1cm4gcn0pKCkiLCIndXNlIHN0cmljdCc7XG5cbmZ1bmN0aW9uIEZGVChzaXplKSB7XG4gIHRoaXMuc2l6ZSA9IHNpemUgfCAwO1xuICBpZiAodGhpcy5zaXplIDw9IDEgfHwgKHRoaXMuc2l6ZSAmICh0aGlzLnNpemUgLSAxKSkgIT09IDApXG4gICAgdGhyb3cgbmV3IEVycm9yKCdGRlQgc2l6ZSBtdXN0IGJlIGEgcG93ZXIgb2YgdHdvIGFuZCBiaWdnZXIgdGhhbiAxJyk7XG5cbiAgdGhpcy5fY3NpemUgPSBzaXplIDw8IDE7XG5cbiAgLy8gTk9URTogVXNlIG9mIGB2YXJgIGlzIGludGVudGlvbmFsIGZvciBvbGQgVjggdmVyc2lvbnNcbiAgdmFyIHRhYmxlID0gbmV3IEFycmF5KHRoaXMuc2l6ZSAqIDIpO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IHRhYmxlLmxlbmd0aDsgaSArPSAyKSB7XG4gICAgY29uc3QgYW5nbGUgPSBNYXRoLlBJICogaSAvIHRoaXMuc2l6ZTtcbiAgICB0YWJsZVtpXSA9IE1hdGguY29zKGFuZ2xlKTtcbiAgICB0YWJsZVtpICsgMV0gPSAtTWF0aC5zaW4oYW5nbGUpO1xuICB9XG4gIHRoaXMudGFibGUgPSB0YWJsZTtcblxuICAvLyBGaW5kIHNpemUncyBwb3dlciBvZiB0d29cbiAgdmFyIHBvd2VyID0gMDtcbiAgZm9yICh2YXIgdCA9IDE7IHRoaXMuc2l6ZSA+IHQ7IHQgPDw9IDEpXG4gICAgcG93ZXIrKztcblxuICAvLyBDYWxjdWxhdGUgaW5pdGlhbCBzdGVwJ3Mgd2lkdGg6XG4gIC8vICAgKiBJZiB3ZSBhcmUgZnVsbCByYWRpeC00IC0gaXQgaXMgMnggc21hbGxlciB0byBnaXZlIGluaXRhbCBsZW49OFxuICAvLyAgICogT3RoZXJ3aXNlIGl0IGlzIHRoZSBzYW1lIGFzIGBwb3dlcmAgdG8gZ2l2ZSBsZW49NFxuICB0aGlzLl93aWR0aCA9IHBvd2VyICUgMiA9PT0gMCA/IHBvd2VyIC0gMSA6IHBvd2VyO1xuXG4gIC8vIFByZS1jb21wdXRlIGJpdC1yZXZlcnNhbCBwYXR0ZXJuc1xuICB0aGlzLl9iaXRyZXYgPSBuZXcgQXJyYXkoMSA8PCB0aGlzLl93aWR0aCk7XG4gIGZvciAodmFyIGogPSAwOyBqIDwgdGhpcy5fYml0cmV2Lmxlbmd0aDsgaisrKSB7XG4gICAgdGhpcy5fYml0cmV2W2pdID0gMDtcbiAgICBmb3IgKHZhciBzaGlmdCA9IDA7IHNoaWZ0IDwgdGhpcy5fd2lkdGg7IHNoaWZ0ICs9IDIpIHtcbiAgICAgIHZhciByZXZTaGlmdCA9IHRoaXMuX3dpZHRoIC0gc2hpZnQgLSAyO1xuICAgICAgdGhpcy5fYml0cmV2W2pdIHw9ICgoaiA+Pj4gc2hpZnQpICYgMykgPDwgcmV2U2hpZnQ7XG4gICAgfVxuICB9XG5cbiAgdGhpcy5fb3V0ID0gbnVsbDtcbiAgdGhpcy5fZGF0YSA9IG51bGw7XG4gIHRoaXMuX2ludiA9IDA7XG59XG5tb2R1bGUuZXhwb3J0cyA9IEZGVDtcblxuRkZULnByb3RvdHlwZS5mcm9tQ29tcGxleEFycmF5ID0gZnVuY3Rpb24gZnJvbUNvbXBsZXhBcnJheShjb21wbGV4LCBzdG9yYWdlKSB7XG4gIHZhciByZXMgPSBzdG9yYWdlIHx8IG5ldyBBcnJheShjb21wbGV4Lmxlbmd0aCA+Pj4gMSk7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgY29tcGxleC5sZW5ndGg7IGkgKz0gMilcbiAgICByZXNbaSA+Pj4gMV0gPSBjb21wbGV4W2ldO1xuICByZXR1cm4gcmVzO1xufTtcblxuRkZULnByb3RvdHlwZS5jcmVhdGVDb21wbGV4QXJyYXkgPSBmdW5jdGlvbiBjcmVhdGVDb21wbGV4QXJyYXkoKSB7XG4gIGNvbnN0IHJlcyA9IG5ldyBBcnJheSh0aGlzLl9jc2l6ZSk7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgcmVzLmxlbmd0aDsgaSsrKVxuICAgIHJlc1tpXSA9IDA7XG4gIHJldHVybiByZXM7XG59O1xuXG5GRlQucHJvdG90eXBlLnRvQ29tcGxleEFycmF5ID0gZnVuY3Rpb24gdG9Db21wbGV4QXJyYXkoaW5wdXQsIHN0b3JhZ2UpIHtcbiAgdmFyIHJlcyA9IHN0b3JhZ2UgfHwgdGhpcy5jcmVhdGVDb21wbGV4QXJyYXkoKTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCByZXMubGVuZ3RoOyBpICs9IDIpIHtcbiAgICByZXNbaV0gPSBpbnB1dFtpID4+PiAxXTtcbiAgICByZXNbaSArIDFdID0gMDtcbiAgfVxuICByZXR1cm4gcmVzO1xufTtcblxuRkZULnByb3RvdHlwZS5jb21wbGV0ZVNwZWN0cnVtID0gZnVuY3Rpb24gY29tcGxldGVTcGVjdHJ1bShzcGVjdHJ1bSkge1xuICB2YXIgc2l6ZSA9IHRoaXMuX2NzaXplO1xuICB2YXIgaGFsZiA9IHNpemUgPj4+IDE7XG4gIGZvciAodmFyIGkgPSAyOyBpIDwgaGFsZjsgaSArPSAyKSB7XG4gICAgc3BlY3RydW1bc2l6ZSAtIGldID0gc3BlY3RydW1baV07XG4gICAgc3BlY3RydW1bc2l6ZSAtIGkgKyAxXSA9IC1zcGVjdHJ1bVtpICsgMV07XG4gIH1cbn07XG5cbkZGVC5wcm90b3R5cGUudHJhbnNmb3JtID0gZnVuY3Rpb24gdHJhbnNmb3JtKG91dCwgZGF0YSkge1xuICBpZiAob3V0ID09PSBkYXRhKVxuICAgIHRocm93IG5ldyBFcnJvcignSW5wdXQgYW5kIG91dHB1dCBidWZmZXJzIG11c3QgYmUgZGlmZmVyZW50Jyk7XG5cbiAgdGhpcy5fb3V0ID0gb3V0O1xuICB0aGlzLl9kYXRhID0gZGF0YTtcbiAgdGhpcy5faW52ID0gMDtcbiAgdGhpcy5fdHJhbnNmb3JtNCgpO1xuICB0aGlzLl9vdXQgPSBudWxsO1xuICB0aGlzLl9kYXRhID0gbnVsbDtcbn07XG5cbkZGVC5wcm90b3R5cGUucmVhbFRyYW5zZm9ybSA9IGZ1bmN0aW9uIHJlYWxUcmFuc2Zvcm0ob3V0LCBkYXRhKSB7XG4gIGlmIChvdXQgPT09IGRhdGEpXG4gICAgdGhyb3cgbmV3IEVycm9yKCdJbnB1dCBhbmQgb3V0cHV0IGJ1ZmZlcnMgbXVzdCBiZSBkaWZmZXJlbnQnKTtcblxuICB0aGlzLl9vdXQgPSBvdXQ7XG4gIHRoaXMuX2RhdGEgPSBkYXRhO1xuICB0aGlzLl9pbnYgPSAwO1xuICB0aGlzLl9yZWFsVHJhbnNmb3JtNCgpO1xuICB0aGlzLl9vdXQgPSBudWxsO1xuICB0aGlzLl9kYXRhID0gbnVsbDtcbn07XG5cbkZGVC5wcm90b3R5cGUuaW52ZXJzZVRyYW5zZm9ybSA9IGZ1bmN0aW9uIGludmVyc2VUcmFuc2Zvcm0ob3V0LCBkYXRhKSB7XG4gIGlmIChvdXQgPT09IGRhdGEpXG4gICAgdGhyb3cgbmV3IEVycm9yKCdJbnB1dCBhbmQgb3V0cHV0IGJ1ZmZlcnMgbXVzdCBiZSBkaWZmZXJlbnQnKTtcblxuICB0aGlzLl9vdXQgPSBvdXQ7XG4gIHRoaXMuX2RhdGEgPSBkYXRhO1xuICB0aGlzLl9pbnYgPSAxO1xuICB0aGlzLl90cmFuc2Zvcm00KCk7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgb3V0Lmxlbmd0aDsgaSsrKVxuICAgIG91dFtpXSAvPSB0aGlzLnNpemU7XG4gIHRoaXMuX291dCA9IG51bGw7XG4gIHRoaXMuX2RhdGEgPSBudWxsO1xufTtcblxuLy8gcmFkaXgtNCBpbXBsZW1lbnRhdGlvblxuLy9cbi8vIE5PVEU6IFVzZXMgb2YgYHZhcmAgYXJlIGludGVudGlvbmFsIGZvciBvbGRlciBWOCB2ZXJzaW9uIHRoYXQgZG8gbm90XG4vLyBzdXBwb3J0IGJvdGggYGxldCBjb21wb3VuZCBhc3NpZ25tZW50c2AgYW5kIGBjb25zdCBwaGlgXG5GRlQucHJvdG90eXBlLl90cmFuc2Zvcm00ID0gZnVuY3Rpb24gX3RyYW5zZm9ybTQoKSB7XG4gIHZhciBvdXQgPSB0aGlzLl9vdXQ7XG4gIHZhciBzaXplID0gdGhpcy5fY3NpemU7XG5cbiAgLy8gSW5pdGlhbCBzdGVwIChwZXJtdXRlIGFuZCB0cmFuc2Zvcm0pXG4gIHZhciB3aWR0aCA9IHRoaXMuX3dpZHRoO1xuICB2YXIgc3RlcCA9IDEgPDwgd2lkdGg7XG4gIHZhciBsZW4gPSAoc2l6ZSAvIHN0ZXApIDw8IDE7XG5cbiAgdmFyIG91dE9mZjtcbiAgdmFyIHQ7XG4gIHZhciBiaXRyZXYgPSB0aGlzLl9iaXRyZXY7XG4gIGlmIChsZW4gPT09IDQpIHtcbiAgICBmb3IgKG91dE9mZiA9IDAsIHQgPSAwOyBvdXRPZmYgPCBzaXplOyBvdXRPZmYgKz0gbGVuLCB0KyspIHtcbiAgICAgIGNvbnN0IG9mZiA9IGJpdHJldlt0XTtcbiAgICAgIHRoaXMuX3NpbmdsZVRyYW5zZm9ybTIob3V0T2ZmLCBvZmYsIHN0ZXApO1xuICAgIH1cbiAgfSBlbHNlIHtcbiAgICAvLyBsZW4gPT09IDhcbiAgICBmb3IgKG91dE9mZiA9IDAsIHQgPSAwOyBvdXRPZmYgPCBzaXplOyBvdXRPZmYgKz0gbGVuLCB0KyspIHtcbiAgICAgIGNvbnN0IG9mZiA9IGJpdHJldlt0XTtcbiAgICAgIHRoaXMuX3NpbmdsZVRyYW5zZm9ybTQob3V0T2ZmLCBvZmYsIHN0ZXApO1xuICAgIH1cbiAgfVxuXG4gIC8vIExvb3AgdGhyb3VnaCBzdGVwcyBpbiBkZWNyZWFzaW5nIG9yZGVyXG4gIHZhciBpbnYgPSB0aGlzLl9pbnYgPyAtMSA6IDE7XG4gIHZhciB0YWJsZSA9IHRoaXMudGFibGU7XG4gIGZvciAoc3RlcCA+Pj0gMjsgc3RlcCA+PSAyOyBzdGVwID4+PSAyKSB7XG4gICAgbGVuID0gKHNpemUgLyBzdGVwKSA8PCAxO1xuICAgIHZhciBxdWFydGVyTGVuID0gbGVuID4+PiAyO1xuXG4gICAgLy8gTG9vcCB0aHJvdWdoIG9mZnNldHMgaW4gdGhlIGRhdGFcbiAgICBmb3IgKG91dE9mZiA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4pIHtcbiAgICAgIC8vIEZ1bGwgY2FzZVxuICAgICAgdmFyIGxpbWl0ID0gb3V0T2ZmICsgcXVhcnRlckxlbjtcbiAgICAgIGZvciAodmFyIGkgPSBvdXRPZmYsIGsgPSAwOyBpIDwgbGltaXQ7IGkgKz0gMiwgayArPSBzdGVwKSB7XG4gICAgICAgIGNvbnN0IEEgPSBpO1xuICAgICAgICBjb25zdCBCID0gQSArIHF1YXJ0ZXJMZW47XG4gICAgICAgIGNvbnN0IEMgPSBCICsgcXVhcnRlckxlbjtcbiAgICAgICAgY29uc3QgRCA9IEMgKyBxdWFydGVyTGVuO1xuXG4gICAgICAgIC8vIE9yaWdpbmFsIHZhbHVlc1xuICAgICAgICBjb25zdCBBciA9IG91dFtBXTtcbiAgICAgICAgY29uc3QgQWkgPSBvdXRbQSArIDFdO1xuICAgICAgICBjb25zdCBCciA9IG91dFtCXTtcbiAgICAgICAgY29uc3QgQmkgPSBvdXRbQiArIDFdO1xuICAgICAgICBjb25zdCBDciA9IG91dFtDXTtcbiAgICAgICAgY29uc3QgQ2kgPSBvdXRbQyArIDFdO1xuICAgICAgICBjb25zdCBEciA9IG91dFtEXTtcbiAgICAgICAgY29uc3QgRGkgPSBvdXRbRCArIDFdO1xuXG4gICAgICAgIC8vIE1pZGRsZSB2YWx1ZXNcbiAgICAgICAgY29uc3QgTUFyID0gQXI7XG4gICAgICAgIGNvbnN0IE1BaSA9IEFpO1xuXG4gICAgICAgIGNvbnN0IHRhYmxlQnIgPSB0YWJsZVtrXTtcbiAgICAgICAgY29uc3QgdGFibGVCaSA9IGludiAqIHRhYmxlW2sgKyAxXTtcbiAgICAgICAgY29uc3QgTUJyID0gQnIgKiB0YWJsZUJyIC0gQmkgKiB0YWJsZUJpO1xuICAgICAgICBjb25zdCBNQmkgPSBCciAqIHRhYmxlQmkgKyBCaSAqIHRhYmxlQnI7XG5cbiAgICAgICAgY29uc3QgdGFibGVDciA9IHRhYmxlWzIgKiBrXTtcbiAgICAgICAgY29uc3QgdGFibGVDaSA9IGludiAqIHRhYmxlWzIgKiBrICsgMV07XG4gICAgICAgIGNvbnN0IE1DciA9IENyICogdGFibGVDciAtIENpICogdGFibGVDaTtcbiAgICAgICAgY29uc3QgTUNpID0gQ3IgKiB0YWJsZUNpICsgQ2kgKiB0YWJsZUNyO1xuXG4gICAgICAgIGNvbnN0IHRhYmxlRHIgPSB0YWJsZVszICoga107XG4gICAgICAgIGNvbnN0IHRhYmxlRGkgPSBpbnYgKiB0YWJsZVszICogayArIDFdO1xuICAgICAgICBjb25zdCBNRHIgPSBEciAqIHRhYmxlRHIgLSBEaSAqIHRhYmxlRGk7XG4gICAgICAgIGNvbnN0IE1EaSA9IERyICogdGFibGVEaSArIERpICogdGFibGVEcjtcblxuICAgICAgICAvLyBQcmUtRmluYWwgdmFsdWVzXG4gICAgICAgIGNvbnN0IFQwciA9IE1BciArIE1DcjtcbiAgICAgICAgY29uc3QgVDBpID0gTUFpICsgTUNpO1xuICAgICAgICBjb25zdCBUMXIgPSBNQXIgLSBNQ3I7XG4gICAgICAgIGNvbnN0IFQxaSA9IE1BaSAtIE1DaTtcbiAgICAgICAgY29uc3QgVDJyID0gTUJyICsgTURyO1xuICAgICAgICBjb25zdCBUMmkgPSBNQmkgKyBNRGk7XG4gICAgICAgIGNvbnN0IFQzciA9IGludiAqIChNQnIgLSBNRHIpO1xuICAgICAgICBjb25zdCBUM2kgPSBpbnYgKiAoTUJpIC0gTURpKTtcblxuICAgICAgICAvLyBGaW5hbCB2YWx1ZXNcbiAgICAgICAgY29uc3QgRkFyID0gVDByICsgVDJyO1xuICAgICAgICBjb25zdCBGQWkgPSBUMGkgKyBUMmk7XG5cbiAgICAgICAgY29uc3QgRkNyID0gVDByIC0gVDJyO1xuICAgICAgICBjb25zdCBGQ2kgPSBUMGkgLSBUMmk7XG5cbiAgICAgICAgY29uc3QgRkJyID0gVDFyICsgVDNpO1xuICAgICAgICBjb25zdCBGQmkgPSBUMWkgLSBUM3I7XG5cbiAgICAgICAgY29uc3QgRkRyID0gVDFyIC0gVDNpO1xuICAgICAgICBjb25zdCBGRGkgPSBUMWkgKyBUM3I7XG5cbiAgICAgICAgb3V0W0FdID0gRkFyO1xuICAgICAgICBvdXRbQSArIDFdID0gRkFpO1xuICAgICAgICBvdXRbQl0gPSBGQnI7XG4gICAgICAgIG91dFtCICsgMV0gPSBGQmk7XG4gICAgICAgIG91dFtDXSA9IEZDcjtcbiAgICAgICAgb3V0W0MgKyAxXSA9IEZDaTtcbiAgICAgICAgb3V0W0RdID0gRkRyO1xuICAgICAgICBvdXRbRCArIDFdID0gRkRpO1xuICAgICAgfVxuICAgIH1cbiAgfVxufTtcblxuLy8gcmFkaXgtMiBpbXBsZW1lbnRhdGlvblxuLy9cbi8vIE5PVEU6IE9ubHkgY2FsbGVkIGZvciBsZW49NFxuRkZULnByb3RvdHlwZS5fc2luZ2xlVHJhbnNmb3JtMiA9IGZ1bmN0aW9uIF9zaW5nbGVUcmFuc2Zvcm0yKG91dE9mZiwgb2ZmLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHN0ZXApIHtcbiAgY29uc3Qgb3V0ID0gdGhpcy5fb3V0O1xuICBjb25zdCBkYXRhID0gdGhpcy5fZGF0YTtcblxuICBjb25zdCBldmVuUiA9IGRhdGFbb2ZmXTtcbiAgY29uc3QgZXZlbkkgPSBkYXRhW29mZiArIDFdO1xuICBjb25zdCBvZGRSID0gZGF0YVtvZmYgKyBzdGVwXTtcbiAgY29uc3Qgb2RkSSA9IGRhdGFbb2ZmICsgc3RlcCArIDFdO1xuXG4gIGNvbnN0IGxlZnRSID0gZXZlblIgKyBvZGRSO1xuICBjb25zdCBsZWZ0SSA9IGV2ZW5JICsgb2RkSTtcbiAgY29uc3QgcmlnaHRSID0gZXZlblIgLSBvZGRSO1xuICBjb25zdCByaWdodEkgPSBldmVuSSAtIG9kZEk7XG5cbiAgb3V0W291dE9mZl0gPSBsZWZ0UjtcbiAgb3V0W291dE9mZiArIDFdID0gbGVmdEk7XG4gIG91dFtvdXRPZmYgKyAyXSA9IHJpZ2h0UjtcbiAgb3V0W291dE9mZiArIDNdID0gcmlnaHRJO1xufTtcblxuLy8gcmFkaXgtNFxuLy9cbi8vIE5PVEU6IE9ubHkgY2FsbGVkIGZvciBsZW49OFxuRkZULnByb3RvdHlwZS5fc2luZ2xlVHJhbnNmb3JtNCA9IGZ1bmN0aW9uIF9zaW5nbGVUcmFuc2Zvcm00KG91dE9mZiwgb2ZmLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHN0ZXApIHtcbiAgY29uc3Qgb3V0ID0gdGhpcy5fb3V0O1xuICBjb25zdCBkYXRhID0gdGhpcy5fZGF0YTtcbiAgY29uc3QgaW52ID0gdGhpcy5faW52ID8gLTEgOiAxO1xuICBjb25zdCBzdGVwMiA9IHN0ZXAgKiAyO1xuICBjb25zdCBzdGVwMyA9IHN0ZXAgKiAzO1xuXG4gIC8vIE9yaWdpbmFsIHZhbHVlc1xuICBjb25zdCBBciA9IGRhdGFbb2ZmXTtcbiAgY29uc3QgQWkgPSBkYXRhW29mZiArIDFdO1xuICBjb25zdCBCciA9IGRhdGFbb2ZmICsgc3RlcF07XG4gIGNvbnN0IEJpID0gZGF0YVtvZmYgKyBzdGVwICsgMV07XG4gIGNvbnN0IENyID0gZGF0YVtvZmYgKyBzdGVwMl07XG4gIGNvbnN0IENpID0gZGF0YVtvZmYgKyBzdGVwMiArIDFdO1xuICBjb25zdCBEciA9IGRhdGFbb2ZmICsgc3RlcDNdO1xuICBjb25zdCBEaSA9IGRhdGFbb2ZmICsgc3RlcDMgKyAxXTtcblxuICAvLyBQcmUtRmluYWwgdmFsdWVzXG4gIGNvbnN0IFQwciA9IEFyICsgQ3I7XG4gIGNvbnN0IFQwaSA9IEFpICsgQ2k7XG4gIGNvbnN0IFQxciA9IEFyIC0gQ3I7XG4gIGNvbnN0IFQxaSA9IEFpIC0gQ2k7XG4gIGNvbnN0IFQyciA9IEJyICsgRHI7XG4gIGNvbnN0IFQyaSA9IEJpICsgRGk7XG4gIGNvbnN0IFQzciA9IGludiAqIChCciAtIERyKTtcbiAgY29uc3QgVDNpID0gaW52ICogKEJpIC0gRGkpO1xuXG4gIC8vIEZpbmFsIHZhbHVlc1xuICBjb25zdCBGQXIgPSBUMHIgKyBUMnI7XG4gIGNvbnN0IEZBaSA9IFQwaSArIFQyaTtcblxuICBjb25zdCBGQnIgPSBUMXIgKyBUM2k7XG4gIGNvbnN0IEZCaSA9IFQxaSAtIFQzcjtcblxuICBjb25zdCBGQ3IgPSBUMHIgLSBUMnI7XG4gIGNvbnN0IEZDaSA9IFQwaSAtIFQyaTtcblxuICBjb25zdCBGRHIgPSBUMXIgLSBUM2k7XG4gIGNvbnN0IEZEaSA9IFQxaSArIFQzcjtcblxuICBvdXRbb3V0T2ZmXSA9IEZBcjtcbiAgb3V0W291dE9mZiArIDFdID0gRkFpO1xuICBvdXRbb3V0T2ZmICsgMl0gPSBGQnI7XG4gIG91dFtvdXRPZmYgKyAzXSA9IEZCaTtcbiAgb3V0W291dE9mZiArIDRdID0gRkNyO1xuICBvdXRbb3V0T2ZmICsgNV0gPSBGQ2k7XG4gIG91dFtvdXRPZmYgKyA2XSA9IEZEcjtcbiAgb3V0W291dE9mZiArIDddID0gRkRpO1xufTtcblxuLy8gUmVhbCBpbnB1dCByYWRpeC00IGltcGxlbWVudGF0aW9uXG5GRlQucHJvdG90eXBlLl9yZWFsVHJhbnNmb3JtNCA9IGZ1bmN0aW9uIF9yZWFsVHJhbnNmb3JtNCgpIHtcbiAgdmFyIG91dCA9IHRoaXMuX291dDtcbiAgdmFyIHNpemUgPSB0aGlzLl9jc2l6ZTtcblxuICAvLyBJbml0aWFsIHN0ZXAgKHBlcm11dGUgYW5kIHRyYW5zZm9ybSlcbiAgdmFyIHdpZHRoID0gdGhpcy5fd2lkdGg7XG4gIHZhciBzdGVwID0gMSA8PCB3aWR0aDtcbiAgdmFyIGxlbiA9IChzaXplIC8gc3RlcCkgPDwgMTtcblxuICB2YXIgb3V0T2ZmO1xuICB2YXIgdDtcbiAgdmFyIGJpdHJldiA9IHRoaXMuX2JpdHJldjtcbiAgaWYgKGxlbiA9PT0gNCkge1xuICAgIGZvciAob3V0T2ZmID0gMCwgdCA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4sIHQrKykge1xuICAgICAgY29uc3Qgb2ZmID0gYml0cmV2W3RdO1xuICAgICAgdGhpcy5fc2luZ2xlUmVhbFRyYW5zZm9ybTIob3V0T2ZmLCBvZmYgPj4+IDEsIHN0ZXAgPj4+IDEpO1xuICAgIH1cbiAgfSBlbHNlIHtcbiAgICAvLyBsZW4gPT09IDhcbiAgICBmb3IgKG91dE9mZiA9IDAsIHQgPSAwOyBvdXRPZmYgPCBzaXplOyBvdXRPZmYgKz0gbGVuLCB0KyspIHtcbiAgICAgIGNvbnN0IG9mZiA9IGJpdHJldlt0XTtcbiAgICAgIHRoaXMuX3NpbmdsZVJlYWxUcmFuc2Zvcm00KG91dE9mZiwgb2ZmID4+PiAxLCBzdGVwID4+PiAxKTtcbiAgICB9XG4gIH1cblxuICAvLyBMb29wIHRocm91Z2ggc3RlcHMgaW4gZGVjcmVhc2luZyBvcmRlclxuICB2YXIgaW52ID0gdGhpcy5faW52ID8gLTEgOiAxO1xuICB2YXIgdGFibGUgPSB0aGlzLnRhYmxlO1xuICBmb3IgKHN0ZXAgPj49IDI7IHN0ZXAgPj0gMjsgc3RlcCA+Pj0gMikge1xuICAgIGxlbiA9IChzaXplIC8gc3RlcCkgPDwgMTtcbiAgICB2YXIgaGFsZkxlbiA9IGxlbiA+Pj4gMTtcbiAgICB2YXIgcXVhcnRlckxlbiA9IGhhbGZMZW4gPj4+IDE7XG4gICAgdmFyIGhxdWFydGVyTGVuID0gcXVhcnRlckxlbiA+Pj4gMTtcblxuICAgIC8vIExvb3AgdGhyb3VnaCBvZmZzZXRzIGluIHRoZSBkYXRhXG4gICAgZm9yIChvdXRPZmYgPSAwOyBvdXRPZmYgPCBzaXplOyBvdXRPZmYgKz0gbGVuKSB7XG4gICAgICBmb3IgKHZhciBpID0gMCwgayA9IDA7IGkgPD0gaHF1YXJ0ZXJMZW47IGkgKz0gMiwgayArPSBzdGVwKSB7XG4gICAgICAgIHZhciBBID0gb3V0T2ZmICsgaTtcbiAgICAgICAgdmFyIEIgPSBBICsgcXVhcnRlckxlbjtcbiAgICAgICAgdmFyIEMgPSBCICsgcXVhcnRlckxlbjtcbiAgICAgICAgdmFyIEQgPSBDICsgcXVhcnRlckxlbjtcblxuICAgICAgICAvLyBPcmlnaW5hbCB2YWx1ZXNcbiAgICAgICAgdmFyIEFyID0gb3V0W0FdO1xuICAgICAgICB2YXIgQWkgPSBvdXRbQSArIDFdO1xuICAgICAgICB2YXIgQnIgPSBvdXRbQl07XG4gICAgICAgIHZhciBCaSA9IG91dFtCICsgMV07XG4gICAgICAgIHZhciBDciA9IG91dFtDXTtcbiAgICAgICAgdmFyIENpID0gb3V0W0MgKyAxXTtcbiAgICAgICAgdmFyIERyID0gb3V0W0RdO1xuICAgICAgICB2YXIgRGkgPSBvdXRbRCArIDFdO1xuXG4gICAgICAgIC8vIE1pZGRsZSB2YWx1ZXNcbiAgICAgICAgdmFyIE1BciA9IEFyO1xuICAgICAgICB2YXIgTUFpID0gQWk7XG5cbiAgICAgICAgdmFyIHRhYmxlQnIgPSB0YWJsZVtrXTtcbiAgICAgICAgdmFyIHRhYmxlQmkgPSBpbnYgKiB0YWJsZVtrICsgMV07XG4gICAgICAgIHZhciBNQnIgPSBCciAqIHRhYmxlQnIgLSBCaSAqIHRhYmxlQmk7XG4gICAgICAgIHZhciBNQmkgPSBCciAqIHRhYmxlQmkgKyBCaSAqIHRhYmxlQnI7XG5cbiAgICAgICAgdmFyIHRhYmxlQ3IgPSB0YWJsZVsyICoga107XG4gICAgICAgIHZhciB0YWJsZUNpID0gaW52ICogdGFibGVbMiAqIGsgKyAxXTtcbiAgICAgICAgdmFyIE1DciA9IENyICogdGFibGVDciAtIENpICogdGFibGVDaTtcbiAgICAgICAgdmFyIE1DaSA9IENyICogdGFibGVDaSArIENpICogdGFibGVDcjtcblxuICAgICAgICB2YXIgdGFibGVEciA9IHRhYmxlWzMgKiBrXTtcbiAgICAgICAgdmFyIHRhYmxlRGkgPSBpbnYgKiB0YWJsZVszICogayArIDFdO1xuICAgICAgICB2YXIgTURyID0gRHIgKiB0YWJsZURyIC0gRGkgKiB0YWJsZURpO1xuICAgICAgICB2YXIgTURpID0gRHIgKiB0YWJsZURpICsgRGkgKiB0YWJsZURyO1xuXG4gICAgICAgIC8vIFByZS1GaW5hbCB2YWx1ZXNcbiAgICAgICAgdmFyIFQwciA9IE1BciArIE1DcjtcbiAgICAgICAgdmFyIFQwaSA9IE1BaSArIE1DaTtcbiAgICAgICAgdmFyIFQxciA9IE1BciAtIE1DcjtcbiAgICAgICAgdmFyIFQxaSA9IE1BaSAtIE1DaTtcbiAgICAgICAgdmFyIFQyciA9IE1CciArIE1EcjtcbiAgICAgICAgdmFyIFQyaSA9IE1CaSArIE1EaTtcbiAgICAgICAgdmFyIFQzciA9IGludiAqIChNQnIgLSBNRHIpO1xuICAgICAgICB2YXIgVDNpID0gaW52ICogKE1CaSAtIE1EaSk7XG5cbiAgICAgICAgLy8gRmluYWwgdmFsdWVzXG4gICAgICAgIHZhciBGQXIgPSBUMHIgKyBUMnI7XG4gICAgICAgIHZhciBGQWkgPSBUMGkgKyBUMmk7XG5cbiAgICAgICAgdmFyIEZCciA9IFQxciArIFQzaTtcbiAgICAgICAgdmFyIEZCaSA9IFQxaSAtIFQzcjtcblxuICAgICAgICBvdXRbQV0gPSBGQXI7XG4gICAgICAgIG91dFtBICsgMV0gPSBGQWk7XG4gICAgICAgIG91dFtCXSA9IEZCcjtcbiAgICAgICAgb3V0W0IgKyAxXSA9IEZCaTtcblxuICAgICAgICAvLyBPdXRwdXQgZmluYWwgbWlkZGxlIHBvaW50XG4gICAgICAgIGlmIChpID09PSAwKSB7XG4gICAgICAgICAgdmFyIEZDciA9IFQwciAtIFQycjtcbiAgICAgICAgICB2YXIgRkNpID0gVDBpIC0gVDJpO1xuICAgICAgICAgIG91dFtDXSA9IEZDcjtcbiAgICAgICAgICBvdXRbQyArIDFdID0gRkNpO1xuICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gRG8gbm90IG92ZXJ3cml0ZSBvdXJzZWx2ZXNcbiAgICAgICAgaWYgKGkgPT09IGhxdWFydGVyTGVuKVxuICAgICAgICAgIGNvbnRpbnVlO1xuXG4gICAgICAgIC8vIEluIHRoZSBmbGlwcGVkIGNhc2U6XG4gICAgICAgIC8vIE1BaSA9IC1NQWlcbiAgICAgICAgLy8gTUJyPS1NQmksIE1CaT0tTUJyXG4gICAgICAgIC8vIE1Dcj0tTUNyXG4gICAgICAgIC8vIE1Ecj1NRGksIE1EaT1NRHJcbiAgICAgICAgdmFyIFNUMHIgPSBUMXI7XG4gICAgICAgIHZhciBTVDBpID0gLVQxaTtcbiAgICAgICAgdmFyIFNUMXIgPSBUMHI7XG4gICAgICAgIHZhciBTVDFpID0gLVQwaTtcbiAgICAgICAgdmFyIFNUMnIgPSAtaW52ICogVDNpO1xuICAgICAgICB2YXIgU1QyaSA9IC1pbnYgKiBUM3I7XG4gICAgICAgIHZhciBTVDNyID0gLWludiAqIFQyaTtcbiAgICAgICAgdmFyIFNUM2kgPSAtaW52ICogVDJyO1xuXG4gICAgICAgIHZhciBTRkFyID0gU1QwciArIFNUMnI7XG4gICAgICAgIHZhciBTRkFpID0gU1QwaSArIFNUMmk7XG5cbiAgICAgICAgdmFyIFNGQnIgPSBTVDFyICsgU1QzaTtcbiAgICAgICAgdmFyIFNGQmkgPSBTVDFpIC0gU1QzcjtcblxuICAgICAgICB2YXIgU0EgPSBvdXRPZmYgKyBxdWFydGVyTGVuIC0gaTtcbiAgICAgICAgdmFyIFNCID0gb3V0T2ZmICsgaGFsZkxlbiAtIGk7XG5cbiAgICAgICAgb3V0W1NBXSA9IFNGQXI7XG4gICAgICAgIG91dFtTQSArIDFdID0gU0ZBaTtcbiAgICAgICAgb3V0W1NCXSA9IFNGQnI7XG4gICAgICAgIG91dFtTQiArIDFdID0gU0ZCaTtcbiAgICAgIH1cbiAgICB9XG4gIH1cbn07XG5cbi8vIHJhZGl4LTIgaW1wbGVtZW50YXRpb25cbi8vXG4vLyBOT1RFOiBPbmx5IGNhbGxlZCBmb3IgbGVuPTRcbkZGVC5wcm90b3R5cGUuX3NpbmdsZVJlYWxUcmFuc2Zvcm0yID0gZnVuY3Rpb24gX3NpbmdsZVJlYWxUcmFuc2Zvcm0yKG91dE9mZixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9mZixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHN0ZXApIHtcbiAgY29uc3Qgb3V0ID0gdGhpcy5fb3V0O1xuICBjb25zdCBkYXRhID0gdGhpcy5fZGF0YTtcblxuICBjb25zdCBldmVuUiA9IGRhdGFbb2ZmXTtcbiAgY29uc3Qgb2RkUiA9IGRhdGFbb2ZmICsgc3RlcF07XG5cbiAgY29uc3QgbGVmdFIgPSBldmVuUiArIG9kZFI7XG4gIGNvbnN0IHJpZ2h0UiA9IGV2ZW5SIC0gb2RkUjtcblxuICBvdXRbb3V0T2ZmXSA9IGxlZnRSO1xuICBvdXRbb3V0T2ZmICsgMV0gPSAwO1xuICBvdXRbb3V0T2ZmICsgMl0gPSByaWdodFI7XG4gIG91dFtvdXRPZmYgKyAzXSA9IDA7XG59O1xuXG4vLyByYWRpeC00XG4vL1xuLy8gTk9URTogT25seSBjYWxsZWQgZm9yIGxlbj04XG5GRlQucHJvdG90eXBlLl9zaW5nbGVSZWFsVHJhbnNmb3JtNCA9IGZ1bmN0aW9uIF9zaW5nbGVSZWFsVHJhbnNmb3JtNChvdXRPZmYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBvZmYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdGVwKSB7XG4gIGNvbnN0IG91dCA9IHRoaXMuX291dDtcbiAgY29uc3QgZGF0YSA9IHRoaXMuX2RhdGE7XG4gIGNvbnN0IGludiA9IHRoaXMuX2ludiA/IC0xIDogMTtcbiAgY29uc3Qgc3RlcDIgPSBzdGVwICogMjtcbiAgY29uc3Qgc3RlcDMgPSBzdGVwICogMztcblxuICAvLyBPcmlnaW5hbCB2YWx1ZXNcbiAgY29uc3QgQXIgPSBkYXRhW29mZl07XG4gIGNvbnN0IEJyID0gZGF0YVtvZmYgKyBzdGVwXTtcbiAgY29uc3QgQ3IgPSBkYXRhW29mZiArIHN0ZXAyXTtcbiAgY29uc3QgRHIgPSBkYXRhW29mZiArIHN0ZXAzXTtcblxuICAvLyBQcmUtRmluYWwgdmFsdWVzXG4gIGNvbnN0IFQwciA9IEFyICsgQ3I7XG4gIGNvbnN0IFQxciA9IEFyIC0gQ3I7XG4gIGNvbnN0IFQyciA9IEJyICsgRHI7XG4gIGNvbnN0IFQzciA9IGludiAqIChCciAtIERyKTtcblxuICAvLyBGaW5hbCB2YWx1ZXNcbiAgY29uc3QgRkFyID0gVDByICsgVDJyO1xuXG4gIGNvbnN0IEZCciA9IFQxcjtcbiAgY29uc3QgRkJpID0gLVQzcjtcblxuICBjb25zdCBGQ3IgPSBUMHIgLSBUMnI7XG5cbiAgY29uc3QgRkRyID0gVDFyO1xuICBjb25zdCBGRGkgPSBUM3I7XG5cbiAgb3V0W291dE9mZl0gPSBGQXI7XG4gIG91dFtvdXRPZmYgKyAxXSA9IDA7XG4gIG91dFtvdXRPZmYgKyAyXSA9IEZCcjtcbiAgb3V0W291dE9mZiArIDNdID0gRkJpO1xuICBvdXRbb3V0T2ZmICsgNF0gPSBGQ3I7XG4gIG91dFtvdXRPZmYgKyA1XSA9IDA7XG4gIG91dFtvdXRPZmYgKyA2XSA9IEZEcjtcbiAgb3V0W291dE9mZiArIDddID0gRkRpO1xufTtcbiIsIlwidXNlIHN0cmljdFwiO1xuY29uc3QgRkZUID0gcmVxdWlyZSgnZmZ0LmpzJyk7XG5jb25zdCBXRUJBVURJT19CTE9DS19TSVpFID0gMTI4O1xuXG5mdW5jdGlvbiBnZW5IYW5uV2luZG93KGxlbmd0aCkge1xuICAgIGxldCB3aW4gPSBuZXcgRmxvYXQzMkFycmF5KGxlbmd0aCk7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBsZW5ndGg7IGkrKykge1xuICAgICAgICB3aW5baV0gPSAwLjUgKiAoMSAtIE1hdGguY29zKDIgKiBNYXRoLlBJICogaSAvIGxlbmd0aCkpO1xuICAgIH1cbiAgICByZXR1cm4gd2luO1xufVxuXG4vKiogRkZUIE5vZGUgKi9cbmNsYXNzIEZGVFByb2Nlc3NvciBleHRlbmRzIEF1ZGlvV29ya2xldFByb2Nlc3NvciB7XG4gICAgY29uc3RydWN0b3Iob3B0aW9ucykge1xuICAgICAgICBzdXBlcihvcHRpb25zKTtcblxuICAgICAgICB0aGlzLm5iSW5wdXRzID0gb3B0aW9ucy5udW1iZXJPZklucHV0cztcbiAgICAgICAgdGhpcy5uYk91dHB1dHMgPSBvcHRpb25zLm51bWJlck9mT3V0cHV0cztcblxuICAgICAgICB0aGlzLmJsb2NrU2l6ZSA9IG9wdGlvbnMucHJvY2Vzc29yT3B0aW9ucy5ibG9ja1NpemU7XG4gICAgICAgICAvLyBUT0RPIGZvciBub3csIHRoZSBvbmx5IHN1cHBvcnQgaG9wIHNpemUgaXMgdGhlIHNpemUgb2YgYSB3ZWIgYXVkaW8gYmxvY2tcbiAgICAgICAgdGhpcy5ob3BTaXplID0gV0VCQVVESU9fQkxPQ0tfU0laRTtcblxuICAgICAgICB0aGlzLm5iT3ZlcmxhcHMgPSB0aGlzLmJsb2NrU2l6ZSAvIHRoaXMuaG9wU2l6ZTtcblxuICAgICAgICAvLyBwcmUtYWxsb2NhdGUgaW5wdXQgYnVmZmVycyAod2lsbCBiZSByZWFsbG9jYXRlZCBpZiBuZWVkZWQpXG4gICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzID0gbmV3IEFycmF5KHRoaXMubmJJbnB1dHMpO1xuICAgICAgICB0aGlzLmlucHV0QnVmZmVyc0hlYWQgPSBuZXcgQXJyYXkodGhpcy5uYklucHV0cyk7XG4gICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzVG9TZW5kID0gbmV3IEFycmF5KHRoaXMubmJJbnB1dHMpO1xuICAgICAgICAvLyBkZWZhdWx0IHRvIDEgY2hhbm5lbCBwZXIgaW5wdXQgdW50aWwgd2Uga25vdyBtb3JlXG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5uYklucHV0czsgaSsrKSB7XG4gICAgICAgICAgICB0aGlzLmFsbG9jYXRlSW5wdXRDaGFubmVscyhpLCAxKTtcbiAgICAgICAgfVxuICAgICAgICAvLyBwcmUtYWxsb2NhdGUgaW5wdXQgYnVmZmVycyAod2lsbCBiZSByZWFsbG9jYXRlZCBpZiBuZWVkZWQpXG4gICAgICAgIHRoaXMub3V0cHV0QnVmZmVycyA9IG5ldyBBcnJheSh0aGlzLm5iT3V0cHV0cyk7XG4gICAgICAgIHRoaXMub3V0cHV0QnVmZmVyc1RvUmV0cmlldmUgPSBuZXcgQXJyYXkodGhpcy5uYk91dHB1dHMpO1xuICAgICAgICAvLyBkZWZhdWx0IHRvIDEgY2hhbm5lbCBwZXIgb3V0cHV0IHVudGlsIHdlIGtub3cgbW9yZVxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHRoaXMubmJPdXRwdXRzOyBpKyspIHtcbiAgICAgICAgICAgIHRoaXMuYWxsb2NhdGVPdXRwdXRDaGFubmVscyhpLCAxKTtcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIHByZXBhcmUgRkZUIGFuZCBwcmUtYWxsb2NhdGUgYnVmZmVyc1xuICAgICAgICB0aGlzLmZmdFNpemUgPSB0aGlzLmJsb2NrU2l6ZTtcbiAgICAgICAgdGhpcy5udW1CaW5zID0gdGhpcy5mZnRTaXplIC8gMiArIDE7XG4gICAgICAgIHRoaXMuZmZ0ID0gbmV3IEZGVCh0aGlzLmZmdFNpemUpO1xuICAgICAgICB0aGlzLmZyZXFDb21wbGV4QnVmZmVyID0gdGhpcy5mZnQuY3JlYXRlQ29tcGxleEFycmF5KCk7XG4gICAgICAgIHRoaXMudGltZUNvbXBsZXhCdWZmZXIgPSB0aGlzLmZmdC5jcmVhdGVDb21wbGV4QXJyYXkoKTtcblxuICAgICAgICB0aGlzLmhhbm5XaW5kb3cgPSBnZW5IYW5uV2luZG93KHRoaXMuYmxvY2tTaXplKTtcblxuICAgIH1cblxuICAgIC8qKiBIYW5kbGVzIGR5bmFtaWMgcmVhbGxvY2F0aW9uIG9mIGlucHV0L291dHB1dCBjaGFubmVscyBidWZmZXJcbiAgICAgKGNoYW5uZWwgbnVtYmVycyBtYXkgdmFyeSBkdXJpbmcgbGlmZWN5Y2xlKSAqKi9cbiAgICByZWFsbG9jYXRlQ2hhbm5lbHNJZk5lZWRlZChpbnB1dHMsIG91dHB1dHMpIHtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm5iSW5wdXRzOyBpKyspIHtcbiAgICAgICAgICAgIGxldCBuYkNoYW5uZWxzID0gaW5wdXRzW2ldLmxlbmd0aDtcbiAgICAgICAgICAgIGlmIChuYkNoYW5uZWxzICE9IHRoaXMuaW5wdXRCdWZmZXJzW2ldLmxlbmd0aCkge1xuICAgICAgICAgICAgICAgIHRoaXMuYWxsb2NhdGVJbnB1dENoYW5uZWxzKGksIG5iQ2hhbm5lbHMpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm5iT3V0cHV0czsgaSsrKSB7XG4gICAgICAgICAgICBsZXQgbmJDaGFubmVscyA9IG91dHB1dHNbaV0ubGVuZ3RoO1xuICAgICAgICAgICAgaWYgKG5iQ2hhbm5lbHMgIT0gdGhpcy5vdXRwdXRCdWZmZXJzW2ldLmxlbmd0aCkge1xuICAgICAgICAgICAgICAgIHRoaXMuYWxsb2NhdGVPdXRwdXRDaGFubmVscyhpLCBuYkNoYW5uZWxzKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIGFsbG9jYXRlSW5wdXRDaGFubmVscyhpbnB1dEluZGV4LCBuYkNoYW5uZWxzKSB7XG4gICAgICAgIC8vIGFsbG9jYXRlIGlucHV0IGJ1ZmZlcnNcblxuICAgICAgICB0aGlzLmlucHV0QnVmZmVyc1tpbnB1dEluZGV4XSA9IG5ldyBBcnJheShuYkNoYW5uZWxzKTtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBuYkNoYW5uZWxzOyBpKyspIHtcbiAgICAgICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzW2lucHV0SW5kZXhdW2ldID0gbmV3IEZsb2F0MzJBcnJheSh0aGlzLmJsb2NrU2l6ZSArIFdFQkFVRElPX0JMT0NLX1NJWkUpO1xuICAgICAgICAgICAgdGhpcy5pbnB1dEJ1ZmZlcnNbaW5wdXRJbmRleF1baV0uZmlsbCgwKTtcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIGFsbG9jYXRlIGlucHV0IGJ1ZmZlcnMgdG8gc2VuZCBhbmQgaGVhZCBwb2ludGVycyB0byBjb3B5IGZyb21cbiAgICAgICAgLy8gKGNhbm5vdCBkaXJlY3RseSBzZW5kIGEgcG9pbnRlci9zdWJhcnJheSBiZWNhdXNlIGlucHV0IG1heSBiZSBtb2RpZmllZClcbiAgICAgICAgdGhpcy5pbnB1dEJ1ZmZlcnNIZWFkW2lucHV0SW5kZXhdID0gbmV3IEFycmF5KG5iQ2hhbm5lbHMpO1xuICAgICAgICB0aGlzLmlucHV0QnVmZmVyc1RvU2VuZFtpbnB1dEluZGV4XSA9IG5ldyBBcnJheShuYkNoYW5uZWxzKTtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBuYkNoYW5uZWxzOyBpKyspIHtcbiAgICAgICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzSGVhZFtpbnB1dEluZGV4XVtpXSA9IHRoaXMuaW5wdXRCdWZmZXJzW2lucHV0SW5kZXhdW2ldIC5zdWJhcnJheSgwLCB0aGlzLmJsb2NrU2l6ZSk7XG4gICAgICAgICAgICB0aGlzLmlucHV0QnVmZmVyc1RvU2VuZFtpbnB1dEluZGV4XVtpXSA9IG5ldyBGbG9hdDMyQXJyYXkodGhpcy5ibG9ja1NpemUpO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgYWxsb2NhdGVPdXRwdXRDaGFubmVscyhvdXRwdXRJbmRleCwgbmJDaGFubmVscykge1xuICAgICAgICAvLyBhbGxvY2F0ZSBvdXRwdXQgYnVmZmVyc1xuICAgICAgICB0aGlzLm91dHB1dEJ1ZmZlcnNbb3V0cHV0SW5kZXhdID0gbmV3IEFycmF5KG5iQ2hhbm5lbHMpO1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG5iQ2hhbm5lbHM7IGkrKykge1xuICAgICAgICAgICAgdGhpcy5vdXRwdXRCdWZmZXJzW291dHB1dEluZGV4XVtpXSA9IG5ldyBGbG9hdDMyQXJyYXkodGhpcy5ibG9ja1NpemUpO1xuICAgICAgICAgICAgdGhpcy5vdXRwdXRCdWZmZXJzW291dHB1dEluZGV4XVtpXS5maWxsKDApO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gYWxsb2NhdGUgb3V0cHV0IGJ1ZmZlcnMgdG8gcmV0cmlldmVcbiAgICAgICAgLy8gKGNhbm5vdCBzZW5kIGEgcG9pbnRlci9zdWJhcnJheSBiZWNhdXNlIG5ldyBvdXRwdXQgaGFzIHRvIGJlIGFkZCB0byBleGlzaW5nIG91dHB1dClcbiAgICAgICAgdGhpcy5vdXRwdXRCdWZmZXJzVG9SZXRyaWV2ZVtvdXRwdXRJbmRleF0gPSBuZXcgQXJyYXkobmJDaGFubmVscyk7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbmJDaGFubmVsczsgaSsrKSB7XG4gICAgICAgICAgICB0aGlzLm91dHB1dEJ1ZmZlcnNUb1JldHJpZXZlW291dHB1dEluZGV4XVtpXSA9IG5ldyBGbG9hdDMyQXJyYXkodGhpcy5ibG9ja1NpemUpO1xuICAgICAgICAgICAgdGhpcy5vdXRwdXRCdWZmZXJzVG9SZXRyaWV2ZVtvdXRwdXRJbmRleF1baV0uZmlsbCgwKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKiBSZWFkIG5leHQgd2ViIGF1ZGlvIGJsb2NrIHRvIGlucHV0IGJ1ZmZlcnMgKiovXG4gICAgcmVhZElucHV0cyhpbnB1dHMpIHtcbiAgICAgICAgLy8gd2hlbiBwbGF5YmFjayBpcyBwYXVzZWQsIHdlIG1heSBzdG9wIHJlY2VpdmluZyBuZXcgc2FtcGxlc1xuICAgICAgICBpZiAoaW5wdXRzWzBdLmxlbmd0aCA9PSAwIHx8IGlucHV0c1swXVswXS5sZW5ndGggPT0gMCkge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm5iSW5wdXRzOyBpKyspIHtcbiAgICAgICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IHRoaXMuaW5wdXRCdWZmZXJzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzW2ldW2pdLmZpbGwoMCwgdGhpcy5ibG9ja1NpemUpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuXG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5uYklucHV0czsgaSsrKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IHRoaXMuaW5wdXRCdWZmZXJzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgbGV0IHdlYkF1ZGlvQmxvY2sgPSBpbnB1dHNbaV1bal07XG4gICAgICAgICAgICAgICAgdGhpcy5pbnB1dEJ1ZmZlcnNbaV1bal0uc2V0KHdlYkF1ZGlvQmxvY2ssIHRoaXMuYmxvY2tTaXplKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIC8qKiBXcml0ZSBuZXh0IHdlYiBhdWRpbyBibG9jayBmcm9tIG91dHB1dCBidWZmZXJzICoqL1xuICAgIHdyaXRlT3V0cHV0cyhvdXRwdXRzKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5uYklucHV0czsgaSsrKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IHRoaXMuaW5wdXRCdWZmZXJzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgbGV0IHdlYkF1ZGlvQmxvY2sgPSB0aGlzLm91dHB1dEJ1ZmZlcnNbaV1bal0uc3ViYXJyYXkoMCwgV0VCQVVESU9fQkxPQ0tfU0laRSk7XG4gICAgICAgICAgICAgICAgb3V0cHV0c1tpXVtqXS5zZXQod2ViQXVkaW9CbG9jayk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKiogU2hpZnQgbGVmdCBjb250ZW50IG9mIGlucHV0IGJ1ZmZlcnMgdG8gcmVjZWl2ZSBuZXcgd2ViIGF1ZGlvIGJsb2NrICoqL1xuICAgIHNoaWZ0SW5wdXRCdWZmZXJzKCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHRoaXMubmJJbnB1dHM7IGkrKykge1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCB0aGlzLmlucHV0QnVmZmVyc1tpXS5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIHRoaXMuaW5wdXRCdWZmZXJzW2ldW2pdLmNvcHlXaXRoaW4oMCwgV0VCQVVESU9fQkxPQ0tfU0laRSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKiogU2hpZnQgbGVmdCBjb250ZW50IG9mIG91dHB1dCBidWZmZXJzIHRvIHJlY2VpdmUgbmV3IHdlYiBhdWRpbyBibG9jayAqKi9cbiAgICBzaGlmdE91dHB1dEJ1ZmZlcnMoKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5uYk91dHB1dHM7IGkrKykge1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCB0aGlzLm91dHB1dEJ1ZmZlcnNbaV0ubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICB0aGlzLm91dHB1dEJ1ZmZlcnNbaV1bal0uY29weVdpdGhpbigwLCBXRUJBVURJT19CTE9DS19TSVpFKTtcbiAgICAgICAgICAgICAgICB0aGlzLm91dHB1dEJ1ZmZlcnNbaV1bal0uc3ViYXJyYXkodGhpcy5ibG9ja1NpemUgLSBXRUJBVURJT19CTE9DS19TSVpFKS5maWxsKDApO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfVxuXG4gICAgLyoqIENvcHkgY29udGVudHMgb2YgaW5wdXQgYnVmZmVycyB0byBidWZmZXIgYWN0dWFsbHkgc2VudCB0byBwcm9jZXNzICoqL1xuICAgIHByZXBhcmVJbnB1dEJ1ZmZlcnNUb1NlbmQoKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5uYklucHV0czsgaSsrKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMDsgaiA8IHRoaXMuaW5wdXRCdWZmZXJzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgdGhpcy5pbnB1dEJ1ZmZlcnNUb1NlbmRbaV1bal0uc2V0KHRoaXMuaW5wdXRCdWZmZXJzSGVhZFtpXVtqXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICAvKiogQWRkIGNvbnRlbnRzIG9mIG91dHB1dCBidWZmZXJzIGp1c3QgcHJvY2Vzc2VkIHRvIG91dHB1dCBidWZmZXJzICoqL1xuICAgIGhhbmRsZU91dHB1dEJ1ZmZlcnNUb1JldHJpZXZlKCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHRoaXMubmJPdXRwdXRzOyBpKyspIHtcbiAgICAgICAgICAgIGZvciAodmFyIGogPSAwOyBqIDwgdGhpcy5vdXRwdXRCdWZmZXJzW2ldLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgayA9IDA7IGsgPCB0aGlzLmJsb2NrU2l6ZTsgaysrKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMub3V0cHV0QnVmZmVyc1tpXVtqXVtrXSArPSB0aGlzLm91dHB1dEJ1ZmZlcnNUb1JldHJpZXZlW2ldW2pdW2tdIC8gdGhpcy5uYk92ZXJsYXBzO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIHByb2Nlc3MoaW5wdXRzLCBvdXRwdXRzLCBwYXJhbXMpIHtcbiAgICAgICAgdGhpcy5yZWFsbG9jYXRlQ2hhbm5lbHNJZk5lZWRlZChpbnB1dHMsIG91dHB1dHMpO1xuXG4gICAgICAgIHRoaXMucmVhZElucHV0cyhpbnB1dHMpO1xuICAgICAgICB0aGlzLnNoaWZ0SW5wdXRCdWZmZXJzKCk7XG4gICAgICAgIHRoaXMucHJlcGFyZUlucHV0QnVmZmVyc1RvU2VuZCgpXG5cbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm5iSW5wdXRzOyBpKyspIHtcbiAgICAgICAgICAgIGZvciAodmFyIGogPSAwOyBqIDwgdGhpcy5pbnB1dEJ1ZmZlcnNUb1NlbmRbaV0ubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICAvLyBiaWcgYXNzdW1wdGlvbiBoZXJlOiBvdXRwdXQgaXMgc3ltZXRyaWMgdG8gaW5wdXRcbiAgICAgICAgICAgICAgICB2YXIgaW5wdXQgPSB0aGlzLmlucHV0QnVmZmVyc1RvU2VuZFtpXVtqXTtcbiAgICAgICAgICAgICAgICB2YXIgb3V0cHV0ID0gdGhpcy5vdXRwdXRCdWZmZXJzVG9SZXRyaWV2ZVtpXVtqXTtcblxuICAgICAgICAgICAgICAgIHRoaXMuYXBwbHlIYW5uV2luZG93KGlucHV0KTtcbiAgICAgICAgICAgICAgICB0aGlzLmZmdC5yZWFsVHJhbnNmb3JtKHRoaXMuZnJlcUNvbXBsZXhCdWZmZXIsIGlucHV0KTtcbiAgICAgICAgICAgICAgICB0aGlzLnByb2Nlc3NGRlQodGhpcy5mcmVxQ29tcGxleEJ1ZmZlciwgb3V0cHV0LCBwYXJhbXMpO1xuXG4gICAgICAgICAgICAgICAgdGhpcy5mZnQuY29tcGxldGVTcGVjdHJ1bSh0aGlzLmZyZXFDb21wbGV4QnVmZmVyKTtcbiAgICAgICAgICAgICAgICB0aGlzLmZmdC5pbnZlcnNlVHJhbnNmb3JtKHRoaXMudGltZUNvbXBsZXhCdWZmZXIsIHRoaXMuZnJlcUNvbXBsZXhCdWZmZXIpO1xuICAgICAgICAgICAgICAgIHRoaXMuZmZ0LmZyb21Db21wbGV4QXJyYXkodGhpcy50aW1lQ29tcGxleEJ1ZmZlciwgb3V0cHV0KTtcblxuICAgICAgICAgICAgICAgIHRoaXMuYXBwbHlIYW5uV2luZG93KG91dHB1dCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuXG4gICAgICAgIHRoaXMuaGFuZGxlT3V0cHV0QnVmZmVyc1RvUmV0cmlldmUoKTtcbiAgICAgICAgdGhpcy53cml0ZU91dHB1dHMob3V0cHV0cyk7XG4gICAgICAgIHRoaXMuc2hpZnRPdXRwdXRCdWZmZXJzKCk7XG5cbiAgICAgICAgcmV0dXJuIHRydWU7XG4gICAgfVxuXG4gICAgcHJvY2Vzc0ZGVChjb21wbGV4U3BlY3RydW0sIG91dHB1dHMsIHBhcmFtcykge1xuICAgICAgICBjb25zb2xlLmFzc2VydChmYWxzZSwgXCJOb3Qgb3ZlcnJpZGVuXCIpO1xuICAgIH1cblxuICAgIC8qKiBBcHBseSBIYW5uIHdpbmRvdyBpbi1wbGFjZSAqL1xuICAgIGFwcGx5SGFubldpbmRvdyhpbnB1dCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IHRoaXMuYmxvY2tTaXplOyBpKyspIHtcbiAgICAgICAgICAgIGlucHV0W2ldID0gaW5wdXRbaV0gKiB0aGlzLmhhbm5XaW5kb3dbaV07XG4gICAgICAgIH1cbiAgICB9XG59XG5cbm1vZHVsZS5leHBvcnRzID0gRkZUUHJvY2Vzc29yO1xuIiwiXCJ1c2Ugc3RyaWN0XCI7XG5cbmNvbnN0IEZGVFByb2Nlc3NvciA9IHJlcXVpcmUoJy4vZmZ0LXByb2Nlc3Nvci5qcycpO1xuXG5jb25zdCBCVUZGRVJFRF9CTE9DS19TSVpFID0gNTEyO1xuY29uc3QgREJfTUlOID0gLTE4MDtcblxuY2xhc3MgTVNEIHtcblxuICAgIGNvbnN0cnVjdG9yKG51bUJpbnMsIGhpc3RvcnlTaXplKXtcbiAgICAgICAgdGhpcy5tc2QgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmYgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeSA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyAqIGhpc3RvcnlTaXplKTtcbiAgICAgICAgdGhpcy5sYXN0TWFnbml0dWRlcyA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubGFzdE1hZ0RpZmYgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmTm9ybWFsaXplID0gMSAvIGhpc3RvcnlTaXplO1xuICAgIH1cblxuICAgIGFkZEZvckJpbihiaW5JbmRleCwgbWFnRGIsIHNtb290aGluZykge1xuICAgICAgICB0aGlzLm1hZ0RpZmZbYmluSW5kZXhdID0gbWFnRGIgLSB0aGlzLmxhc3RNYWduaXR1ZGVzW2JpbkluZGV4XTtcbiAgICAgICAgY29uc3QgbWFnRGlmZkRpZmYgPSBNYXRoLnBvdyh0aGlzLm1hZ0RpZmZbYmluSW5kZXhdIC0gdGhpcy5sYXN0TWFnRGlmZltiaW5JbmRleF0sMik7XG4gICAgICAgIHRoaXMubXNkW2JpbkluZGV4XSArPSAoMS1zbW9vdGhpbmcpICogbWFnRGlmZkRpZmYgKyBzbW9vdGhpbmcgKiB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeVtiaW5JbmRleF0gO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeVtiaW5JbmRleF0gPSBtYWdEaWZmRGlmZjtcblxuICAgICAgICB0aGlzLmxhc3RNYWdEaWZmW2JpbkluZGV4XSA9IHRoaXMubWFnRGlmZltiaW5JbmRleF07XG4gICAgICAgIHRoaXMubGFzdE1hZ25pdHVkZXNbYmluSW5kZXhdID0gdGhpcy5tYWdEYltiaW5JbmRleF07XG4gICAgfVxuXG59XG5cbmNsYXNzIE1hZ25pdHVkZXNIaXN0b3J5IHtcblxuICAgIGNvbnN0cnVjdG9yKG51bUJpbnMsIGhpc3RvcnlTaXplKSB7XG4gICAgICAgIHRoaXMubnVtQmlucyA9IG51bUJpbnM7XG4gICAgICAgIHRoaXMuck51bUJpbnMgPSAxLjAgLyBudW1CaW5zO1xuICAgICAgICB0aGlzLmhpc3RvcnlTaXplID0gaGlzdG9yeVNpemU7XG4gICAgICAgIHRoaXMubWFnRGIgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ1NjYWxlID0gMiAvIHRoaXMubnVtQmlucztcblxuICAgICAgICB0aGlzLm1zZCA9IG5ldyBNU0QodGhpcy5udW1CaW5zLCBoaXN0b3J5U2l6ZSk7XG5cbiAgICAgICAgdGhpcy5zbW9vdGhpbmcgPSAwLjk5O1xuICAgICAgICB0aGlzLnNtb290aGVkTWFnRGIgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLmF2ZXJhZ2VEYiA9IDA7XG4gICAgICAgIHRoaXMuYXZnVGhyID0gMC41O1xuICAgICAgICB0aGlzLnNtb290aGVkTWF4ID0gREJfTUlOO1xuXG4gICAgICAgIHRoaXMubWF4UGVha3MgPSAxMDtcbiAgICAgICAgdGhpcy5uYlBlYWtzID0gMDtcbiAgICAgICAgdGhpcy5taW5QZWFrVGhyID0gLTQwO1xuICAgICAgICB0aGlzLnBlYWtUaHIgPSAwO1xuICAgICAgICB0aGlzLnBlYWtJbmRleGVzID0gbmV3IEludDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMucGVha0hpc3RvcnlTaXplID0gMTAwO1xuICAgICAgICB0aGlzLnBlYWtIaXN0b3J5ID0gbmV3IEZsb2F0MzJBcnJheShudW1CaW5zICogdGhpcy5wZWFrSGlzdG9yeVNpemUpO1xuICAgICAgICB0aGlzLnBlYWtQZXJzaXN0ZW5jZSA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubWluUGVha1BlcnNpc3RlbmNlID0gdGhpcy5wZWFrSGlzdG9yeVNpemUgKiAwLjc1O1xuICAgICAgICB0aGlzLnBlYWtXaWR0aCA9IDExO1xuXG4gICAgICAgIHRoaXMuZmJJbmRleGVzID0gbmV3IEludDMyQXJyYXkodGhpcy5udW1CaW5zKTtcbiAgICAgICAgdGhpcy5uYkZiID0gMDtcbiAgICAgICAgdGhpcy5mYkhpc3RvcnkgPSBuZXcgSW50MzJBcnJheSgyNSk7XG4gICAgICAgIHRoaXMuZmJIaXN0b3J5UG9zID0gMDtcblxuICAgICAgICB0aGlzLm1hZ1JlZHVjdGlvbnMgPSBuZXcgRmxvYXQzMkFycmF5KHRoaXMubnVtQmlucyk7XG5cbiAgICAgICAgY29uc29sZS5sb2coXCJGRlQgY29uc3RydWN0ZWRcIilcbiAgICB9XG5cbiAgICBzaGlmdEhpc3RvcnkoKSB7XG4gICAgICAgIC8vY29uc3QgbGFzdE1lbUJsb2NrID0gIHRoaXMubnVtQmlucyAqICh0aGlzLmhpc3RvcnlTaXplIC0gMSk7XG4gICAgICAgIGNvbnN0IGxhc3RQZWFrTWVtQmxvY2sgPSB0aGlzLm51bUJpbnMgKiAodGhpcy5wZWFrSGlzdG9yeVNpemUgLSAxKTtcbiAgICAgICAgZm9yKGxldCBiaW4gPSAwOyBiaW4gPCB0aGlzLm51bUJpbnM7ICsrYmluKSB7XG4gICAgICAgICAgICAvLyB0aGlzLm1zZC5tc2RbYmluXSAtPSB0aGlzLm1zZC5tYWdEaWZmRGlmZkhpc3RvcnlbbGFzdE1lbUJsb2NrICsgYmluXTtcbiAgICAgICAgICAgIHRoaXMucGVha1BlcnNpc3RlbmNlW2Jpbl0gLT0gdGhpcy5wZWFrSGlzdG9yeVtsYXN0UGVha01lbUJsb2NrICsgYmluXTsgIFxuICAgICAgICAgICAgaWYodGhpcy5tYWdSZWR1Y3Rpb25zW2Jpbl08MCkgdGhpcy5tYWdSZWR1Y3Rpb25zW2Jpbl0gKz0gMC4xXG4gICAgICAgIH1cblxuICAgICAgICB0aGlzLnBlYWtUaHIgPSB0aGlzLmF2ZXJhZ2VEYiAqIHRoaXMuck51bUJpbnM7XG4gICAgICAgIHRoaXMucGVha1RociArPSAodGhpcy5zbW9vdGhlZE1heCAtIHRoaXMucGVha1RocikgKiB0aGlzLmF2Z1RocjtcbiAgICAgICAgaWYodGhpcy5wZWFrVGhyIDwgdGhpcy5taW5QZWFrVGhyKSB0aGlzLnBlYWtUaHIgPSB0aGlzLm1pblBlYWtUaHJcbiAgICAgICAgdGhpcy5hdmVyYWdlRGIgPSAwO1xuICAgICAgICB0aGlzLnNtb290aGVkTWF4ID0gREJfTUlOO1xuICAgICAgICB0aGlzLm5iUGVha3MgPSAwO1xuICAgICAgICAvLyB0aGlzLm1zZC5tYWdEaWZmRGlmZkhpc3RvcnkuY29weVdpdGhpbih0aGlzLm51bUJpbnMsMCk7XG4gICAgICAgIHRoaXMucGVha0hpc3RvcnkuY29weVdpdGhpbih0aGlzLm51bUJpbnMsMCk7XG5cbiAgICB9XG5cblxuICAgIGFkZEZvckJpbihiaW5JbmRleCwgbWFnKSB7XG4gICAgICAgIGxldCBtYWdEYiA9ICBtYWcgPT0gMCA/IERCX01JTiA6IDIwICogTWF0aC5sb2cxMChtYWcgKiB0aGlzLm1hZ1NjYWxlKTtcbiAgICAgICAgdGhpcy5zbW9vdGhlZE1hZ0RiW2JpbkluZGV4XSA9ICgxLXRoaXMuc21vb3RoaW5nKSAqIG1hZ0RiICsgdGhpcy5zbW9vdGhpbmcgKiB0aGlzLnNtb290aGVkTWFnRGJbYmluSW5kZXhdO1xuICAgICAgICB0aGlzLm1hZ0RiW2JpbkluZGV4XSA9IG1hZ0RiO1xuICAgICAgICAvLyB0aGlzLm1zZC5hZGRGb3JCaW4oYmluSW5kZXgsIG1hZ0RiLCB0aGlzLnNtb290aGluZylcblxuICAgICAgICB0aGlzLmF2ZXJhZ2VEYiArPSB0aGlzLnNtb290aGVkTWFnRGJbYmluSW5kZXhdO1xuICAgICAgICBpZih0aGlzLnNtb290aGVkTWFnRGJbYmluSW5kZXhdID4gdGhpcy5zbW9vdGhlZE1heCkgdGhpcy5zbW9vdGhlZE1heCA9IHRoaXMuc21vb3RoZWRNYWdEYltiaW5JbmRleF1cbiAgICAgICAgdGhpcy5wZWFrSGlzdG9yeVtiaW5JbmRleF0gPSAwO1xuICAgICAgICAvL3RoaXMuY2hlY2tQZWFrKGJpbkluZGV4KTtcbiAgICB9XG5cbiAgICBmaW5kUGVha3MoKXtcbiAgICAgICAgY29uc3QgbWFncyA9IHRoaXMuc21vb3RoZWRNYWdEYjtcbiAgICAgICAgY29uc3QgdGhyID0gdGhpcy5wZWFrVGhyO1xuICAgICAgICBsZXQgYmxvY2tNYXggPSAtSW5maW5pdHk7XG4gICAgICAgIGxldCBibG9ja01heEJpbiA9IC0xO1xuICAgICAgICBsZXQgbF9iaW4gPSAwOyBsZXQgcl9iaW4gPSAwOyBcblxuICAgICAgICBjb25zdCBjaGVja0Jsb2NrID0gKCk9PntcbiAgICAgICAgICAgIHJfYmluID0gbF9iaW47XG4gICAgICAgICAgICBibG9ja01heCA9IC1JbmZpbml0eTtcbiAgICAgICAgICAgIGJsb2NrTWF4QmluID0gLTE7XG4gICAgICAgICAgICB3aGlsZShyX2JpbiA8IGxfYmluICsgdGhpcy5wZWFrV2lkdGggKiAyKSB7XG4gICAgICAgICAgICAgICAgaWYobWFnc1tyX2Jpbl0gPj0gdGhyICYmIG1hZ3Nbcl9iaW5dID4gYmxvY2tNYXgpIHsgYmxvY2tNYXggPSBtYWdzW3JfYmluXTsgYmxvY2tNYXhCaW4gPSByX2JpbiB9XG4gICAgICAgICAgICAgICAgcl9iaW4rKztcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGNoZWNrQmxvY2soKTtcbiAgICAgICAgbGV0IGJpbiA9IHJfYmluIC0gdGhpcy5wZWFrV2lkdGg7XG4gICAgICAgIGlmKGJsb2NrTWF4QmluID09IGJpbikgdGhpcy5hZGRQZWFrKGJpbik7XG4gICAgICAgIHJfYmluKys7IGxfYmluKys7IGJpbisrO1xuXG4gICAgICAgIHdoaWxlKHJfYmluIDwgdGhpcy5udW1CaW5zKXtcbiAgICAgICAgICAgIGlmKG1hZ3Nbcl9iaW5dID49IHRociAmJiBtYWdzW3JfYmluXSA+IGJsb2NrTWF4KSB7IFxuICAgICAgICAgICAgICAgIGJsb2NrTWF4ID0gbWFnc1tyX2Jpbl07IGJsb2NrTWF4QmluID0gcl9iaW5cbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgaWYoYmxvY2tNYXhCaW4gPiAwICYmIGJsb2NrTWF4QmluIDwgbF9iaW4pIGNoZWNrQmxvY2soKTtcbiAgICAgICAgICAgICAgICBpZihibG9ja01heEJpbiA9PSBiaW4pIHRoaXMuYWRkUGVhayhiaW4pO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcl9iaW4rKzsgbF9iaW4rKzsgYmluKys7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBhZGRQZWFrKGJpbkluZGV4KXtcbiAgICAgICAgdGhpcy5wZWFrSW5kZXhlc1t0aGlzLm5iUGVha3NdID0gYmluSW5kZXg7XG4gICAgICAgIHRoaXMucGVha0hpc3RvcnlbYmluSW5kZXhdID0gMTtcbiAgICAgICAgdGhpcy5wZWFrUGVyc2lzdGVuY2VbYmluSW5kZXhdKys7XG4gICAgICAgIHRoaXMubmJQZWFrcysrO1xuICAgIH1cbiAgICBcbiAgICBmaW5kRmVlZGJhY2tDYW5kaWRhdGVzKCkge1xuICAgICAgICB0aGlzLm5iRmIgPSAwO1xuICAgICAgICBsZXQgbWluTXNkQmluID0gLTE7XG4gICAgICAgIGxldCBtaW5Nc2QgPSBJbmZpbml0eTtcblxuICAgICAgICAvL2NvbnNvbGUubG9nKHRoaXMucGVha0hpc3RvcnkpXG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgdGhpcy5uYlBlYWtzOyArK2kpIHtcbiAgICAgICAgICAgIGNvbnN0IGJpbiA9IHRoaXMucGVha0luZGV4ZXNbaV07XG4gICAgICAgICAgICBpZih0aGlzLnBlYWtQZXJzaXN0ZW5jZVtiaW5dID49IHRoaXMubWluUGVha1BlcnNpc3RlbmNlKSB7XG4gICAgICAgICAgICAvKmlmKHRoaXMubXNkW2Jpbl0gPCBtaW5Nc2QgJiYgdGhpcy5zbW9vdGhlZE1hZ0RiW2Jpbi0xXSA8IHRoaXMuc21vb3RoZWRNYWdEYltiaW5dICYmIHRoaXMuc21vb3RoZWRNYWdEYltiaW4rMV0gPCB0aGlzLnNtb290aGVkTWFnRGJbYmluXSkge1xuICAgICAgICAgICAgICAgIG1pbk1zZCA9IHRoaXMubXNkW2Jpbl07XG4gICAgICAgICAgICAgICAgbWluTXNkQmluID0gYmluO1xuICAgICAgICAgICAgfSovXG4gICAgICAgICAgICB0aGlzLmZiSW5kZXhlc1t0aGlzLm5iRmIrK10gPSBiaW47XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB0aGlzLmNvcnJlY3RNYWduaXR1ZGUoYmluKVxuXG4gICAgICAgIH1cblxuICAgICAgICBcbiAgICAgICAgLyppZihtaW5Nc2RCaW4gPj0gMCAmJiBtaW5Nc2QgPD0gMC4wMSkge1xuICAgICAgICAgICAgdGhpcy5mYkhpc3RvcnlbdGhpcy5mYkhpc3RvcnlQb3MrK10gPSBtaW5Nc2RCaW47XG4gICAgICAgICAgICBpZih0aGlzLmZiSGlzdG9yeS5yZWR1Y2UoKHQsdik9PnY9PW1pbk1zZEJpbj90KzE6dCwwKSA+PSB0aGlzLmZiSGlzdG9yeS5sZW5ndGggKiAwLjUpIHtcbiAgICAgICAgICAgICAgICB0aGlzLmZiSW5kZXhlc1t0aGlzLm5iRmIrK10gPSBtaW5Nc2RCaW47XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0aGlzLmZiSGlzdG9yeVt0aGlzLmZiSGlzdG9yeVBvcysrXSA9IC0xO1xuICAgICAgICB9XG4gICAgICAgIHRoaXMuZmJIaXN0b3J5LmZpbHRlcigodixpLGEpPT52ID49IDAgJiYgYS5pbmRleE9mKHYpPT09aSkuZm9yRWFjaChiaW49PntcbiAgICAgICAgICAgIHRoaXMuZmJJbmRleGVzW3RoaXMubmJGYl0gPSBiaW47XG4gICAgICAgICAgICB0aGlzLm5iRmIrKztcbiAgICAgICAgfSlcbiAgICAgICAgaWYodGhpcy5mYkhpc3RvcnlQb3MgPiB0aGlzLmZiSGlzdG9yeS5sZW5ndGgpIHRoaXMuZmJIaXN0b3J5UG9zID0gMDtcbiAgICAgICAgKi9cbiAgICB9XG5cbiAgICBjb3JyZWN0TWFnbml0dWRlKGJpbikge1xuICAgICAgICB0aGlzLm1hZ1JlZHVjdGlvbnNbYmluXSArPSAodGhpcy5wZWFrVGhyIC0gdGhpcy5tYWdEYltiaW5dKSAvIDJcbiAgICB9XG5cbn1cblxuY2xhc3MgUGhhc2VWb2NvZGVyUHJvY2Vzc29yIGV4dGVuZHMgRkZUUHJvY2Vzc29yIHtcblxuICAgIGNvbnN0cnVjdG9yKG9wdGlvbnMpIHtcbiAgICAgICAgb3B0aW9ucy5wcm9jZXNzb3JPcHRpb25zID0ge1xuICAgICAgICAgICAgYmxvY2tTaXplOiBvcHRpb25zLnByb2Nlc3Nvck9wdGlvbnMuYmxvY2tTaXplIHx8IEJVRkZFUkVEX0JMT0NLX1NJWkUsXG4gICAgICAgIH07XG4gICAgICAgIHN1cGVyKG9wdGlvbnMpO1xuICAgICAgICB0aGlzLm1hZ0hpc3RvcnkgPSBuZXcgTWFnbml0dWRlc0hpc3RvcnkodGhpcy5udW1CaW5zLCAxMSk7XG4gICAgICAgIHRoaXMucGhhc2VzID0gbmV3IEZsb2F0MzJBcnJheSh0aGlzLm51bUJpbnMpXG4gICAgICAgIHRoaXMubWFnbml0dWRlcyA9IG5ldyBGbG9hdDMyQXJyYXkodGhpcy5udW1CaW5zKVxuXG4gICAgICAgIHRoaXMubGlzdGVuVG9Ob2RlTWVzc2FnZXMoKTtcbiAgICB9XG5cbiAgICBsaXN0ZW5Ub05vZGVNZXNzYWdlcygpIHtcbiAgICAgICAgdGhpcy5yZXNwb25kZXJzID0ge1xuICAgICAgICAgICAgZ2V0U3BlY3RydW06ICgpPT50aGlzLnNlbmRNYWduaXR1ZGVzKCksXG4gICAgICAgICAgICBnZXRTbG9wZTogKCk9PnRoaXMuc2VuZFNsb3BlKCksXG4gICAgICAgICAgICBnZXRQZWFrczogKCk9PnRoaXMuc2VuZFBlYWtzKCksXG4gICAgICAgICAgICBnZXRIaXN0bzogKCk9PnRoaXMuc2VuZEhpc3RvKCksXG4gICAgICAgICAgICBnZXRGYkNhbmRpZGF0ZXM6ICgpPT50aGlzLnNlbmRGYigpLFxuICAgICAgICAgICAgZ2V0TWF4RGI6ICgpPT50aGlzLnNlbmRNYXhEYigpLFxuICAgICAgICAgICAgZ2V0TWFnUmVkdWN0aW9uczogKCk9PnRoaXMuc2VuZE1hZ1JlZHVjdGlvbnMoKVxuICAgICAgICB9O1xuICAgICAgICB0aGlzLnBvcnQub25tZXNzYWdlID0gZXZlbnQgPT4ge1xuICAgICAgICAgICAgY29uc3QgbWV0aG9kID0gZXZlbnQuZGF0YTtcbiAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKFwicmVxdWVzdGVkXCIsIG1ldGhvZClcbiAgICAgICAgICAgIGNvbnN0IHJlc3BvbmRlciA9IHRoaXMucmVzcG9uZGVyc1ttZXRob2RdO1xuICAgICAgICAgICAgaWYocmVzcG9uZGVyKSByZXNwb25kZXIoKTtcbiAgICAgICAgfTtcbiAgICAgICAgY29uc29sZS5sb2coXCJQcm9jZXNzb3IgbGlzdGVuaW5nXCIpXG4gICAgfTtcbiAgICBzZW5kTWFnbml0dWRlcyhkYXRhKSB7IHRoaXMucG9ydC5wb3N0TWVzc2FnZShbXCJzcGVjdHJ1bVwiLCB0aGlzLm1hZ0hpc3Rvcnkuc21vb3RoZWRNYWdEYl0pIH1cbiAgICBzZW5kU2xvcGUoKSB7IHRoaXMucG9ydC5wb3N0TWVzc2FnZShbXCJzbG9wZVwiLCB0aGlzLm1hZ0hpc3RvcnkubXNkLm1zZF0pIH1cbiAgICBzZW5kSGlzdG8oKSB7IHRoaXMucG9ydC5wb3N0TWVzc2FnZShbXCJoaXN0b1wiLCB0aGlzLm1hZ0hpc3RvcnkubXNkLm1hZ0RpZmZEaWZmSGlzdG9yeV0pIH1cbiAgICBzZW5kUGVha3MoKSB7ICB0aGlzLnBvcnQucG9zdE1lc3NhZ2UoW1wicGVha3NcIiwgdGhpcy5tYWdIaXN0b3J5Lm5iUGVha3MsIHRoaXMubWFnSGlzdG9yeS5wZWFrSW5kZXhlc10pIH1cbiAgICBzZW5kRmIoKSB7ICB0aGlzLnBvcnQucG9zdE1lc3NhZ2UoW1wiZmJcIix0aGlzLm1hZ0hpc3RvcnkubmJGYiwgdGhpcy5tYWdIaXN0b3J5LmZiSW5kZXhlc10pIH1cbiAgICBzZW5kTWF4RGIoKSB7ICB0aGlzLnBvcnQucG9zdE1lc3NhZ2UoW1wibWF4RGJcIix0aGlzLm1hZ0hpc3Rvcnkuc21vb3RoZWRNYXhdKSB9XG4gICAgc2VuZE1hZ1JlZHVjdGlvbnMoKSB7IHRoaXMucG9ydC5wb3N0TWVzc2FnZShbXCJyZWR1Y3Rpb25zXCIsIHRoaXMubWFnSGlzdG9yeS5tYWdSZWR1Y3Rpb25zXSkgfVxuXG4gICAgcHJvY2Vzc0ZGVChjb21wbGV4U3BlY3RydW0sIG91dHB1dCwgcGFyYW1ldGVycykge1xuICAgICAgICB0aGlzLm1hZ0hpc3Rvcnkuc2hpZnRIaXN0b3J5KCk7XG4gICAgICAgIHRoaXMuY29tcHV0ZU1hZ25pdHVkZXMoY29tcGxleFNwZWN0cnVtKTtcbiAgICAgICAgdGhpcy5tYWdIaXN0b3J5LmZpbmRQZWFrcygpO1xuICAgICAgICB0aGlzLm1hZ0hpc3RvcnkuZmluZEZlZWRiYWNrQ2FuZGlkYXRlcygpO1xuICAgICAgICBsZXQgaiA9IDA7XG4gICAgICAgIGZvcihsZXQgaT0wOyBpPHRoaXMubnVtQmluczsgaSsrLCBqKz0yKSB7XG4gICAgICAgICAgICBjb25zdCBjb3JyZWN0aW9uID0gMTAgKiogKHRoaXMubWFnSGlzdG9yeS5tYWdSZWR1Y3Rpb25zW2ldICogMC4wNSk7XG4gICAgICAgICAgICB0aGlzLmZyZXFDb21wbGV4QnVmZmVyW2pdICo9IGNvcnJlY3Rpb247XG4gICAgICAgICAgICB0aGlzLmZyZXFDb21wbGV4QnVmZmVyW2orMV0gKj0gY29ycmVjdGlvbjtcbiAgICAgICAgfVxuXG4gICAgICAgIFxuICAgIH1cblxuICAgIC8qKiBDb21wdXRlIHNxdWFyZWQgbWFnbml0dWRlcyBmb3IgcGVhayBmaW5kaW5nICoqL1xuICAgIGNvbXB1dGVNYWduaXR1ZGVzKGZmdCkge1xuICAgICAgICB2YXIgaSA9IDAsIGogPSAwO1xuICAgICAgICB3aGlsZSAoaSA8IHRoaXMubnVtQmlucykge1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBsZXQgcmVhbCA9IGZmdFtqXTtcbiAgICAgICAgICAgIGxldCBpbWFnID0gZmZ0W2ogKyAxXTtcbiAgICAgICAgICAgIHRoaXMucGhhc2VzW2ldID0gaW1hZy9yZWFsO1xuICAgICAgICAgICAgdGhpcy5tYWdIaXN0b3J5LmFkZEZvckJpbihpLCBNYXRoLnNxcnQocmVhbCAqKiAyICsgaW1hZyAqKiAyKSlcbiAgICAgICAgICAgIGkrPTE7XG4gICAgICAgICAgICBqKz0yO1xuICAgICAgICB9XG4gICAgfVxufVxuXG5yZWdpc3RlclByb2Nlc3NvcihcInBoYXNlLXZvY29kZXItcHJvY2Vzc29yXCIsIFBoYXNlVm9jb2RlclByb2Nlc3Nvcik7XG5cbiJdfQ==
