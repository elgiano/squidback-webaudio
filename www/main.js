(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
module.exports = class Spline {
  constructor(xs, ys) {
    this.xs = xs;
    this.ys = ys;
    this.ks = this.getNaturalKs(new Float64Array(this.xs.length));
  }

  getNaturalKs(ks) {
    const n = this.xs.length - 1;
    const A = zerosMat(n + 1, n + 2);

    for (
      let i = 1;
      i < n;
      i++ // rows
    ) {
      A[i][i - 1] = 1 / (this.xs[i] - this.xs[i - 1]);
      A[i][i] =
        2 *
        (1 / (this.xs[i] - this.xs[i - 1]) + 1 / (this.xs[i + 1] - this.xs[i]));
      A[i][i + 1] = 1 / (this.xs[i + 1] - this.xs[i]);
      A[i][n + 1] =
        3 *
        ((this.ys[i] - this.ys[i - 1]) /
          ((this.xs[i] - this.xs[i - 1]) * (this.xs[i] - this.xs[i - 1])) +
          (this.ys[i + 1] - this.ys[i]) /
            ((this.xs[i + 1] - this.xs[i]) * (this.xs[i + 1] - this.xs[i])));
    }

    A[0][0] = 2 / (this.xs[1] - this.xs[0]);
    A[0][1] = 1 / (this.xs[1] - this.xs[0]);
    A[0][n + 1] =
      (3 * (this.ys[1] - this.ys[0])) /
      ((this.xs[1] - this.xs[0]) * (this.xs[1] - this.xs[0]));

    A[n][n - 1] = 1 / (this.xs[n] - this.xs[n - 1]);
    A[n][n] = 2 / (this.xs[n] - this.xs[n - 1]);
    A[n][n + 1] =
      (3 * (this.ys[n] - this.ys[n - 1])) /
      ((this.xs[n] - this.xs[n - 1]) * (this.xs[n] - this.xs[n - 1]));

    return solve(A, ks);
  }

  /**
   * inspired by https://stackoverflow.com/a/40850313/4417327
   */
  getIndexBefore(target) {
    let low = 0;
    let high = this.xs.length;
    let mid = 0;
    while (low < high) {
      mid = Math.floor((low + high) / 2);
      if (this.xs[mid] < target && mid !== low) {
        low = mid;
      } else if (this.xs[mid] >= target && mid !== high) {
        high = mid;
      } else {
        high = low;
      }
    }
    return low + 1;
  }

  at(x) {
    let i = this.getIndexBefore(x);
    const t = (x - this.xs[i - 1]) / (this.xs[i] - this.xs[i - 1]);
    const a =
      this.ks[i - 1] * (this.xs[i] - this.xs[i - 1]) -
      (this.ys[i] - this.ys[i - 1]);
    const b =
      -this.ks[i] * (this.xs[i] - this.xs[i - 1]) +
      (this.ys[i] - this.ys[i - 1]);
    const q =
      (1 - t) * this.ys[i - 1] +
      t * this.ys[i] +
      t * (1 - t) * (a * (1 - t) + b * t);
    return q;
  }
};

function solve(A, ks) {
  const m = A.length;
  let h = 0;
  let k = 0;
  while (h < m && k <= m) {
    let i_max = 0;
    let max = -Infinity;
    for (let i = h; i < m; i++) {
      const v = Math.abs(A[i][k]);
      if (v > max) {
        i_max = i;
        max = v;
      }
    }

    if (A[i_max][k] === 0) {
      k++;
    } else {
      swapRows(A, h, i_max);
      for (let i = h + 1; i < m; i++) {
        const f = A[i][k] / A[h][k];
        A[i][k] = 0;
        for (let j = k + 1; j <= m; j++) A[i][j] -= A[h][j] * f;
      }
      h++;
      k++;
    }
  }

  for (
    let i = m - 1;
    i >= 0;
    i-- // rows = columns
  ) {
    var v = 0;
    if (A[i][i]) {
      v = A[i][m] / A[i][i];
    }
    ks[i] = v;
    for (
      let j = i - 1;
      j >= 0;
      j-- // rows
    ) {
      A[j][m] -= A[j][i] * v;
      A[j][i] = 0;
    }
  }
  return ks;
}

function zerosMat(r, c) {
  const A = [];
  for (let i = 0; i < r; i++) A.push(new Float64Array(c));
  return A;
}

function swapRows(m, k, l) {
  let p = m[k];
  m[k] = m[l];
  m[l] = p;
}

},{}],2:[function(require,module,exports){
(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (factory((global.d3_peaks = global.d3_peaks || {})));
}(this, function (exports) { 'use strict';

  /**
   * See https://en.wikipedia.org/wiki/Mexican_hat_wavelet
   */
  function ricker() {
    var σ = 1;
    
    var ricker = function(t) {
      var t2 = t*t,
          variance = σ*σ;
      
      var C = 2.0 / ( Math.sqrt(3 * σ) * (Math.pow(Math.PI, 0.25)) );
      var norm = (1.0 - (t2)/(variance));
      var gauss = Math.exp( -(t2) / (2*variance) );
      
      return C*norm*gauss;
    }
    
    ricker.std = function(_) {
      return arguments.length ? (σ = _, ricker) : σ;
    }
    
    /**
     * Range of points to sample from the wavelet. [-reach, reach]
     */
    ricker.reach = function() {
      return 5 * σ;
    }
    
    return ricker;
  };

  function convolve() {
    var kernel = ricker();
    
    /**
     * y[n] = Sum_k{x[k] * h[n-k]}
     * y: output
     * x: input
     * h: smoother
     */
    var convolve = function(signal) {
      var size = signal.length,
          n = -1,
          convolution = new Array(size);
          
      while (++n < size) {
        var y = 0;
        
        var box = boundingBox(n, kernel.reach(), 0, size - 1);
        box.forEach(function(δ) {
          var k = n + δ;
          y += signal[k] * kernel(δ);
        });
        convolution[n] = y;
      }
      
      return convolution;
    };
    
    convolve.kernel = function(_) {
      return arguments.length ? (kernel = _, convolve) : kernel;
    }
    
    function range(reach) {
      reach = +reach;
      var i = -1,
          n = 2*reach + 1,
          range = new Array(n);
      while(++i < n) {
        range[i] = (-reach) + i;
      }
      return range;
    }
    
    function boundingBox(n, reach, lo, hi) {
      for (var i = 1; i <= reach; i++) {
        var left  = n - i,
            right = n + i;
        if (left >= lo && right <= hi) continue;
        return range(i - 1);
      }
      return range(reach);
    }
    
    return convolve;
  };

  function isLocalMaxima(arr, index) {
    var current = arr[index],
        left = arr[index - 1],
        right = arr[index + 1];
        
    if (left !== undefined && right !== undefined) {
      if (current > left && current > right) { return true; }
      else if (current >= left && current > right) { return true; }
      else if (current > left && current >= right) { return true; }
    }
    else if (left !== undefined && current > left) { return true; }
    else if (right !== undefined && current > right) { return true; }
    
    return false;
  }

  /**
   * @param {arr} row in the CWT matrix.
   * @return Array of indices with relative maximas.
   */
  function maximas(arr) {
    var maximas = [];
    arr.forEach(function(value, index) {
      if (isLocalMaxima(arr, index)) maximas.push({x: index, y: value});
    });
    return maximas;
  };

  function nearestNeighbor(line, maximas, window) {
    var cache = {};
    maximas.forEach(function(d) {
      cache[d.x] = d.y;
    });
    
    var point = line.top();
    for (var i = 0; i <= window; i++) {
      var left = point.x + i;
      var right = point.x - i;
      
      if ( (left in cache) && (right in cache) ) {
        if (cache[left] > cache[right]) {
          return left;
        }
        return right;
      }
      else if (left in cache) {
        return left;
      }
      else if (right in cache) {
        return right;
      }
    }
    return null;
  }

  function percentile(arr, perc) {
    var length = arr.length;
    var index = Math.min(length - 1, Math.ceil(perc * length));
    
    arr.sort(function(a, b) { return a - b; });
    return arr[index];
  }

  function Point(x, y, width) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.snr = undefined;
  }

  Point.prototype.SNR = function(conv) {
    var smoothingFactor = 0.00001;
    var signal = this.y;
    
    var lowerBound = Math.max(0, this.x - this.width);
    var upperBound = Math.min(conv.length, this.x + this.width + 1);
    var neighbors = conv.slice(lowerBound, upperBound);
    var noise = percentile(neighbors, 0.95);
    
    signal += smoothingFactor;
    noise += smoothingFactor;
    this.snr = signal/noise;
    return this.snr;
  }

  Point.prototype.serialize = function() {
    return {index: this.x, width: this.width, snr: this.snr};
  }

  function RidgeLine() {
    this.points = [];
    this.gap = 0;
  }

  /**
   * If the point is valid append it to the ridgeline, and reset the gap.
   * Otherwise, increment the gap and do nothing.
   * 
   * @param {point} Point object.
   */
  RidgeLine.prototype.add = function(point) {
    if (point === null || point === undefined) {
      this.gap += 1;
      return;
    } else {
      this.points.push(point);
      this.gap = 0;
    }
  }

  /**
   * @return {Point} Last point added into the ridgeline.
   */
  RidgeLine.prototype.top = function() {
    return this.points[this.points.length - 1];
  }

  /**
   * @return {number} Length of points on the ridgeline.
   */
  RidgeLine.prototype.length = function() {
    return this.points.length;
  }

  /**
   * @return {boolean} True if the gap in the line is above a threshold. False otherwise.
   */
  RidgeLine.prototype.isDisconnected = function (threshold) {
    return this.gap > threshold;
  }

  /**
   * @param {Array} Smallest scale in the convolution matrix
   */
  RidgeLine.prototype.SNR = function(conv) {
    var maxSnr = Number.NEGATIVE_INFINITY;
    this.points.forEach(function(point) {
      var snr = point.SNR(conv);
      if (snr > maxSnr) maxSnr = snr;
    });
    return maxSnr;
  }

  function findPeaks() {
    var kernel = ricker,
        gapThreshold = 1,
        minLineLength = 1,
        minSNR = 1.0,
        widths = [1];
    
    var findPeaks = function(signal) {
      var M = CWT(signal);
      
      var ridgeLines = initializeRidgeLines(M);
      ridgeLines = connectRidgeLines(M, ridgeLines);
      ridgeLines = filterRidgeLines(signal, ridgeLines);
      
      return peaks(signal, ridgeLines);
    };
    
    /**
     * Smoothing function.
     */
    findPeaks.kernel = function(_) {
      return arguments.length ? (kernel = _, findPeaks) : kernel;
    }
    
    /**
     * Expected widths of the peaks.
     */
    findPeaks.widths = function(_) {
      _.sort(function(a, b) { return a - b; });
      return arguments.length ? (widths = _, findPeaks) : widths;
    }
    
    /**
     * Number of gaps that we allow in the ridge lines.
     */
    findPeaks.gapThreshold = function(_) {
      return arguments.length ? (gapThreshold = _, findPeaks) : gapThreshold;
    }
    
    /**
     * Minimum ridge line length.
     */
    findPeaks.minLineLength = function(_) {
      return arguments.length ? (minLineLength = _, findPeaks) : minLineLength;
    }
    
    /**
     * Minimum signal to noise ratio for the peaks.
     */
    findPeaks.minSNR = function(_) {
      return arguments.length ? (minSNR = _, findPeaks) : minSNR;
    }
    
    var CWT = function(signal) {
      var M = new Array(widths.length);
      widths.forEach(function(width, i) {
        var smoother = kernel()
          .std(width);
        var transform = convolve()
          .kernel(smoother);
        
        var convolution = transform(signal);
        M[i] = convolution;
      });
      return M;
    }
    
    
    var initializeRidgeLines = function(M) {
      var n = widths.length;
      var locals = maximas(M[n - 1], widths[n - 1]);
      var ridgeLines = [];
      locals.forEach(function(d) {
        var point = new Point(d.x, d.y, widths[n - 1]);
        var line = new RidgeLine();
        line.add(point);
        ridgeLines.push(line);
      });
      return ridgeLines;
    }
    
    var connectRidgeLines = function(M, ridgeLines) {
      var n = widths.length;
      for (var row = n - 2; row >= 0; row--) {
        var locals = maximas(M[row], widths[row]);
        var addedLocals = [];
        
        // Find nearest neighbor at next scale and add to the line
        ridgeLines.forEach(function(line, i) {
          var x = nearestNeighbor(line, locals, widths[row]);
          line.add(x === null ? null : new Point(x, M[row][x], widths[row]));
          
          if (x !== null) {
            addedLocals.push(x);
          }
        });
        
        // Remove lines that has exceeded the gap threshold
        ridgeLines = ridgeLines.filter(function(line) {
          return !line.isDisconnected(gapThreshold);
        });
        
        // Add all the unitialized ridge lines
        locals.forEach(function(d) {
          if (addedLocals.indexOf(d.x) !== -1) return;
          
          var point = new Point(d.x, d.y, widths[row]);
          var ridgeLine = new RidgeLine();
          ridgeLine.add(point);
          ridgeLines.push(ridgeLine);
        });
      }
      return ridgeLines;
    }
    
    var filterRidgeLines = function(signal, ridgeLines) {
      var smoother = kernel()
          .std(1.0);
      var transform = convolve()
        .kernel(smoother);
      var convolution = transform(signal);
        
      ridgeLines = ridgeLines.filter(function(line) {
        var snr = line.SNR(convolution);
        return (snr >= minSNR) && (line.length() >= minLineLength);
      });
      return ridgeLines
    }
    
    /**
     * Pick the point with the highest y value within that range.
     */
    var peaks = function(signal, ridgeLines) {
      var peaks = ridgeLines.map(function(line) {
        var points = line.points;
        var maxValue = Number.NEGATIVE_INFINITY,
            maxPoint = undefined;
        points.forEach(function(point) {
          var y = signal[point.x];
          if (y > maxValue) {
            maxPoint = point;
            maxValue = y;
          }
        });
        return maxPoint.serialize();
      });
      return peaks;
    }
    
    return findPeaks;
  };

  var version = "0.0.1";

  exports.version = version;
  exports.ricker = ricker;
  exports.convolve = convolve;
  exports.findPeaks = findPeaks;

}));
},{}],3:[function(require,module,exports){
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

},{}],4:[function(require,module,exports){
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
},{"d3-peaks":2}],5:[function(require,module,exports){
class AutoGain {
    constructor(audioContext) {
        this.gainIncrement = 0.001;
        this.gainDecrement = 0.001;
        this.minVolume = -10;
        this.maxVolume = -10;
        this.maxGain = 1e5;
        this.minGain = -1e5;

        this.gain = audioContext.createGain();
        this.limiter = audioContext.createDynamicsCompressor()
        this.limiter.threshold.value = -40
        this.limiter.ratio.value = 20
        this.limiter.attack.value = 0.001
        this.limiter.release.value = 0.01

        this.now = ()=>audioContext.currentTime
    }

    updateGain(maxDb) {
        if(maxDb <= -180) return
        const currentGain = 20 * Math.log10(this.gain.gain.value);
        this.gain.gain.cancelAndHoldAtTime(this.now())

        let nextGain = currentGain;
        if(maxDb > this.maxVolume) {
            nextGain += (this.minVolume - maxDb) * this.gainDecrement
        } else {
            nextGain += (this.maxVolume - maxDb) * this.gainIncrement
        }
        if(nextGain != currentGain) {
            if(this.isGainValid(nextGain)) {
                nextGain = Math.pow(10, nextGain*0.05);
                this.gain.gain.linearRampToValueAtTime(nextGain, this.now())
            }
        }
    }

    isGainValid(gain) {
        return gain > this.minGain && gain < this.maxGain;
    }

    getCurrentValue() { return this.gain.gain.value }
    connect(destination) { return this.gain.connect(this.limiter).connect(destination) }
}

module.exports = { AutoGain }
},{}],6:[function(require,module,exports){
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
},{"fft.js":3}],7:[function(require,module,exports){
module.exports = {
    ...require("./anal.js"),
    ...require("./autoGain.js"),
    ...require("./fftConvFilter.js"),
    ...require("./notchFilterBank.js"),
}
},{"./anal.js":4,"./autoGain.js":5,"./fftConvFilter.js":6,"./notchFilterBank.js":8}],8:[function(require,module,exports){
const {calculateQs, chromaticFilter} = require('../scales.js')

class NotchFilterBank {
    constructor(audioContext, scale, q) {
        this.audioContext = audioContext
        if(!scale) {
            const scaleObj = chromaticFilter();
            scale = scaleObj.scale;
            q = scaleObj.qs;
        }
        this.scale = scale;
        this.filters = [];
        this.inMeters = [];
        this.gains = [];
        this.outMeters = [];

        this.q = q || calculateQs(this.scale); //console.log(q);

        this.makeFiltersFromScale(this.scale, this.q);
    }

    makeFiltersFromScale(scale, q) {
        this.scale.forEach((freq,i)=>{
          const filter = this.audioContext.createBiquadFilter();
          filter.type = "notch"; filter.frequency.setValueAtTime(freq, this.audioContext.currentTime);
          filter.gain.value = 1;
          // TODO: need to calculate Q more accurately
          filter.Q.value = typeof q == 'object' ? q[i] : q ;
          this.filters.push(filter);
        })

        this.filters.forEach((f,n) => { if(n<this.filters.length-1)f.connect(this.filters[n+1]) });
    }

    getFrequencyResponse(freqs, response) {
        const tmp = new Float32Array(freqs.length);
        const _ = new Float32Array(freqs.length);
        this.filters.forEach(f=>{
            f.getFrequencyResponse(freqs,tmp,_)
            tmp.forEach((t,n)=>response[n]*=t)
        })
    }

    connectInput(input) {
        input.connect(this.filters[0])
    }
    connect(output) {
        this.filters[this.filters.length-1].connect(output)
        return output
    }

}

module.exports = { NotchFilterBank }
},{"../scales.js":12}],9:[function(require,module,exports){
const SquidbackGraph = require('./graph.js')
const Graph = new SquidbackGraph();
const { 
    AutoGain, 
    FFTConvFilter,
    MagnitudesHistory
} = require("./audio-processors/index.js")

class SquidbackFFTProcess {

    constructor(audioContext, fftSize = 1024) {
        this.audioContext = audioContext;
    }

    async start() {
        this.initFFTNode();
        this.initBuffers();
        this.initConv();
        this.initExtremesFilters();
        this.initGain();
        await this.initInput();
        if(this.input) this.connectAll();
        requestAnimationFrame(()=>this.update());
    }

    initGain() {
        this.autoGain = new AutoGain(this.audioContext)
    }

    initBuffers() {
        const numBins = this.fftNode.frequencyBinCount;
        this.numBins = numBins;
        this.binToFreq = this.audioContext.sampleRate / 2 / numBins;
        this.fftBuffer = new Float32Array(numBins);
        this.outFftBuffer = new Float32Array(numBins);
        this.fbBuffer = new Int32Array(numBins);
        this.numFb = 0;

        this.maxDb = -180;

        this.anal = new MagnitudesHistory(numBins, 10, -60);
    }

    initConv() {
        this.conv = new FFTConvFilter(this.audioContext, this.numBins)
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter();
        hipass.type = "highpass"; hipass.frequency.setValueAtTime(this.binToFreq * 2, this.now());
        hipass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        const lopass = this.audioContext.createBiquadFilter();
        lopass.type = "lowpass"; lopass.frequency.setValueAtTime(15000, this.now());
        lopass.Q.setValueAtTime(0.707, this.audioContext.currentTime);

        this.hipass = hipass; this.lopass = lopass;
    }

    initFFTNode() {
        this.fftNode = this.audioContext.createAnalyser();
        this.fftNode.fftSize = 512;
        this.fftNode.smoothingTimeConstant = 0.9;

        this.outFftNode = this.audioContext.createAnalyser();
        this.outFftNode.fftSize = this.fftNode.fftSize;
        this.outFftNode.smoothingTimeConstant = this.fftNode.smoothingTimeConstant;
    }

    async initInput() {
        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await navigator.mediaDevices.getUserMedia(constraints);
           this.input = this.audioContext.createMediaStreamSource(stream);
           console.log("[Input] Got access to input device")
         } catch(err) {
            console.error("[Input] Can't access user input device")
         }
    }

    connectAll() {
        const signal = this.input
        .connect(this.hipass)
        .connect(this.lopass)

        signal.connect(this.fftNode)

        signal.connect(this.conv.node).connect(this.autoGain.gain).connect(this.autoGain.limiter).connect(this.audioContext.destination)
        this.conv.connect(this.outFftNode);
    }

    update(){
        requestAnimationFrame(()=>this.update())
        this.updateSpectrum();
        this.autoGain.updateGain(this.anal.maxDb);
        this.conv.updateKernel(this.anal.magReductionsAmp);
        this.drawSpectrum();
    }

    updateSpectrum() {
        this.fftNode.getFloatFrequencyData(this.fftBuffer);
        this.outFftNode.getFloatFrequencyData(this.outFftBuffer);
        this.maxDb = this.fftNode.maxDecibels;
        this.anal.analyseSpectrum(this.fftBuffer, this.maxDb);
    }

    // graphs
    drawSpectrum() {
        //const pitchI = this.anal.fbIndexes[0];
        const pitchI = this.fftBuffer.indexOf(this.anal.maxDb);
        const bgColor = this.anal.nbPeaks > 0 ? Graph.pitchAmpToHSLA(pitchI * this.binToFreq, this.fftBuffer[pitchI]) : 'rgba(100,100,100,0.3)';
        //Graph.drawBg("canvas#spectrum", 'rgba(0,0,0,0.5)')
        Graph.drawBg("canvas#spectrum", bgColor)
        Graph.drawGain("canvas#spectrum", "red", this.autoGain.getCurrentValue(), -100, 80)

        //Graph.drawSmoothLogLine("canvas#spectrum", this.fftBuffer, -100, 0);
        Graph.drawLogLine("canvas#spectrum", this.fftBuffer, -100, 0);

        // Graph.drawSmoothLogLine("canvas#spectrum", this.outFftBuffer, -100, 0, "rgba(0,0,0,0.1)");
        
        //Graph.drawSmoothLogLineInverted("canvas#spectrum", this.anal.magReductions, -100, 0, "rgba(0,0,0,0.5)");
        Graph.drawLogLineInverted("canvas#spectrum", this.anal.magReductions, -100, 0, "rgba(0,0,0,0.5)");

        /*
        Graph.drawHLine("canvas#spectrum",  this.anal.averageDb, -100, 0, "green");
        //Graph.drawHLine("canvas#spectrum", this.anal.maxDb, -100, 0, "red");
        Graph.drawHLine("canvas#spectrum", this.anal.peakThr, -100, 0, "blue");
        */

        Graph.drawVerticals("canvas#spectrum", this.anal.nbPeaks, this.anal.peakIndexes, -80, 0, "rgba(255,255,255,0.01)")
        Graph.drawVerticals("canvas#spectrum", this.numFb, this.fbBuffer, -80, 0, "rgba(255,0,0,0.5)")

        //Graph.drawBg("canvas#slope", "white")
        //Graph.drawLine("canvas#slope", this.conv.irBuffer.getChannelData(0), -1, 1, "rgba(255,0,0,1)");
    }
}

module.exports = SquidbackFFTProcess
},{"./audio-processors/index.js":7,"./graph.js":10}],10:[function(require,module,exports){
const Spline = require("cubic-spline")

class SquidbackGraph {

    drawBg(selector, color) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        canvasCtx.beginPath();
        canvasCtx.fillStyle = color;
        canvasCtx.fillRect(0, 0, canvas.width, canvas.height);
        canvasCtx.closePath();
    }

    drawGain(selector, color, amp=0, min=0, max=1) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const gainHeight = (20 * Math.log10(amp) - min) / (max-min);
        //console.log(amp, gainHeight, min, max)
        canvasCtx.beginPath();
        canvasCtx.fillStyle = color;
        canvasCtx.fillRect(canvas.width-2, canvas.height, canvas.width, -canvas.height * gainHeight);
        canvasCtx.closePath();
    }

    pitchAmpToHSLA(pitch,amp=1,alpha=0.1){
        let octaves = Math.log2(pitch/440);
        let pc = 12 * octaves % 12;
        let h = pc/12*360;
        let l = 10 * (octaves + 5.5);
        let s = 255 + amp;
        let code = `hsla(${h},${s<0?0:s>100?100:s}%,${l<5?5:l>100?100:l}%,${alpha})`;
        //console.log(pitch,pc,octaves,s,h);
        return code;
    }

    drawLine(selector, data, min, max, color="rgba(255,255,255,1)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / data.length);
        let barHeight;
        let x = 0;
        canvasCtx.beginPath();
        canvasCtx.moveTo(x, canvas.height - (data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - (min) ) / (max-min) * canvas.height;
            canvasCtx.lineTo(x + barWidth, canvas.height-barHeight);
            x += barWidth;
        }
        canvasCtx.strokeStyle = color;
        canvasCtx.strokeWidth = 0.1;
        canvasCtx.stroke();
    }

    drawLogLine(selector, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length)
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height);
        canvasCtx.lineTo(0, canvas.height - (data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawSmoothLogLine(selector, data, min, max, color="rgba(255,255,255,0.5)", upsample=3, stroke = false) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length);
        const xs = [...data.keys()];
        data = new Spline(xs, data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height);
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = data.at(i); // data[i]
//console.log(i)
            barHeight = ( val - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        if(stroke) {
            canvasCtx.strokeStyle = color
            canvasCtx.stroke();
        } else {
            canvasCtx.fillStyle = color
            canvasCtx.fill();
        }

    }

    drawLogLineInverted(selector, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length - 1)
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, 0);
        canvasCtx.lineTo(0, canvas.height-(data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawSmoothLogLineInverted(selector, data, min, max, color="rgba(255,255,255,0.5)", upsample = 3) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length - 1)
        const xs = [...data.keys()];
        data = new Spline(xs, data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0.5, 0);
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = data.at(i); // data[i]
            barHeight = (val - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawVerticals(selector, numPoints, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        const logIToX = canvas.width / Math.log(data.length - 1);
        canvasCtx.fillStyle = color;
        canvasCtx.beginPath();
        
        for(let n = 0; n < numPoints; ++n){
            const i = data[n];
            //console.log(i);
            const x = Math.log(i) * logIToX;
            canvasCtx.fillRect(x, 0, Math.log(i+1) * logIToX - x, canvas.height)
            //canvasCtx.moveTo(x + barWidth*0.5, canvas.height)
            //canvasCtx.lineTo(x + barWidth*0.5,0);
        }
        //canvasCtx.strokeStyle = color
        //canvasCtx.lineWidth = 1;
        //canvasCtx.stroke();
        canvasCtx.fill();
    }

    drawHLine(selector, h, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const y = canvas.height * (1 - (h - min) / (max-min));
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, y);
        canvasCtx.lineTo(canvas.width, y);
        canvasCtx.strokeStyle = color
        canvasCtx.lineWidth = 0.5;
        canvasCtx.stroke();
    }

}

module.exports = SquidbackGraph
},{"cubic-spline":1}],11:[function(require,module,exports){
"use strict";

const SquidbackWorkletProcess = require('./workletProcess.js')
const SquidbackFFTProcess = require('./fftProcess.js')

async function init() {
    const audioContext = new AudioContext()
    if (audioContext.audioWorklet === undefined) {
        handleNoWorklet();
        return;
    }
    let process = new SquidbackFFTProcess(audioContext, 512);
    await process.start()
}

function handleNoWorklet() {
    let $noWorklet = document.querySelector("#no-worklet");
    $noWorklet.style.display = 'block';
    let $timeline = document.querySelector(".timeline");
    $timeline.style.display = 'none';
    let $controls = document.querySelector(".controls");
    $controls.style.display = 'none';
}

function resize(canvas) {
  // Lookup the size the browser is displaying the canvas.
  var displayWidth  = canvas.clientWidth;
  var displayHeight = canvas.clientHeight;
 
  // Check if the canvas is not the same size.
  if (canvas.width  != displayWidth ||
      canvas.height != displayHeight) {
 
    // Make the canvas the same size
    canvas.width  = displayWidth;
    canvas.height = displayHeight;
  }
}

function resizeAllCanvas() {
    document.querySelectorAll("canvas").forEach(canvas=>resize(canvas))
}


window.addEventListener('load', ()=>{
    const button = document.querySelector("button#start")
    button.addEventListener("click", ()=>{
        init();
        button.classList.add("started")
    });
    resizeAllCanvas();
});

window.addEventListener('resize', ()=>{
    resizeAllCanvas()
});

},{"./fftProcess.js":9,"./workletProcess.js":13}],12:[function(require,module,exports){
const partchScale7 = [30,30.375,31.5,32,33.33333,33.75,34.28571,35,35.55556,36,37.5,38.57143,39.375,40,40.5,42,42.85714,44.44444,45,46.66667,47.14286,48,50,50.625,51.42857,52.5,53.33333,54,56.25,57.14286,59.25926,60,60.75,63,64,66.66667,67.5,68.57143,70,71.11111,72,75,77.14286,78.75,80,81,84,85.71429,88.88889,90,93.33333,94.28571,96,100,101.25,102.8571,105,106.6667,108,112.5,114.2857,118.5185,120,121.5,126,128,133.3333,135,137.1429,140,142.2222,144,150,154.2857,157.5,160,162,168,171.4286,177.7778,180,186.6667,188.5714,192,200,202.5,205.7143,210,213.3333,216,225,228.5714,237.037,240,243,252,256,266.6667,270,274.2857,280,284.4444,288,300,308.5714,315,320,324,336,342.8571,355.5556,360,373.3333,377.1429,384,400,405,411.4286,420,426.6667,432,450,457.1429,474.0741,480,486,504,512,533.3333,540,548.5714,560,568.8889,576,600,617.1429,630,640,648,672,685.7143,711.1111,720,746.6667,754.2857,768,800,810,822.8571,840,853.3333,864,900,914.2857,948.1481,960,972,1008,1024,1066.667,1080,1097.143,1120,1137.778,1152,1200,1234.286,1260,1280,1296,1344,1371.429,1422.222,1440,1493.333,1508.571,1536,1600,1620,1645.714,1680,1706.667,1728,1800,1828.571,1896.296,1920,1944,2016,2048,2133.333,2160,2194.286,2240,2275.556,2304,2400,2468.571,2520,2560,2592,2688,2742.857,2844.444,2880,2986.667,3017.143,3072,3200,3240,3291.429,3360,3413.333,3456,3600,3657.143,3792.593,3840,3888,4032,4096,4266.667,4320,4388.571,4480,4551.111,4608,4800,4937.143,5040,5120,5184,5376,5485.714,5688.889,5760,5973.333,6034.286,6144,6400,6480,6582.857,6720,6826.667,6912,7200,7314.286,7585.185,7680,7776,8064,8192,8533.333,8640,8777.143,8960,9102.222,9216,9600,9874.286,10080,10240,10368,10752,10971.43,11377.78,11520,11946.67,12068.57,12288,12800,12960,13165.71,13440,13653.33,13824,14400,14628.57,15170.37,15360,15552,16128,16384,17066.67,17280,17554.29,17920,18204.44];

function chromatic(min=30,max=20000,div=12){
  var n_octaves = Math.log2(max/min);
  var n_freqs = Math.ceil(n_octaves * div)+1;
  var scale = [min];
  for(let i=1; i<n_freqs; i++){
    scale = [...scale, scale[scale.length-1] * Math.pow(2, 1/div) ];
  };
  return scale;
}

function chromaticFilter(min=30,max=20000,div=12){
    const scale = chromatic(min, max, div);
    const qs = calculateQs(scale);
    const avgs = scale.map((freq,i)=>{
        if(i == 0 || i == scale.length - 1) return freq
        return scale[i-1] + scale[i+1] / 2
    })

    return {scale, qs}
}

function bins(numBins, sr=44100) {
    const binToFreq = sr / 2 / numBins;
    var scale = new Array(numBins);
    var qs = new Array(numBins);
    for(let i = 0; i < numBins; ++i) {
        scale[i] = binToFreq * i;
        qs[i] = i;
    }
    return {scale, qs}

}

function calculateQs(scale, overlap=0){
    return scale.map((freq,freq_i)=>{
        let f0, f1, delta_f;
        if(freq_i==0){
            f0 = freq; f1 = scale[freq_i+1];
        } else if(freq_i == (scale.length-1)) {
            f1 = freq; f0 = scale[freq_i-1];
        } else {
            f0 = scale[freq_i-1];  f1 = scale[freq_i+1];
            /*if(overlap > 0){
              delta_f = Math.max(delta_f0,delta_f1) * 2 * overlap;
            }else{
              delta_f = Math.min(delta_f0,delta_f1) * 2;
            }*/
        }
        delta_f = f1 - f0;
        freq =  (f1+f0) / 2;
        return freq/delta_f
    });
}

module.exports = { partchScale7, chromatic, chromaticFilter, calculateQs, bins };

},{}],13:[function(require,module,exports){
const SquidbackGraph = require('./graph.js')
const { chromaticFilter } = require('./scales.js')
const Graph = new SquidbackGraph();

class SquidbackWorkletProcess {

    constructor(audioContext, fftSize = 1024) {
        this.audioContext = audioContext;
        this.fftSize = fftSize;
        this.gainIncrement = 0.001;
        this.gainDecrement = 0.001;
        this.minVolume = -10;
        this.maxVolume = -10;

        this.initBuffers();
    }

    async start() {
        await this.audioContext.audioWorklet.addModule('phase-vocoder.js');
        this.initExtremesFilters();
        this.initFFTNode();
        this.initGain();
        //this.initFilter();
        await this.initInput();
        if(this.input) this.connectAll();
    }

    initFilter() {
        /*const scale = chromaticFilter(this.binToFreq, 22000, 12)
        this.filter = new NotchFilterBank(this.audioContext, scale.scale, 100);
        this.freqResp = new Float32Array(this.fftSize/2+1);
        this.freqResp.fill(1);
        this.filter.getFrequencyResponse(freqs, this.freqResp);*/

        const freqs = new Float32Array(this.fftSize/2+1);
        freqs.forEach((f,n)=>freqs[n]=n*this.binToFreq);
        // console.log(scale.scale)
    }

    initGain() {
        this.mainGain = this.audioContext.createGain();
        this.limiter = this.audioContext.createDynamicsCompressor()
        this.limiter.threshold.value = -40
        this.limiter.ratio.value = 20
        this.limiter.attack.value = 0.001
        this.limiter.release.value = 0.01
    }

    initBuffers() {
        this.binToFreq = this.audioContext.sampleRate / this.fftSize;
        const numBins = this.fftSize/2 + 1;
        this.fftBuffer = new Float32Array(numBins);
        this.peaksBuffer = new Int32Array(numBins);
        this.numPeaks = 0;
        this.fbBuffer = new Int32Array(numBins);
        this.numFb = 0;
        this.slopeBuffer = new Float32Array(numBins);
        this.filterBuffer = new Int32Array(numBins);

        this.maxDb = -180;
    }

    now() {
        return this.audioContext.currentTime
    }

    initExtremesFilters() {
        const hipass = this.audioContext.createBiquadFilter();
        hipass.type = "highpass"; hipass.frequency.setValueAtTime(this.binToFreq, this.now());
        hipass.Q.setValueAtTime(0.707, this.now());

        const lopass = this.audioContext.createBiquadFilter();
        lopass.type = "lowpass"; lopass.frequency.setValueAtTime(18000, this.now());

        this.hipass = hipass; this.lopass = lopass;
    }

    initFFTNode() {
        this.fftNode = new AudioWorkletNode(this.audioContext, 'phase-vocoder-processor', {processorOptions: {blockSize: this.fftSize}});
        this.initWorkletPort();
    }

    initWorkletPort() {
        this.fftNode.port.onmessage = ({data}) => {
            switch(data[0]) {
                case "spectrum": this.fftBuffer.set(data[1]); break;
                case "slope": this.slopeBuffer.set(data[1]); break;
                case "histo": this.histoBuffer.set(data[1]); break;
                case "peaks":
                    this.numPeaks = data[1]
                    this.peaksBuffer.set(data[2])
                    break;
                case "fb":
                    this.numFb = data[1]
                    this.fbBuffer.set(data[2])
                    this.onFbData();
                    break;
                case 'maxDb': this.maxDb = data[1]; break;
                case 'reductions': this.filterBuffer = data[1]; break;
            }
            //const sum = Math.max(...data);
            // if(data[0] > 0) console.log("Fb candidates:", data[0], data[1])
            // console.log(data)
        };
        requestAnimationFrame(()=>this.showSpectrum());
        requestAnimationFrame(()=>this.showSlope());
    }

    connectAll() {
        this.input
        .connect(this.hipass)
        .connect(this.lopass)
        .connect(this.fftNode)

        //this.filter.connectInput(this.lopass)
        //this.filter
        this.lopass.connect(this.mainGain).connect(this.limiter).connect(this.audioContext.destination)

        //this.lopass.connect(this.mainGain).connect(this.audioContext.destination)
    }

    async initInput() {
        try {
           const constraints = {audio: { noiseSuppression: false, echoCancellation: false, autoGainControl: false }};
           const stream = await navigator.mediaDevices.getUserMedia(constraints);
           this.input = this.audioContext.createMediaStreamSource(stream);
           console.log("GOT INPUT")
         } catch(err) {
           /* handle the error */
            console.error("CANT GET USER INPUT")
         }
    }
    

    onFbData() {
        this.updateGain()
    }

    updateGain() {
        if(this.maxDb <= -180) return
        const gain = this.mainGain.gain;
        gain.cancelAndHoldAtTime(this.now())
        const currentGain = 20 * Math.log10(this.mainGain.gain.value);
        let nextGain = currentGain;
        if(this.maxDb > this.maxVolume) {
            nextGain += (this.minVolume - this.maxDb) * this.gainDecrement
            //console.log("gain--", this.maxDb, nextGain)
        } else {
            nextGain += (this.maxVolume - this.maxDb) * this.gainIncrement
            //console.log("gain++",this.maxDb, nextGain)
        }
        //console.log(this.maxDb,nextGain)
        if(nextGain != currentGain) {
            nextGain = Math.pow(10, nextGain*0.05);
            if(isFinite(nextGain) && nextGain > 0)
                gain.linearRampToValueAtTime(nextGain, this.now())
        }

    }

    // graphs
    showSpectrum(){
        requestAnimationFrame(()=>this.showSpectrum())
        this.fftNode.port.postMessage("getSpectrum");
        this.fftNode.port.postMessage("getPeaks");
        this.fftNode.port.postMessage("getFbCandidates");
        this.fftNode.port.postMessage("getMaxDb");
        this.fftNode.port.postMessage("getMagReductions");

        const bgColor = this.numFb > 0 ? Graph.pitchAmpToHSLA(this.fbBuffer[0] * this.binToFreq, this.fftBuffer[this.fbBuffer[0]]) : 'rgba(100,100,100,1)';


        //Graph.drawBg("canvas#spectrum", 'rgba(0,0,0,0.5)')
        Graph.drawBg("canvas#spectrum", bgColor)
        Graph.drawGain("canvas#spectrum", "red", this.mainGain.gain.value, -100, 80)

        Graph.drawLogLine("canvas#spectrum", this.fftBuffer, -100, 0);
        //Graph.drawLogLine("canvas#spectrum", this.freqResp, 0, 1, "red");

        Graph.drawLogLineInverted("canvas#spectrum", this.filterBuffer, -100, 0, "rgba(0,0,0,0.5)");

        //Graph.drawVerticals("canvas#spectrum", this.numPeaks, this.peaksBuffer, -80, 0, "rgba(255,255,255,1)")
        //Graph.drawVerticals("canvas#spectrum", this.numFb, this.fbBuffer, -80, 0, "rgba(255,0,0,0.5)")
    }

    showSlope(node){
        requestAnimationFrame(()=>this.showSlope())
        this.fftNode.port.postMessage("getSlope");
        //console.log(slopeBuffer)
        Graph.drawBg("canvas#slope", 'white')
        Graph.drawLogLine("canvas#slope", this.slopeBuffer, 0, 500, "black");
    }
}

module.exports = SquidbackWorkletProcess
},{"./graph.js":10,"./scales.js":12}]},{},[11])
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJub2RlX21vZHVsZXMvY3ViaWMtc3BsaW5lL2luZGV4LmpzIiwibm9kZV9tb2R1bGVzL2QzLXBlYWtzL2J1aWxkL2QzLXBlYWtzLmpzIiwibm9kZV9tb2R1bGVzL2ZmdC5qcy9saWIvZmZ0LmpzIiwic3JjL2F1ZGlvLXByb2Nlc3NvcnMvYW5hbC5qcyIsInNyYy9hdWRpby1wcm9jZXNzb3JzL2F1dG9HYWluLmpzIiwic3JjL2F1ZGlvLXByb2Nlc3NvcnMvZmZ0Q29udkZpbHRlci5qcyIsInNyYy9hdWRpby1wcm9jZXNzb3JzL2luZGV4LmpzIiwic3JjL2F1ZGlvLXByb2Nlc3NvcnMvbm90Y2hGaWx0ZXJCYW5rLmpzIiwic3JjL2ZmdFByb2Nlc3MuanMiLCJzcmMvZ3JhcGguanMiLCJzcmMvbWFpbi5qcyIsInNyYy9zY2FsZXMuanMiLCJzcmMvd29ya2xldFByb2Nlc3MuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDL0lBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzWUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM2ZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqTkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDTEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDeklBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6TEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3hEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24oKXtmdW5jdGlvbiByKGUsbix0KXtmdW5jdGlvbiBvKGksZil7aWYoIW5baV0pe2lmKCFlW2ldKXt2YXIgYz1cImZ1bmN0aW9uXCI9PXR5cGVvZiByZXF1aXJlJiZyZXF1aXJlO2lmKCFmJiZjKXJldHVybiBjKGksITApO2lmKHUpcmV0dXJuIHUoaSwhMCk7dmFyIGE9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitpK1wiJ1wiKTt0aHJvdyBhLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsYX12YXIgcD1uW2ldPXtleHBvcnRzOnt9fTtlW2ldWzBdLmNhbGwocC5leHBvcnRzLGZ1bmN0aW9uKHIpe3ZhciBuPWVbaV1bMV1bcl07cmV0dXJuIG8obnx8cil9LHAscC5leHBvcnRzLHIsZSxuLHQpfXJldHVybiBuW2ldLmV4cG9ydHN9Zm9yKHZhciB1PVwiZnVuY3Rpb25cIj09dHlwZW9mIHJlcXVpcmUmJnJlcXVpcmUsaT0wO2k8dC5sZW5ndGg7aSsrKW8odFtpXSk7cmV0dXJuIG99cmV0dXJuIHJ9KSgpIiwibW9kdWxlLmV4cG9ydHMgPSBjbGFzcyBTcGxpbmUge1xuICBjb25zdHJ1Y3Rvcih4cywgeXMpIHtcbiAgICB0aGlzLnhzID0geHM7XG4gICAgdGhpcy55cyA9IHlzO1xuICAgIHRoaXMua3MgPSB0aGlzLmdldE5hdHVyYWxLcyhuZXcgRmxvYXQ2NEFycmF5KHRoaXMueHMubGVuZ3RoKSk7XG4gIH1cblxuICBnZXROYXR1cmFsS3Moa3MpIHtcbiAgICBjb25zdCBuID0gdGhpcy54cy5sZW5ndGggLSAxO1xuICAgIGNvbnN0IEEgPSB6ZXJvc01hdChuICsgMSwgbiArIDIpO1xuXG4gICAgZm9yIChcbiAgICAgIGxldCBpID0gMTtcbiAgICAgIGkgPCBuO1xuICAgICAgaSsrIC8vIHJvd3NcbiAgICApIHtcbiAgICAgIEFbaV1baSAtIDFdID0gMSAvICh0aGlzLnhzW2ldIC0gdGhpcy54c1tpIC0gMV0pO1xuICAgICAgQVtpXVtpXSA9XG4gICAgICAgIDIgKlxuICAgICAgICAoMSAvICh0aGlzLnhzW2ldIC0gdGhpcy54c1tpIC0gMV0pICsgMSAvICh0aGlzLnhzW2kgKyAxXSAtIHRoaXMueHNbaV0pKTtcbiAgICAgIEFbaV1baSArIDFdID0gMSAvICh0aGlzLnhzW2kgKyAxXSAtIHRoaXMueHNbaV0pO1xuICAgICAgQVtpXVtuICsgMV0gPVxuICAgICAgICAzICpcbiAgICAgICAgKCh0aGlzLnlzW2ldIC0gdGhpcy55c1tpIC0gMV0pIC9cbiAgICAgICAgICAoKHRoaXMueHNbaV0gLSB0aGlzLnhzW2kgLSAxXSkgKiAodGhpcy54c1tpXSAtIHRoaXMueHNbaSAtIDFdKSkgK1xuICAgICAgICAgICh0aGlzLnlzW2kgKyAxXSAtIHRoaXMueXNbaV0pIC9cbiAgICAgICAgICAgICgodGhpcy54c1tpICsgMV0gLSB0aGlzLnhzW2ldKSAqICh0aGlzLnhzW2kgKyAxXSAtIHRoaXMueHNbaV0pKSk7XG4gICAgfVxuXG4gICAgQVswXVswXSA9IDIgLyAodGhpcy54c1sxXSAtIHRoaXMueHNbMF0pO1xuICAgIEFbMF1bMV0gPSAxIC8gKHRoaXMueHNbMV0gLSB0aGlzLnhzWzBdKTtcbiAgICBBWzBdW24gKyAxXSA9XG4gICAgICAoMyAqICh0aGlzLnlzWzFdIC0gdGhpcy55c1swXSkpIC9cbiAgICAgICgodGhpcy54c1sxXSAtIHRoaXMueHNbMF0pICogKHRoaXMueHNbMV0gLSB0aGlzLnhzWzBdKSk7XG5cbiAgICBBW25dW24gLSAxXSA9IDEgLyAodGhpcy54c1tuXSAtIHRoaXMueHNbbiAtIDFdKTtcbiAgICBBW25dW25dID0gMiAvICh0aGlzLnhzW25dIC0gdGhpcy54c1tuIC0gMV0pO1xuICAgIEFbbl1bbiArIDFdID1cbiAgICAgICgzICogKHRoaXMueXNbbl0gLSB0aGlzLnlzW24gLSAxXSkpIC9cbiAgICAgICgodGhpcy54c1tuXSAtIHRoaXMueHNbbiAtIDFdKSAqICh0aGlzLnhzW25dIC0gdGhpcy54c1tuIC0gMV0pKTtcblxuICAgIHJldHVybiBzb2x2ZShBLCBrcyk7XG4gIH1cblxuICAvKipcbiAgICogaW5zcGlyZWQgYnkgaHR0cHM6Ly9zdGFja292ZXJmbG93LmNvbS9hLzQwODUwMzEzLzQ0MTczMjdcbiAgICovXG4gIGdldEluZGV4QmVmb3JlKHRhcmdldCkge1xuICAgIGxldCBsb3cgPSAwO1xuICAgIGxldCBoaWdoID0gdGhpcy54cy5sZW5ndGg7XG4gICAgbGV0IG1pZCA9IDA7XG4gICAgd2hpbGUgKGxvdyA8IGhpZ2gpIHtcbiAgICAgIG1pZCA9IE1hdGguZmxvb3IoKGxvdyArIGhpZ2gpIC8gMik7XG4gICAgICBpZiAodGhpcy54c1ttaWRdIDwgdGFyZ2V0ICYmIG1pZCAhPT0gbG93KSB7XG4gICAgICAgIGxvdyA9IG1pZDtcbiAgICAgIH0gZWxzZSBpZiAodGhpcy54c1ttaWRdID49IHRhcmdldCAmJiBtaWQgIT09IGhpZ2gpIHtcbiAgICAgICAgaGlnaCA9IG1pZDtcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIGhpZ2ggPSBsb3c7XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVybiBsb3cgKyAxO1xuICB9XG5cbiAgYXQoeCkge1xuICAgIGxldCBpID0gdGhpcy5nZXRJbmRleEJlZm9yZSh4KTtcbiAgICBjb25zdCB0ID0gKHggLSB0aGlzLnhzW2kgLSAxXSkgLyAodGhpcy54c1tpXSAtIHRoaXMueHNbaSAtIDFdKTtcbiAgICBjb25zdCBhID1cbiAgICAgIHRoaXMua3NbaSAtIDFdICogKHRoaXMueHNbaV0gLSB0aGlzLnhzW2kgLSAxXSkgLVxuICAgICAgKHRoaXMueXNbaV0gLSB0aGlzLnlzW2kgLSAxXSk7XG4gICAgY29uc3QgYiA9XG4gICAgICAtdGhpcy5rc1tpXSAqICh0aGlzLnhzW2ldIC0gdGhpcy54c1tpIC0gMV0pICtcbiAgICAgICh0aGlzLnlzW2ldIC0gdGhpcy55c1tpIC0gMV0pO1xuICAgIGNvbnN0IHEgPVxuICAgICAgKDEgLSB0KSAqIHRoaXMueXNbaSAtIDFdICtcbiAgICAgIHQgKiB0aGlzLnlzW2ldICtcbiAgICAgIHQgKiAoMSAtIHQpICogKGEgKiAoMSAtIHQpICsgYiAqIHQpO1xuICAgIHJldHVybiBxO1xuICB9XG59O1xuXG5mdW5jdGlvbiBzb2x2ZShBLCBrcykge1xuICBjb25zdCBtID0gQS5sZW5ndGg7XG4gIGxldCBoID0gMDtcbiAgbGV0IGsgPSAwO1xuICB3aGlsZSAoaCA8IG0gJiYgayA8PSBtKSB7XG4gICAgbGV0IGlfbWF4ID0gMDtcbiAgICBsZXQgbWF4ID0gLUluZmluaXR5O1xuICAgIGZvciAobGV0IGkgPSBoOyBpIDwgbTsgaSsrKSB7XG4gICAgICBjb25zdCB2ID0gTWF0aC5hYnMoQVtpXVtrXSk7XG4gICAgICBpZiAodiA+IG1heCkge1xuICAgICAgICBpX21heCA9IGk7XG4gICAgICAgIG1heCA9IHY7XG4gICAgICB9XG4gICAgfVxuXG4gICAgaWYgKEFbaV9tYXhdW2tdID09PSAwKSB7XG4gICAgICBrKys7XG4gICAgfSBlbHNlIHtcbiAgICAgIHN3YXBSb3dzKEEsIGgsIGlfbWF4KTtcbiAgICAgIGZvciAobGV0IGkgPSBoICsgMTsgaSA8IG07IGkrKykge1xuICAgICAgICBjb25zdCBmID0gQVtpXVtrXSAvIEFbaF1ba107XG4gICAgICAgIEFbaV1ba10gPSAwO1xuICAgICAgICBmb3IgKGxldCBqID0gayArIDE7IGogPD0gbTsgaisrKSBBW2ldW2pdIC09IEFbaF1bal0gKiBmO1xuICAgICAgfVxuICAgICAgaCsrO1xuICAgICAgaysrO1xuICAgIH1cbiAgfVxuXG4gIGZvciAoXG4gICAgbGV0IGkgPSBtIC0gMTtcbiAgICBpID49IDA7XG4gICAgaS0tIC8vIHJvd3MgPSBjb2x1bW5zXG4gICkge1xuICAgIHZhciB2ID0gMDtcbiAgICBpZiAoQVtpXVtpXSkge1xuICAgICAgdiA9IEFbaV1bbV0gLyBBW2ldW2ldO1xuICAgIH1cbiAgICBrc1tpXSA9IHY7XG4gICAgZm9yIChcbiAgICAgIGxldCBqID0gaSAtIDE7XG4gICAgICBqID49IDA7XG4gICAgICBqLS0gLy8gcm93c1xuICAgICkge1xuICAgICAgQVtqXVttXSAtPSBBW2pdW2ldICogdjtcbiAgICAgIEFbal1baV0gPSAwO1xuICAgIH1cbiAgfVxuICByZXR1cm4ga3M7XG59XG5cbmZ1bmN0aW9uIHplcm9zTWF0KHIsIGMpIHtcbiAgY29uc3QgQSA9IFtdO1xuICBmb3IgKGxldCBpID0gMDsgaSA8IHI7IGkrKykgQS5wdXNoKG5ldyBGbG9hdDY0QXJyYXkoYykpO1xuICByZXR1cm4gQTtcbn1cblxuZnVuY3Rpb24gc3dhcFJvd3MobSwgaywgbCkge1xuICBsZXQgcCA9IG1ba107XG4gIG1ba10gPSBtW2xdO1xuICBtW2xdID0gcDtcbn1cbiIsIihmdW5jdGlvbiAoZ2xvYmFsLCBmYWN0b3J5KSB7XG4gIHR5cGVvZiBleHBvcnRzID09PSAnb2JqZWN0JyAmJiB0eXBlb2YgbW9kdWxlICE9PSAndW5kZWZpbmVkJyA/IGZhY3RvcnkoZXhwb3J0cykgOlxuICB0eXBlb2YgZGVmaW5lID09PSAnZnVuY3Rpb24nICYmIGRlZmluZS5hbWQgPyBkZWZpbmUoWydleHBvcnRzJ10sIGZhY3RvcnkpIDpcbiAgKGZhY3RvcnkoKGdsb2JhbC5kM19wZWFrcyA9IGdsb2JhbC5kM19wZWFrcyB8fCB7fSkpKTtcbn0odGhpcywgZnVuY3Rpb24gKGV4cG9ydHMpIHsgJ3VzZSBzdHJpY3QnO1xuXG4gIC8qKlxuICAgKiBTZWUgaHR0cHM6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvTWV4aWNhbl9oYXRfd2F2ZWxldFxuICAgKi9cbiAgZnVuY3Rpb24gcmlja2VyKCkge1xuICAgIHZhciDPgyA9IDE7XG4gICAgXG4gICAgdmFyIHJpY2tlciA9IGZ1bmN0aW9uKHQpIHtcbiAgICAgIHZhciB0MiA9IHQqdCxcbiAgICAgICAgICB2YXJpYW5jZSA9IM+DKs+DO1xuICAgICAgXG4gICAgICB2YXIgQyA9IDIuMCAvICggTWF0aC5zcXJ0KDMgKiDPgykgKiAoTWF0aC5wb3coTWF0aC5QSSwgMC4yNSkpICk7XG4gICAgICB2YXIgbm9ybSA9ICgxLjAgLSAodDIpLyh2YXJpYW5jZSkpO1xuICAgICAgdmFyIGdhdXNzID0gTWF0aC5leHAoIC0odDIpIC8gKDIqdmFyaWFuY2UpICk7XG4gICAgICBcbiAgICAgIHJldHVybiBDKm5vcm0qZ2F1c3M7XG4gICAgfVxuICAgIFxuICAgIHJpY2tlci5zdGQgPSBmdW5jdGlvbihfKSB7XG4gICAgICByZXR1cm4gYXJndW1lbnRzLmxlbmd0aCA/ICjPgyA9IF8sIHJpY2tlcikgOiDPgztcbiAgICB9XG4gICAgXG4gICAgLyoqXG4gICAgICogUmFuZ2Ugb2YgcG9pbnRzIHRvIHNhbXBsZSBmcm9tIHRoZSB3YXZlbGV0LiBbLXJlYWNoLCByZWFjaF1cbiAgICAgKi9cbiAgICByaWNrZXIucmVhY2ggPSBmdW5jdGlvbigpIHtcbiAgICAgIHJldHVybiA1ICogz4M7XG4gICAgfVxuICAgIFxuICAgIHJldHVybiByaWNrZXI7XG4gIH07XG5cbiAgZnVuY3Rpb24gY29udm9sdmUoKSB7XG4gICAgdmFyIGtlcm5lbCA9IHJpY2tlcigpO1xuICAgIFxuICAgIC8qKlxuICAgICAqIHlbbl0gPSBTdW1fa3t4W2tdICogaFtuLWtdfVxuICAgICAqIHk6IG91dHB1dFxuICAgICAqIHg6IGlucHV0XG4gICAgICogaDogc21vb3RoZXJcbiAgICAgKi9cbiAgICB2YXIgY29udm9sdmUgPSBmdW5jdGlvbihzaWduYWwpIHtcbiAgICAgIHZhciBzaXplID0gc2lnbmFsLmxlbmd0aCxcbiAgICAgICAgICBuID0gLTEsXG4gICAgICAgICAgY29udm9sdXRpb24gPSBuZXcgQXJyYXkoc2l6ZSk7XG4gICAgICAgICAgXG4gICAgICB3aGlsZSAoKytuIDwgc2l6ZSkge1xuICAgICAgICB2YXIgeSA9IDA7XG4gICAgICAgIFxuICAgICAgICB2YXIgYm94ID0gYm91bmRpbmdCb3gobiwga2VybmVsLnJlYWNoKCksIDAsIHNpemUgLSAxKTtcbiAgICAgICAgYm94LmZvckVhY2goZnVuY3Rpb24ozrQpIHtcbiAgICAgICAgICB2YXIgayA9IG4gKyDOtDtcbiAgICAgICAgICB5ICs9IHNpZ25hbFtrXSAqIGtlcm5lbCjOtCk7XG4gICAgICAgIH0pO1xuICAgICAgICBjb252b2x1dGlvbltuXSA9IHk7XG4gICAgICB9XG4gICAgICBcbiAgICAgIHJldHVybiBjb252b2x1dGlvbjtcbiAgICB9O1xuICAgIFxuICAgIGNvbnZvbHZlLmtlcm5lbCA9IGZ1bmN0aW9uKF8pIHtcbiAgICAgIHJldHVybiBhcmd1bWVudHMubGVuZ3RoID8gKGtlcm5lbCA9IF8sIGNvbnZvbHZlKSA6IGtlcm5lbDtcbiAgICB9XG4gICAgXG4gICAgZnVuY3Rpb24gcmFuZ2UocmVhY2gpIHtcbiAgICAgIHJlYWNoID0gK3JlYWNoO1xuICAgICAgdmFyIGkgPSAtMSxcbiAgICAgICAgICBuID0gMipyZWFjaCArIDEsXG4gICAgICAgICAgcmFuZ2UgPSBuZXcgQXJyYXkobik7XG4gICAgICB3aGlsZSgrK2kgPCBuKSB7XG4gICAgICAgIHJhbmdlW2ldID0gKC1yZWFjaCkgKyBpO1xuICAgICAgfVxuICAgICAgcmV0dXJuIHJhbmdlO1xuICAgIH1cbiAgICBcbiAgICBmdW5jdGlvbiBib3VuZGluZ0JveChuLCByZWFjaCwgbG8sIGhpKSB7XG4gICAgICBmb3IgKHZhciBpID0gMTsgaSA8PSByZWFjaDsgaSsrKSB7XG4gICAgICAgIHZhciBsZWZ0ICA9IG4gLSBpLFxuICAgICAgICAgICAgcmlnaHQgPSBuICsgaTtcbiAgICAgICAgaWYgKGxlZnQgPj0gbG8gJiYgcmlnaHQgPD0gaGkpIGNvbnRpbnVlO1xuICAgICAgICByZXR1cm4gcmFuZ2UoaSAtIDEpO1xuICAgICAgfVxuICAgICAgcmV0dXJuIHJhbmdlKHJlYWNoKTtcbiAgICB9XG4gICAgXG4gICAgcmV0dXJuIGNvbnZvbHZlO1xuICB9O1xuXG4gIGZ1bmN0aW9uIGlzTG9jYWxNYXhpbWEoYXJyLCBpbmRleCkge1xuICAgIHZhciBjdXJyZW50ID0gYXJyW2luZGV4XSxcbiAgICAgICAgbGVmdCA9IGFycltpbmRleCAtIDFdLFxuICAgICAgICByaWdodCA9IGFycltpbmRleCArIDFdO1xuICAgICAgICBcbiAgICBpZiAobGVmdCAhPT0gdW5kZWZpbmVkICYmIHJpZ2h0ICE9PSB1bmRlZmluZWQpIHtcbiAgICAgIGlmIChjdXJyZW50ID4gbGVmdCAmJiBjdXJyZW50ID4gcmlnaHQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICAgIGVsc2UgaWYgKGN1cnJlbnQgPj0gbGVmdCAmJiBjdXJyZW50ID4gcmlnaHQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICAgIGVsc2UgaWYgKGN1cnJlbnQgPiBsZWZ0ICYmIGN1cnJlbnQgPj0gcmlnaHQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICB9XG4gICAgZWxzZSBpZiAobGVmdCAhPT0gdW5kZWZpbmVkICYmIGN1cnJlbnQgPiBsZWZ0KSB7IHJldHVybiB0cnVlOyB9XG4gICAgZWxzZSBpZiAocmlnaHQgIT09IHVuZGVmaW5lZCAmJiBjdXJyZW50ID4gcmlnaHQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICBcbiAgICByZXR1cm4gZmFsc2U7XG4gIH1cblxuICAvKipcbiAgICogQHBhcmFtIHthcnJ9IHJvdyBpbiB0aGUgQ1dUIG1hdHJpeC5cbiAgICogQHJldHVybiBBcnJheSBvZiBpbmRpY2VzIHdpdGggcmVsYXRpdmUgbWF4aW1hcy5cbiAgICovXG4gIGZ1bmN0aW9uIG1heGltYXMoYXJyKSB7XG4gICAgdmFyIG1heGltYXMgPSBbXTtcbiAgICBhcnIuZm9yRWFjaChmdW5jdGlvbih2YWx1ZSwgaW5kZXgpIHtcbiAgICAgIGlmIChpc0xvY2FsTWF4aW1hKGFyciwgaW5kZXgpKSBtYXhpbWFzLnB1c2goe3g6IGluZGV4LCB5OiB2YWx1ZX0pO1xuICAgIH0pO1xuICAgIHJldHVybiBtYXhpbWFzO1xuICB9O1xuXG4gIGZ1bmN0aW9uIG5lYXJlc3ROZWlnaGJvcihsaW5lLCBtYXhpbWFzLCB3aW5kb3cpIHtcbiAgICB2YXIgY2FjaGUgPSB7fTtcbiAgICBtYXhpbWFzLmZvckVhY2goZnVuY3Rpb24oZCkge1xuICAgICAgY2FjaGVbZC54XSA9IGQueTtcbiAgICB9KTtcbiAgICBcbiAgICB2YXIgcG9pbnQgPSBsaW5lLnRvcCgpO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDw9IHdpbmRvdzsgaSsrKSB7XG4gICAgICB2YXIgbGVmdCA9IHBvaW50LnggKyBpO1xuICAgICAgdmFyIHJpZ2h0ID0gcG9pbnQueCAtIGk7XG4gICAgICBcbiAgICAgIGlmICggKGxlZnQgaW4gY2FjaGUpICYmIChyaWdodCBpbiBjYWNoZSkgKSB7XG4gICAgICAgIGlmIChjYWNoZVtsZWZ0XSA+IGNhY2hlW3JpZ2h0XSkge1xuICAgICAgICAgIHJldHVybiBsZWZ0O1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiByaWdodDtcbiAgICAgIH1cbiAgICAgIGVsc2UgaWYgKGxlZnQgaW4gY2FjaGUpIHtcbiAgICAgICAgcmV0dXJuIGxlZnQ7XG4gICAgICB9XG4gICAgICBlbHNlIGlmIChyaWdodCBpbiBjYWNoZSkge1xuICAgICAgICByZXR1cm4gcmlnaHQ7XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVybiBudWxsO1xuICB9XG5cbiAgZnVuY3Rpb24gcGVyY2VudGlsZShhcnIsIHBlcmMpIHtcbiAgICB2YXIgbGVuZ3RoID0gYXJyLmxlbmd0aDtcbiAgICB2YXIgaW5kZXggPSBNYXRoLm1pbihsZW5ndGggLSAxLCBNYXRoLmNlaWwocGVyYyAqIGxlbmd0aCkpO1xuICAgIFxuICAgIGFyci5zb3J0KGZ1bmN0aW9uKGEsIGIpIHsgcmV0dXJuIGEgLSBiOyB9KTtcbiAgICByZXR1cm4gYXJyW2luZGV4XTtcbiAgfVxuXG4gIGZ1bmN0aW9uIFBvaW50KHgsIHksIHdpZHRoKSB7XG4gICAgdGhpcy54ID0geDtcbiAgICB0aGlzLnkgPSB5O1xuICAgIHRoaXMud2lkdGggPSB3aWR0aDtcbiAgICB0aGlzLnNuciA9IHVuZGVmaW5lZDtcbiAgfVxuXG4gIFBvaW50LnByb3RvdHlwZS5TTlIgPSBmdW5jdGlvbihjb252KSB7XG4gICAgdmFyIHNtb290aGluZ0ZhY3RvciA9IDAuMDAwMDE7XG4gICAgdmFyIHNpZ25hbCA9IHRoaXMueTtcbiAgICBcbiAgICB2YXIgbG93ZXJCb3VuZCA9IE1hdGgubWF4KDAsIHRoaXMueCAtIHRoaXMud2lkdGgpO1xuICAgIHZhciB1cHBlckJvdW5kID0gTWF0aC5taW4oY29udi5sZW5ndGgsIHRoaXMueCArIHRoaXMud2lkdGggKyAxKTtcbiAgICB2YXIgbmVpZ2hib3JzID0gY29udi5zbGljZShsb3dlckJvdW5kLCB1cHBlckJvdW5kKTtcbiAgICB2YXIgbm9pc2UgPSBwZXJjZW50aWxlKG5laWdoYm9ycywgMC45NSk7XG4gICAgXG4gICAgc2lnbmFsICs9IHNtb290aGluZ0ZhY3RvcjtcbiAgICBub2lzZSArPSBzbW9vdGhpbmdGYWN0b3I7XG4gICAgdGhpcy5zbnIgPSBzaWduYWwvbm9pc2U7XG4gICAgcmV0dXJuIHRoaXMuc25yO1xuICB9XG5cbiAgUG9pbnQucHJvdG90eXBlLnNlcmlhbGl6ZSA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB7aW5kZXg6IHRoaXMueCwgd2lkdGg6IHRoaXMud2lkdGgsIHNucjogdGhpcy5zbnJ9O1xuICB9XG5cbiAgZnVuY3Rpb24gUmlkZ2VMaW5lKCkge1xuICAgIHRoaXMucG9pbnRzID0gW107XG4gICAgdGhpcy5nYXAgPSAwO1xuICB9XG5cbiAgLyoqXG4gICAqIElmIHRoZSBwb2ludCBpcyB2YWxpZCBhcHBlbmQgaXQgdG8gdGhlIHJpZGdlbGluZSwgYW5kIHJlc2V0IHRoZSBnYXAuXG4gICAqIE90aGVyd2lzZSwgaW5jcmVtZW50IHRoZSBnYXAgYW5kIGRvIG5vdGhpbmcuXG4gICAqIFxuICAgKiBAcGFyYW0ge3BvaW50fSBQb2ludCBvYmplY3QuXG4gICAqL1xuICBSaWRnZUxpbmUucHJvdG90eXBlLmFkZCA9IGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgaWYgKHBvaW50ID09PSBudWxsIHx8IHBvaW50ID09PSB1bmRlZmluZWQpIHtcbiAgICAgIHRoaXMuZ2FwICs9IDE7XG4gICAgICByZXR1cm47XG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMucG9pbnRzLnB1c2gocG9pbnQpO1xuICAgICAgdGhpcy5nYXAgPSAwO1xuICAgIH1cbiAgfVxuXG4gIC8qKlxuICAgKiBAcmV0dXJuIHtQb2ludH0gTGFzdCBwb2ludCBhZGRlZCBpbnRvIHRoZSByaWRnZWxpbmUuXG4gICAqL1xuICBSaWRnZUxpbmUucHJvdG90eXBlLnRvcCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnBvaW50c1t0aGlzLnBvaW50cy5sZW5ndGggLSAxXTtcbiAgfVxuXG4gIC8qKlxuICAgKiBAcmV0dXJuIHtudW1iZXJ9IExlbmd0aCBvZiBwb2ludHMgb24gdGhlIHJpZGdlbGluZS5cbiAgICovXG4gIFJpZGdlTGluZS5wcm90b3R5cGUubGVuZ3RoID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMucG9pbnRzLmxlbmd0aDtcbiAgfVxuXG4gIC8qKlxuICAgKiBAcmV0dXJuIHtib29sZWFufSBUcnVlIGlmIHRoZSBnYXAgaW4gdGhlIGxpbmUgaXMgYWJvdmUgYSB0aHJlc2hvbGQuIEZhbHNlIG90aGVyd2lzZS5cbiAgICovXG4gIFJpZGdlTGluZS5wcm90b3R5cGUuaXNEaXNjb25uZWN0ZWQgPSBmdW5jdGlvbiAodGhyZXNob2xkKSB7XG4gICAgcmV0dXJuIHRoaXMuZ2FwID4gdGhyZXNob2xkO1xuICB9XG5cbiAgLyoqXG4gICAqIEBwYXJhbSB7QXJyYXl9IFNtYWxsZXN0IHNjYWxlIGluIHRoZSBjb252b2x1dGlvbiBtYXRyaXhcbiAgICovXG4gIFJpZGdlTGluZS5wcm90b3R5cGUuU05SID0gZnVuY3Rpb24oY29udikge1xuICAgIHZhciBtYXhTbnIgPSBOdW1iZXIuTkVHQVRJVkVfSU5GSU5JVFk7XG4gICAgdGhpcy5wb2ludHMuZm9yRWFjaChmdW5jdGlvbihwb2ludCkge1xuICAgICAgdmFyIHNuciA9IHBvaW50LlNOUihjb252KTtcbiAgICAgIGlmIChzbnIgPiBtYXhTbnIpIG1heFNuciA9IHNucjtcbiAgICB9KTtcbiAgICByZXR1cm4gbWF4U25yO1xuICB9XG5cbiAgZnVuY3Rpb24gZmluZFBlYWtzKCkge1xuICAgIHZhciBrZXJuZWwgPSByaWNrZXIsXG4gICAgICAgIGdhcFRocmVzaG9sZCA9IDEsXG4gICAgICAgIG1pbkxpbmVMZW5ndGggPSAxLFxuICAgICAgICBtaW5TTlIgPSAxLjAsXG4gICAgICAgIHdpZHRocyA9IFsxXTtcbiAgICBcbiAgICB2YXIgZmluZFBlYWtzID0gZnVuY3Rpb24oc2lnbmFsKSB7XG4gICAgICB2YXIgTSA9IENXVChzaWduYWwpO1xuICAgICAgXG4gICAgICB2YXIgcmlkZ2VMaW5lcyA9IGluaXRpYWxpemVSaWRnZUxpbmVzKE0pO1xuICAgICAgcmlkZ2VMaW5lcyA9IGNvbm5lY3RSaWRnZUxpbmVzKE0sIHJpZGdlTGluZXMpO1xuICAgICAgcmlkZ2VMaW5lcyA9IGZpbHRlclJpZGdlTGluZXMoc2lnbmFsLCByaWRnZUxpbmVzKTtcbiAgICAgIFxuICAgICAgcmV0dXJuIHBlYWtzKHNpZ25hbCwgcmlkZ2VMaW5lcyk7XG4gICAgfTtcbiAgICBcbiAgICAvKipcbiAgICAgKiBTbW9vdGhpbmcgZnVuY3Rpb24uXG4gICAgICovXG4gICAgZmluZFBlYWtzLmtlcm5lbCA9IGZ1bmN0aW9uKF8pIHtcbiAgICAgIHJldHVybiBhcmd1bWVudHMubGVuZ3RoID8gKGtlcm5lbCA9IF8sIGZpbmRQZWFrcykgOiBrZXJuZWw7XG4gICAgfVxuICAgIFxuICAgIC8qKlxuICAgICAqIEV4cGVjdGVkIHdpZHRocyBvZiB0aGUgcGVha3MuXG4gICAgICovXG4gICAgZmluZFBlYWtzLndpZHRocyA9IGZ1bmN0aW9uKF8pIHtcbiAgICAgIF8uc29ydChmdW5jdGlvbihhLCBiKSB7IHJldHVybiBhIC0gYjsgfSk7XG4gICAgICByZXR1cm4gYXJndW1lbnRzLmxlbmd0aCA/ICh3aWR0aHMgPSBfLCBmaW5kUGVha3MpIDogd2lkdGhzO1xuICAgIH1cbiAgICBcbiAgICAvKipcbiAgICAgKiBOdW1iZXIgb2YgZ2FwcyB0aGF0IHdlIGFsbG93IGluIHRoZSByaWRnZSBsaW5lcy5cbiAgICAgKi9cbiAgICBmaW5kUGVha3MuZ2FwVGhyZXNob2xkID0gZnVuY3Rpb24oXykge1xuICAgICAgcmV0dXJuIGFyZ3VtZW50cy5sZW5ndGggPyAoZ2FwVGhyZXNob2xkID0gXywgZmluZFBlYWtzKSA6IGdhcFRocmVzaG9sZDtcbiAgICB9XG4gICAgXG4gICAgLyoqXG4gICAgICogTWluaW11bSByaWRnZSBsaW5lIGxlbmd0aC5cbiAgICAgKi9cbiAgICBmaW5kUGVha3MubWluTGluZUxlbmd0aCA9IGZ1bmN0aW9uKF8pIHtcbiAgICAgIHJldHVybiBhcmd1bWVudHMubGVuZ3RoID8gKG1pbkxpbmVMZW5ndGggPSBfLCBmaW5kUGVha3MpIDogbWluTGluZUxlbmd0aDtcbiAgICB9XG4gICAgXG4gICAgLyoqXG4gICAgICogTWluaW11bSBzaWduYWwgdG8gbm9pc2UgcmF0aW8gZm9yIHRoZSBwZWFrcy5cbiAgICAgKi9cbiAgICBmaW5kUGVha3MubWluU05SID0gZnVuY3Rpb24oXykge1xuICAgICAgcmV0dXJuIGFyZ3VtZW50cy5sZW5ndGggPyAobWluU05SID0gXywgZmluZFBlYWtzKSA6IG1pblNOUjtcbiAgICB9XG4gICAgXG4gICAgdmFyIENXVCA9IGZ1bmN0aW9uKHNpZ25hbCkge1xuICAgICAgdmFyIE0gPSBuZXcgQXJyYXkod2lkdGhzLmxlbmd0aCk7XG4gICAgICB3aWR0aHMuZm9yRWFjaChmdW5jdGlvbih3aWR0aCwgaSkge1xuICAgICAgICB2YXIgc21vb3RoZXIgPSBrZXJuZWwoKVxuICAgICAgICAgIC5zdGQod2lkdGgpO1xuICAgICAgICB2YXIgdHJhbnNmb3JtID0gY29udm9sdmUoKVxuICAgICAgICAgIC5rZXJuZWwoc21vb3RoZXIpO1xuICAgICAgICBcbiAgICAgICAgdmFyIGNvbnZvbHV0aW9uID0gdHJhbnNmb3JtKHNpZ25hbCk7XG4gICAgICAgIE1baV0gPSBjb252b2x1dGlvbjtcbiAgICAgIH0pO1xuICAgICAgcmV0dXJuIE07XG4gICAgfVxuICAgIFxuICAgIFxuICAgIHZhciBpbml0aWFsaXplUmlkZ2VMaW5lcyA9IGZ1bmN0aW9uKE0pIHtcbiAgICAgIHZhciBuID0gd2lkdGhzLmxlbmd0aDtcbiAgICAgIHZhciBsb2NhbHMgPSBtYXhpbWFzKE1bbiAtIDFdLCB3aWR0aHNbbiAtIDFdKTtcbiAgICAgIHZhciByaWRnZUxpbmVzID0gW107XG4gICAgICBsb2NhbHMuZm9yRWFjaChmdW5jdGlvbihkKSB7XG4gICAgICAgIHZhciBwb2ludCA9IG5ldyBQb2ludChkLngsIGQueSwgd2lkdGhzW24gLSAxXSk7XG4gICAgICAgIHZhciBsaW5lID0gbmV3IFJpZGdlTGluZSgpO1xuICAgICAgICBsaW5lLmFkZChwb2ludCk7XG4gICAgICAgIHJpZGdlTGluZXMucHVzaChsaW5lKTtcbiAgICAgIH0pO1xuICAgICAgcmV0dXJuIHJpZGdlTGluZXM7XG4gICAgfVxuICAgIFxuICAgIHZhciBjb25uZWN0UmlkZ2VMaW5lcyA9IGZ1bmN0aW9uKE0sIHJpZGdlTGluZXMpIHtcbiAgICAgIHZhciBuID0gd2lkdGhzLmxlbmd0aDtcbiAgICAgIGZvciAodmFyIHJvdyA9IG4gLSAyOyByb3cgPj0gMDsgcm93LS0pIHtcbiAgICAgICAgdmFyIGxvY2FscyA9IG1heGltYXMoTVtyb3ddLCB3aWR0aHNbcm93XSk7XG4gICAgICAgIHZhciBhZGRlZExvY2FscyA9IFtdO1xuICAgICAgICBcbiAgICAgICAgLy8gRmluZCBuZWFyZXN0IG5laWdoYm9yIGF0IG5leHQgc2NhbGUgYW5kIGFkZCB0byB0aGUgbGluZVxuICAgICAgICByaWRnZUxpbmVzLmZvckVhY2goZnVuY3Rpb24obGluZSwgaSkge1xuICAgICAgICAgIHZhciB4ID0gbmVhcmVzdE5laWdoYm9yKGxpbmUsIGxvY2Fscywgd2lkdGhzW3Jvd10pO1xuICAgICAgICAgIGxpbmUuYWRkKHggPT09IG51bGwgPyBudWxsIDogbmV3IFBvaW50KHgsIE1bcm93XVt4XSwgd2lkdGhzW3Jvd10pKTtcbiAgICAgICAgICBcbiAgICAgICAgICBpZiAoeCAhPT0gbnVsbCkge1xuICAgICAgICAgICAgYWRkZWRMb2NhbHMucHVzaCh4KTtcbiAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgICAgICBcbiAgICAgICAgLy8gUmVtb3ZlIGxpbmVzIHRoYXQgaGFzIGV4Y2VlZGVkIHRoZSBnYXAgdGhyZXNob2xkXG4gICAgICAgIHJpZGdlTGluZXMgPSByaWRnZUxpbmVzLmZpbHRlcihmdW5jdGlvbihsaW5lKSB7XG4gICAgICAgICAgcmV0dXJuICFsaW5lLmlzRGlzY29ubmVjdGVkKGdhcFRocmVzaG9sZCk7XG4gICAgICAgIH0pO1xuICAgICAgICBcbiAgICAgICAgLy8gQWRkIGFsbCB0aGUgdW5pdGlhbGl6ZWQgcmlkZ2UgbGluZXNcbiAgICAgICAgbG9jYWxzLmZvckVhY2goZnVuY3Rpb24oZCkge1xuICAgICAgICAgIGlmIChhZGRlZExvY2Fscy5pbmRleE9mKGQueCkgIT09IC0xKSByZXR1cm47XG4gICAgICAgICAgXG4gICAgICAgICAgdmFyIHBvaW50ID0gbmV3IFBvaW50KGQueCwgZC55LCB3aWR0aHNbcm93XSk7XG4gICAgICAgICAgdmFyIHJpZGdlTGluZSA9IG5ldyBSaWRnZUxpbmUoKTtcbiAgICAgICAgICByaWRnZUxpbmUuYWRkKHBvaW50KTtcbiAgICAgICAgICByaWRnZUxpbmVzLnB1c2gocmlkZ2VMaW5lKTtcbiAgICAgICAgfSk7XG4gICAgICB9XG4gICAgICByZXR1cm4gcmlkZ2VMaW5lcztcbiAgICB9XG4gICAgXG4gICAgdmFyIGZpbHRlclJpZGdlTGluZXMgPSBmdW5jdGlvbihzaWduYWwsIHJpZGdlTGluZXMpIHtcbiAgICAgIHZhciBzbW9vdGhlciA9IGtlcm5lbCgpXG4gICAgICAgICAgLnN0ZCgxLjApO1xuICAgICAgdmFyIHRyYW5zZm9ybSA9IGNvbnZvbHZlKClcbiAgICAgICAgLmtlcm5lbChzbW9vdGhlcik7XG4gICAgICB2YXIgY29udm9sdXRpb24gPSB0cmFuc2Zvcm0oc2lnbmFsKTtcbiAgICAgICAgXG4gICAgICByaWRnZUxpbmVzID0gcmlkZ2VMaW5lcy5maWx0ZXIoZnVuY3Rpb24obGluZSkge1xuICAgICAgICB2YXIgc25yID0gbGluZS5TTlIoY29udm9sdXRpb24pO1xuICAgICAgICByZXR1cm4gKHNuciA+PSBtaW5TTlIpICYmIChsaW5lLmxlbmd0aCgpID49IG1pbkxpbmVMZW5ndGgpO1xuICAgICAgfSk7XG4gICAgICByZXR1cm4gcmlkZ2VMaW5lc1xuICAgIH1cbiAgICBcbiAgICAvKipcbiAgICAgKiBQaWNrIHRoZSBwb2ludCB3aXRoIHRoZSBoaWdoZXN0IHkgdmFsdWUgd2l0aGluIHRoYXQgcmFuZ2UuXG4gICAgICovXG4gICAgdmFyIHBlYWtzID0gZnVuY3Rpb24oc2lnbmFsLCByaWRnZUxpbmVzKSB7XG4gICAgICB2YXIgcGVha3MgPSByaWRnZUxpbmVzLm1hcChmdW5jdGlvbihsaW5lKSB7XG4gICAgICAgIHZhciBwb2ludHMgPSBsaW5lLnBvaW50cztcbiAgICAgICAgdmFyIG1heFZhbHVlID0gTnVtYmVyLk5FR0FUSVZFX0lORklOSVRZLFxuICAgICAgICAgICAgbWF4UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgIHBvaW50cy5mb3JFYWNoKGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgICAgICAgdmFyIHkgPSBzaWduYWxbcG9pbnQueF07XG4gICAgICAgICAgaWYgKHkgPiBtYXhWYWx1ZSkge1xuICAgICAgICAgICAgbWF4UG9pbnQgPSBwb2ludDtcbiAgICAgICAgICAgIG1heFZhbHVlID0geTtcbiAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgICAgICByZXR1cm4gbWF4UG9pbnQuc2VyaWFsaXplKCk7XG4gICAgICB9KTtcbiAgICAgIHJldHVybiBwZWFrcztcbiAgICB9XG4gICAgXG4gICAgcmV0dXJuIGZpbmRQZWFrcztcbiAgfTtcblxuICB2YXIgdmVyc2lvbiA9IFwiMC4wLjFcIjtcblxuICBleHBvcnRzLnZlcnNpb24gPSB2ZXJzaW9uO1xuICBleHBvcnRzLnJpY2tlciA9IHJpY2tlcjtcbiAgZXhwb3J0cy5jb252b2x2ZSA9IGNvbnZvbHZlO1xuICBleHBvcnRzLmZpbmRQZWFrcyA9IGZpbmRQZWFrcztcblxufSkpOyIsIid1c2Ugc3RyaWN0JztcblxuZnVuY3Rpb24gRkZUKHNpemUpIHtcbiAgdGhpcy5zaXplID0gc2l6ZSB8IDA7XG4gIGlmICh0aGlzLnNpemUgPD0gMSB8fCAodGhpcy5zaXplICYgKHRoaXMuc2l6ZSAtIDEpKSAhPT0gMClcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0ZGVCBzaXplIG11c3QgYmUgYSBwb3dlciBvZiB0d28gYW5kIGJpZ2dlciB0aGFuIDEnKTtcblxuICB0aGlzLl9jc2l6ZSA9IHNpemUgPDwgMTtcblxuICAvLyBOT1RFOiBVc2Ugb2YgYHZhcmAgaXMgaW50ZW50aW9uYWwgZm9yIG9sZCBWOCB2ZXJzaW9uc1xuICB2YXIgdGFibGUgPSBuZXcgQXJyYXkodGhpcy5zaXplICogMik7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgdGFibGUubGVuZ3RoOyBpICs9IDIpIHtcbiAgICBjb25zdCBhbmdsZSA9IE1hdGguUEkgKiBpIC8gdGhpcy5zaXplO1xuICAgIHRhYmxlW2ldID0gTWF0aC5jb3MoYW5nbGUpO1xuICAgIHRhYmxlW2kgKyAxXSA9IC1NYXRoLnNpbihhbmdsZSk7XG4gIH1cbiAgdGhpcy50YWJsZSA9IHRhYmxlO1xuXG4gIC8vIEZpbmQgc2l6ZSdzIHBvd2VyIG9mIHR3b1xuICB2YXIgcG93ZXIgPSAwO1xuICBmb3IgKHZhciB0ID0gMTsgdGhpcy5zaXplID4gdDsgdCA8PD0gMSlcbiAgICBwb3dlcisrO1xuXG4gIC8vIENhbGN1bGF0ZSBpbml0aWFsIHN0ZXAncyB3aWR0aDpcbiAgLy8gICAqIElmIHdlIGFyZSBmdWxsIHJhZGl4LTQgLSBpdCBpcyAyeCBzbWFsbGVyIHRvIGdpdmUgaW5pdGFsIGxlbj04XG4gIC8vICAgKiBPdGhlcndpc2UgaXQgaXMgdGhlIHNhbWUgYXMgYHBvd2VyYCB0byBnaXZlIGxlbj00XG4gIHRoaXMuX3dpZHRoID0gcG93ZXIgJSAyID09PSAwID8gcG93ZXIgLSAxIDogcG93ZXI7XG5cbiAgLy8gUHJlLWNvbXB1dGUgYml0LXJldmVyc2FsIHBhdHRlcm5zXG4gIHRoaXMuX2JpdHJldiA9IG5ldyBBcnJheSgxIDw8IHRoaXMuX3dpZHRoKTtcbiAgZm9yICh2YXIgaiA9IDA7IGogPCB0aGlzLl9iaXRyZXYubGVuZ3RoOyBqKyspIHtcbiAgICB0aGlzLl9iaXRyZXZbal0gPSAwO1xuICAgIGZvciAodmFyIHNoaWZ0ID0gMDsgc2hpZnQgPCB0aGlzLl93aWR0aDsgc2hpZnQgKz0gMikge1xuICAgICAgdmFyIHJldlNoaWZ0ID0gdGhpcy5fd2lkdGggLSBzaGlmdCAtIDI7XG4gICAgICB0aGlzLl9iaXRyZXZbal0gfD0gKChqID4+PiBzaGlmdCkgJiAzKSA8PCByZXZTaGlmdDtcbiAgICB9XG4gIH1cblxuICB0aGlzLl9vdXQgPSBudWxsO1xuICB0aGlzLl9kYXRhID0gbnVsbDtcbiAgdGhpcy5faW52ID0gMDtcbn1cbm1vZHVsZS5leHBvcnRzID0gRkZUO1xuXG5GRlQucHJvdG90eXBlLmZyb21Db21wbGV4QXJyYXkgPSBmdW5jdGlvbiBmcm9tQ29tcGxleEFycmF5KGNvbXBsZXgsIHN0b3JhZ2UpIHtcbiAgdmFyIHJlcyA9IHN0b3JhZ2UgfHwgbmV3IEFycmF5KGNvbXBsZXgubGVuZ3RoID4+PiAxKTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBjb21wbGV4Lmxlbmd0aDsgaSArPSAyKVxuICAgIHJlc1tpID4+PiAxXSA9IGNvbXBsZXhbaV07XG4gIHJldHVybiByZXM7XG59O1xuXG5GRlQucHJvdG90eXBlLmNyZWF0ZUNvbXBsZXhBcnJheSA9IGZ1bmN0aW9uIGNyZWF0ZUNvbXBsZXhBcnJheSgpIHtcbiAgY29uc3QgcmVzID0gbmV3IEFycmF5KHRoaXMuX2NzaXplKTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCByZXMubGVuZ3RoOyBpKyspXG4gICAgcmVzW2ldID0gMDtcbiAgcmV0dXJuIHJlcztcbn07XG5cbkZGVC5wcm90b3R5cGUudG9Db21wbGV4QXJyYXkgPSBmdW5jdGlvbiB0b0NvbXBsZXhBcnJheShpbnB1dCwgc3RvcmFnZSkge1xuICB2YXIgcmVzID0gc3RvcmFnZSB8fCB0aGlzLmNyZWF0ZUNvbXBsZXhBcnJheSgpO1xuICBmb3IgKHZhciBpID0gMDsgaSA8IHJlcy5sZW5ndGg7IGkgKz0gMikge1xuICAgIHJlc1tpXSA9IGlucHV0W2kgPj4+IDFdO1xuICAgIHJlc1tpICsgMV0gPSAwO1xuICB9XG4gIHJldHVybiByZXM7XG59O1xuXG5GRlQucHJvdG90eXBlLmNvbXBsZXRlU3BlY3RydW0gPSBmdW5jdGlvbiBjb21wbGV0ZVNwZWN0cnVtKHNwZWN0cnVtKSB7XG4gIHZhciBzaXplID0gdGhpcy5fY3NpemU7XG4gIHZhciBoYWxmID0gc2l6ZSA+Pj4gMTtcbiAgZm9yICh2YXIgaSA9IDI7IGkgPCBoYWxmOyBpICs9IDIpIHtcbiAgICBzcGVjdHJ1bVtzaXplIC0gaV0gPSBzcGVjdHJ1bVtpXTtcbiAgICBzcGVjdHJ1bVtzaXplIC0gaSArIDFdID0gLXNwZWN0cnVtW2kgKyAxXTtcbiAgfVxufTtcblxuRkZULnByb3RvdHlwZS50cmFuc2Zvcm0gPSBmdW5jdGlvbiB0cmFuc2Zvcm0ob3V0LCBkYXRhKSB7XG4gIGlmIChvdXQgPT09IGRhdGEpXG4gICAgdGhyb3cgbmV3IEVycm9yKCdJbnB1dCBhbmQgb3V0cHV0IGJ1ZmZlcnMgbXVzdCBiZSBkaWZmZXJlbnQnKTtcblxuICB0aGlzLl9vdXQgPSBvdXQ7XG4gIHRoaXMuX2RhdGEgPSBkYXRhO1xuICB0aGlzLl9pbnYgPSAwO1xuICB0aGlzLl90cmFuc2Zvcm00KCk7XG4gIHRoaXMuX291dCA9IG51bGw7XG4gIHRoaXMuX2RhdGEgPSBudWxsO1xufTtcblxuRkZULnByb3RvdHlwZS5yZWFsVHJhbnNmb3JtID0gZnVuY3Rpb24gcmVhbFRyYW5zZm9ybShvdXQsIGRhdGEpIHtcbiAgaWYgKG91dCA9PT0gZGF0YSlcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0lucHV0IGFuZCBvdXRwdXQgYnVmZmVycyBtdXN0IGJlIGRpZmZlcmVudCcpO1xuXG4gIHRoaXMuX291dCA9IG91dDtcbiAgdGhpcy5fZGF0YSA9IGRhdGE7XG4gIHRoaXMuX2ludiA9IDA7XG4gIHRoaXMuX3JlYWxUcmFuc2Zvcm00KCk7XG4gIHRoaXMuX291dCA9IG51bGw7XG4gIHRoaXMuX2RhdGEgPSBudWxsO1xufTtcblxuRkZULnByb3RvdHlwZS5pbnZlcnNlVHJhbnNmb3JtID0gZnVuY3Rpb24gaW52ZXJzZVRyYW5zZm9ybShvdXQsIGRhdGEpIHtcbiAgaWYgKG91dCA9PT0gZGF0YSlcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0lucHV0IGFuZCBvdXRwdXQgYnVmZmVycyBtdXN0IGJlIGRpZmZlcmVudCcpO1xuXG4gIHRoaXMuX291dCA9IG91dDtcbiAgdGhpcy5fZGF0YSA9IGRhdGE7XG4gIHRoaXMuX2ludiA9IDE7XG4gIHRoaXMuX3RyYW5zZm9ybTQoKTtcbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBvdXQubGVuZ3RoOyBpKyspXG4gICAgb3V0W2ldIC89IHRoaXMuc2l6ZTtcbiAgdGhpcy5fb3V0ID0gbnVsbDtcbiAgdGhpcy5fZGF0YSA9IG51bGw7XG59O1xuXG4vLyByYWRpeC00IGltcGxlbWVudGF0aW9uXG4vL1xuLy8gTk9URTogVXNlcyBvZiBgdmFyYCBhcmUgaW50ZW50aW9uYWwgZm9yIG9sZGVyIFY4IHZlcnNpb24gdGhhdCBkbyBub3Rcbi8vIHN1cHBvcnQgYm90aCBgbGV0IGNvbXBvdW5kIGFzc2lnbm1lbnRzYCBhbmQgYGNvbnN0IHBoaWBcbkZGVC5wcm90b3R5cGUuX3RyYW5zZm9ybTQgPSBmdW5jdGlvbiBfdHJhbnNmb3JtNCgpIHtcbiAgdmFyIG91dCA9IHRoaXMuX291dDtcbiAgdmFyIHNpemUgPSB0aGlzLl9jc2l6ZTtcblxuICAvLyBJbml0aWFsIHN0ZXAgKHBlcm11dGUgYW5kIHRyYW5zZm9ybSlcbiAgdmFyIHdpZHRoID0gdGhpcy5fd2lkdGg7XG4gIHZhciBzdGVwID0gMSA8PCB3aWR0aDtcbiAgdmFyIGxlbiA9IChzaXplIC8gc3RlcCkgPDwgMTtcblxuICB2YXIgb3V0T2ZmO1xuICB2YXIgdDtcbiAgdmFyIGJpdHJldiA9IHRoaXMuX2JpdHJldjtcbiAgaWYgKGxlbiA9PT0gNCkge1xuICAgIGZvciAob3V0T2ZmID0gMCwgdCA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4sIHQrKykge1xuICAgICAgY29uc3Qgb2ZmID0gYml0cmV2W3RdO1xuICAgICAgdGhpcy5fc2luZ2xlVHJhbnNmb3JtMihvdXRPZmYsIG9mZiwgc3RlcCk7XG4gICAgfVxuICB9IGVsc2Uge1xuICAgIC8vIGxlbiA9PT0gOFxuICAgIGZvciAob3V0T2ZmID0gMCwgdCA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4sIHQrKykge1xuICAgICAgY29uc3Qgb2ZmID0gYml0cmV2W3RdO1xuICAgICAgdGhpcy5fc2luZ2xlVHJhbnNmb3JtNChvdXRPZmYsIG9mZiwgc3RlcCk7XG4gICAgfVxuICB9XG5cbiAgLy8gTG9vcCB0aHJvdWdoIHN0ZXBzIGluIGRlY3JlYXNpbmcgb3JkZXJcbiAgdmFyIGludiA9IHRoaXMuX2ludiA/IC0xIDogMTtcbiAgdmFyIHRhYmxlID0gdGhpcy50YWJsZTtcbiAgZm9yIChzdGVwID4+PSAyOyBzdGVwID49IDI7IHN0ZXAgPj49IDIpIHtcbiAgICBsZW4gPSAoc2l6ZSAvIHN0ZXApIDw8IDE7XG4gICAgdmFyIHF1YXJ0ZXJMZW4gPSBsZW4gPj4+IDI7XG5cbiAgICAvLyBMb29wIHRocm91Z2ggb2Zmc2V0cyBpbiB0aGUgZGF0YVxuICAgIGZvciAob3V0T2ZmID0gMDsgb3V0T2ZmIDwgc2l6ZTsgb3V0T2ZmICs9IGxlbikge1xuICAgICAgLy8gRnVsbCBjYXNlXG4gICAgICB2YXIgbGltaXQgPSBvdXRPZmYgKyBxdWFydGVyTGVuO1xuICAgICAgZm9yICh2YXIgaSA9IG91dE9mZiwgayA9IDA7IGkgPCBsaW1pdDsgaSArPSAyLCBrICs9IHN0ZXApIHtcbiAgICAgICAgY29uc3QgQSA9IGk7XG4gICAgICAgIGNvbnN0IEIgPSBBICsgcXVhcnRlckxlbjtcbiAgICAgICAgY29uc3QgQyA9IEIgKyBxdWFydGVyTGVuO1xuICAgICAgICBjb25zdCBEID0gQyArIHF1YXJ0ZXJMZW47XG5cbiAgICAgICAgLy8gT3JpZ2luYWwgdmFsdWVzXG4gICAgICAgIGNvbnN0IEFyID0gb3V0W0FdO1xuICAgICAgICBjb25zdCBBaSA9IG91dFtBICsgMV07XG4gICAgICAgIGNvbnN0IEJyID0gb3V0W0JdO1xuICAgICAgICBjb25zdCBCaSA9IG91dFtCICsgMV07XG4gICAgICAgIGNvbnN0IENyID0gb3V0W0NdO1xuICAgICAgICBjb25zdCBDaSA9IG91dFtDICsgMV07XG4gICAgICAgIGNvbnN0IERyID0gb3V0W0RdO1xuICAgICAgICBjb25zdCBEaSA9IG91dFtEICsgMV07XG5cbiAgICAgICAgLy8gTWlkZGxlIHZhbHVlc1xuICAgICAgICBjb25zdCBNQXIgPSBBcjtcbiAgICAgICAgY29uc3QgTUFpID0gQWk7XG5cbiAgICAgICAgY29uc3QgdGFibGVCciA9IHRhYmxlW2tdO1xuICAgICAgICBjb25zdCB0YWJsZUJpID0gaW52ICogdGFibGVbayArIDFdO1xuICAgICAgICBjb25zdCBNQnIgPSBCciAqIHRhYmxlQnIgLSBCaSAqIHRhYmxlQmk7XG4gICAgICAgIGNvbnN0IE1CaSA9IEJyICogdGFibGVCaSArIEJpICogdGFibGVCcjtcblxuICAgICAgICBjb25zdCB0YWJsZUNyID0gdGFibGVbMiAqIGtdO1xuICAgICAgICBjb25zdCB0YWJsZUNpID0gaW52ICogdGFibGVbMiAqIGsgKyAxXTtcbiAgICAgICAgY29uc3QgTUNyID0gQ3IgKiB0YWJsZUNyIC0gQ2kgKiB0YWJsZUNpO1xuICAgICAgICBjb25zdCBNQ2kgPSBDciAqIHRhYmxlQ2kgKyBDaSAqIHRhYmxlQ3I7XG5cbiAgICAgICAgY29uc3QgdGFibGVEciA9IHRhYmxlWzMgKiBrXTtcbiAgICAgICAgY29uc3QgdGFibGVEaSA9IGludiAqIHRhYmxlWzMgKiBrICsgMV07XG4gICAgICAgIGNvbnN0IE1EciA9IERyICogdGFibGVEciAtIERpICogdGFibGVEaTtcbiAgICAgICAgY29uc3QgTURpID0gRHIgKiB0YWJsZURpICsgRGkgKiB0YWJsZURyO1xuXG4gICAgICAgIC8vIFByZS1GaW5hbCB2YWx1ZXNcbiAgICAgICAgY29uc3QgVDByID0gTUFyICsgTUNyO1xuICAgICAgICBjb25zdCBUMGkgPSBNQWkgKyBNQ2k7XG4gICAgICAgIGNvbnN0IFQxciA9IE1BciAtIE1DcjtcbiAgICAgICAgY29uc3QgVDFpID0gTUFpIC0gTUNpO1xuICAgICAgICBjb25zdCBUMnIgPSBNQnIgKyBNRHI7XG4gICAgICAgIGNvbnN0IFQyaSA9IE1CaSArIE1EaTtcbiAgICAgICAgY29uc3QgVDNyID0gaW52ICogKE1CciAtIE1Ecik7XG4gICAgICAgIGNvbnN0IFQzaSA9IGludiAqIChNQmkgLSBNRGkpO1xuXG4gICAgICAgIC8vIEZpbmFsIHZhbHVlc1xuICAgICAgICBjb25zdCBGQXIgPSBUMHIgKyBUMnI7XG4gICAgICAgIGNvbnN0IEZBaSA9IFQwaSArIFQyaTtcblxuICAgICAgICBjb25zdCBGQ3IgPSBUMHIgLSBUMnI7XG4gICAgICAgIGNvbnN0IEZDaSA9IFQwaSAtIFQyaTtcblxuICAgICAgICBjb25zdCBGQnIgPSBUMXIgKyBUM2k7XG4gICAgICAgIGNvbnN0IEZCaSA9IFQxaSAtIFQzcjtcblxuICAgICAgICBjb25zdCBGRHIgPSBUMXIgLSBUM2k7XG4gICAgICAgIGNvbnN0IEZEaSA9IFQxaSArIFQzcjtcblxuICAgICAgICBvdXRbQV0gPSBGQXI7XG4gICAgICAgIG91dFtBICsgMV0gPSBGQWk7XG4gICAgICAgIG91dFtCXSA9IEZCcjtcbiAgICAgICAgb3V0W0IgKyAxXSA9IEZCaTtcbiAgICAgICAgb3V0W0NdID0gRkNyO1xuICAgICAgICBvdXRbQyArIDFdID0gRkNpO1xuICAgICAgICBvdXRbRF0gPSBGRHI7XG4gICAgICAgIG91dFtEICsgMV0gPSBGRGk7XG4gICAgICB9XG4gICAgfVxuICB9XG59O1xuXG4vLyByYWRpeC0yIGltcGxlbWVudGF0aW9uXG4vL1xuLy8gTk9URTogT25seSBjYWxsZWQgZm9yIGxlbj00XG5GRlQucHJvdG90eXBlLl9zaW5nbGVUcmFuc2Zvcm0yID0gZnVuY3Rpb24gX3NpbmdsZVRyYW5zZm9ybTIob3V0T2ZmLCBvZmYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgc3RlcCkge1xuICBjb25zdCBvdXQgPSB0aGlzLl9vdXQ7XG4gIGNvbnN0IGRhdGEgPSB0aGlzLl9kYXRhO1xuXG4gIGNvbnN0IGV2ZW5SID0gZGF0YVtvZmZdO1xuICBjb25zdCBldmVuSSA9IGRhdGFbb2ZmICsgMV07XG4gIGNvbnN0IG9kZFIgPSBkYXRhW29mZiArIHN0ZXBdO1xuICBjb25zdCBvZGRJID0gZGF0YVtvZmYgKyBzdGVwICsgMV07XG5cbiAgY29uc3QgbGVmdFIgPSBldmVuUiArIG9kZFI7XG4gIGNvbnN0IGxlZnRJID0gZXZlbkkgKyBvZGRJO1xuICBjb25zdCByaWdodFIgPSBldmVuUiAtIG9kZFI7XG4gIGNvbnN0IHJpZ2h0SSA9IGV2ZW5JIC0gb2RkSTtcblxuICBvdXRbb3V0T2ZmXSA9IGxlZnRSO1xuICBvdXRbb3V0T2ZmICsgMV0gPSBsZWZ0STtcbiAgb3V0W291dE9mZiArIDJdID0gcmlnaHRSO1xuICBvdXRbb3V0T2ZmICsgM10gPSByaWdodEk7XG59O1xuXG4vLyByYWRpeC00XG4vL1xuLy8gTk9URTogT25seSBjYWxsZWQgZm9yIGxlbj04XG5GRlQucHJvdG90eXBlLl9zaW5nbGVUcmFuc2Zvcm00ID0gZnVuY3Rpb24gX3NpbmdsZVRyYW5zZm9ybTQob3V0T2ZmLCBvZmYsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgc3RlcCkge1xuICBjb25zdCBvdXQgPSB0aGlzLl9vdXQ7XG4gIGNvbnN0IGRhdGEgPSB0aGlzLl9kYXRhO1xuICBjb25zdCBpbnYgPSB0aGlzLl9pbnYgPyAtMSA6IDE7XG4gIGNvbnN0IHN0ZXAyID0gc3RlcCAqIDI7XG4gIGNvbnN0IHN0ZXAzID0gc3RlcCAqIDM7XG5cbiAgLy8gT3JpZ2luYWwgdmFsdWVzXG4gIGNvbnN0IEFyID0gZGF0YVtvZmZdO1xuICBjb25zdCBBaSA9IGRhdGFbb2ZmICsgMV07XG4gIGNvbnN0IEJyID0gZGF0YVtvZmYgKyBzdGVwXTtcbiAgY29uc3QgQmkgPSBkYXRhW29mZiArIHN0ZXAgKyAxXTtcbiAgY29uc3QgQ3IgPSBkYXRhW29mZiArIHN0ZXAyXTtcbiAgY29uc3QgQ2kgPSBkYXRhW29mZiArIHN0ZXAyICsgMV07XG4gIGNvbnN0IERyID0gZGF0YVtvZmYgKyBzdGVwM107XG4gIGNvbnN0IERpID0gZGF0YVtvZmYgKyBzdGVwMyArIDFdO1xuXG4gIC8vIFByZS1GaW5hbCB2YWx1ZXNcbiAgY29uc3QgVDByID0gQXIgKyBDcjtcbiAgY29uc3QgVDBpID0gQWkgKyBDaTtcbiAgY29uc3QgVDFyID0gQXIgLSBDcjtcbiAgY29uc3QgVDFpID0gQWkgLSBDaTtcbiAgY29uc3QgVDJyID0gQnIgKyBEcjtcbiAgY29uc3QgVDJpID0gQmkgKyBEaTtcbiAgY29uc3QgVDNyID0gaW52ICogKEJyIC0gRHIpO1xuICBjb25zdCBUM2kgPSBpbnYgKiAoQmkgLSBEaSk7XG5cbiAgLy8gRmluYWwgdmFsdWVzXG4gIGNvbnN0IEZBciA9IFQwciArIFQycjtcbiAgY29uc3QgRkFpID0gVDBpICsgVDJpO1xuXG4gIGNvbnN0IEZCciA9IFQxciArIFQzaTtcbiAgY29uc3QgRkJpID0gVDFpIC0gVDNyO1xuXG4gIGNvbnN0IEZDciA9IFQwciAtIFQycjtcbiAgY29uc3QgRkNpID0gVDBpIC0gVDJpO1xuXG4gIGNvbnN0IEZEciA9IFQxciAtIFQzaTtcbiAgY29uc3QgRkRpID0gVDFpICsgVDNyO1xuXG4gIG91dFtvdXRPZmZdID0gRkFyO1xuICBvdXRbb3V0T2ZmICsgMV0gPSBGQWk7XG4gIG91dFtvdXRPZmYgKyAyXSA9IEZCcjtcbiAgb3V0W291dE9mZiArIDNdID0gRkJpO1xuICBvdXRbb3V0T2ZmICsgNF0gPSBGQ3I7XG4gIG91dFtvdXRPZmYgKyA1XSA9IEZDaTtcbiAgb3V0W291dE9mZiArIDZdID0gRkRyO1xuICBvdXRbb3V0T2ZmICsgN10gPSBGRGk7XG59O1xuXG4vLyBSZWFsIGlucHV0IHJhZGl4LTQgaW1wbGVtZW50YXRpb25cbkZGVC5wcm90b3R5cGUuX3JlYWxUcmFuc2Zvcm00ID0gZnVuY3Rpb24gX3JlYWxUcmFuc2Zvcm00KCkge1xuICB2YXIgb3V0ID0gdGhpcy5fb3V0O1xuICB2YXIgc2l6ZSA9IHRoaXMuX2NzaXplO1xuXG4gIC8vIEluaXRpYWwgc3RlcCAocGVybXV0ZSBhbmQgdHJhbnNmb3JtKVxuICB2YXIgd2lkdGggPSB0aGlzLl93aWR0aDtcbiAgdmFyIHN0ZXAgPSAxIDw8IHdpZHRoO1xuICB2YXIgbGVuID0gKHNpemUgLyBzdGVwKSA8PCAxO1xuXG4gIHZhciBvdXRPZmY7XG4gIHZhciB0O1xuICB2YXIgYml0cmV2ID0gdGhpcy5fYml0cmV2O1xuICBpZiAobGVuID09PSA0KSB7XG4gICAgZm9yIChvdXRPZmYgPSAwLCB0ID0gMDsgb3V0T2ZmIDwgc2l6ZTsgb3V0T2ZmICs9IGxlbiwgdCsrKSB7XG4gICAgICBjb25zdCBvZmYgPSBiaXRyZXZbdF07XG4gICAgICB0aGlzLl9zaW5nbGVSZWFsVHJhbnNmb3JtMihvdXRPZmYsIG9mZiA+Pj4gMSwgc3RlcCA+Pj4gMSk7XG4gICAgfVxuICB9IGVsc2Uge1xuICAgIC8vIGxlbiA9PT0gOFxuICAgIGZvciAob3V0T2ZmID0gMCwgdCA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4sIHQrKykge1xuICAgICAgY29uc3Qgb2ZmID0gYml0cmV2W3RdO1xuICAgICAgdGhpcy5fc2luZ2xlUmVhbFRyYW5zZm9ybTQob3V0T2ZmLCBvZmYgPj4+IDEsIHN0ZXAgPj4+IDEpO1xuICAgIH1cbiAgfVxuXG4gIC8vIExvb3AgdGhyb3VnaCBzdGVwcyBpbiBkZWNyZWFzaW5nIG9yZGVyXG4gIHZhciBpbnYgPSB0aGlzLl9pbnYgPyAtMSA6IDE7XG4gIHZhciB0YWJsZSA9IHRoaXMudGFibGU7XG4gIGZvciAoc3RlcCA+Pj0gMjsgc3RlcCA+PSAyOyBzdGVwID4+PSAyKSB7XG4gICAgbGVuID0gKHNpemUgLyBzdGVwKSA8PCAxO1xuICAgIHZhciBoYWxmTGVuID0gbGVuID4+PiAxO1xuICAgIHZhciBxdWFydGVyTGVuID0gaGFsZkxlbiA+Pj4gMTtcbiAgICB2YXIgaHF1YXJ0ZXJMZW4gPSBxdWFydGVyTGVuID4+PiAxO1xuXG4gICAgLy8gTG9vcCB0aHJvdWdoIG9mZnNldHMgaW4gdGhlIGRhdGFcbiAgICBmb3IgKG91dE9mZiA9IDA7IG91dE9mZiA8IHNpemU7IG91dE9mZiArPSBsZW4pIHtcbiAgICAgIGZvciAodmFyIGkgPSAwLCBrID0gMDsgaSA8PSBocXVhcnRlckxlbjsgaSArPSAyLCBrICs9IHN0ZXApIHtcbiAgICAgICAgdmFyIEEgPSBvdXRPZmYgKyBpO1xuICAgICAgICB2YXIgQiA9IEEgKyBxdWFydGVyTGVuO1xuICAgICAgICB2YXIgQyA9IEIgKyBxdWFydGVyTGVuO1xuICAgICAgICB2YXIgRCA9IEMgKyBxdWFydGVyTGVuO1xuXG4gICAgICAgIC8vIE9yaWdpbmFsIHZhbHVlc1xuICAgICAgICB2YXIgQXIgPSBvdXRbQV07XG4gICAgICAgIHZhciBBaSA9IG91dFtBICsgMV07XG4gICAgICAgIHZhciBCciA9IG91dFtCXTtcbiAgICAgICAgdmFyIEJpID0gb3V0W0IgKyAxXTtcbiAgICAgICAgdmFyIENyID0gb3V0W0NdO1xuICAgICAgICB2YXIgQ2kgPSBvdXRbQyArIDFdO1xuICAgICAgICB2YXIgRHIgPSBvdXRbRF07XG4gICAgICAgIHZhciBEaSA9IG91dFtEICsgMV07XG5cbiAgICAgICAgLy8gTWlkZGxlIHZhbHVlc1xuICAgICAgICB2YXIgTUFyID0gQXI7XG4gICAgICAgIHZhciBNQWkgPSBBaTtcblxuICAgICAgICB2YXIgdGFibGVCciA9IHRhYmxlW2tdO1xuICAgICAgICB2YXIgdGFibGVCaSA9IGludiAqIHRhYmxlW2sgKyAxXTtcbiAgICAgICAgdmFyIE1CciA9IEJyICogdGFibGVCciAtIEJpICogdGFibGVCaTtcbiAgICAgICAgdmFyIE1CaSA9IEJyICogdGFibGVCaSArIEJpICogdGFibGVCcjtcblxuICAgICAgICB2YXIgdGFibGVDciA9IHRhYmxlWzIgKiBrXTtcbiAgICAgICAgdmFyIHRhYmxlQ2kgPSBpbnYgKiB0YWJsZVsyICogayArIDFdO1xuICAgICAgICB2YXIgTUNyID0gQ3IgKiB0YWJsZUNyIC0gQ2kgKiB0YWJsZUNpO1xuICAgICAgICB2YXIgTUNpID0gQ3IgKiB0YWJsZUNpICsgQ2kgKiB0YWJsZUNyO1xuXG4gICAgICAgIHZhciB0YWJsZURyID0gdGFibGVbMyAqIGtdO1xuICAgICAgICB2YXIgdGFibGVEaSA9IGludiAqIHRhYmxlWzMgKiBrICsgMV07XG4gICAgICAgIHZhciBNRHIgPSBEciAqIHRhYmxlRHIgLSBEaSAqIHRhYmxlRGk7XG4gICAgICAgIHZhciBNRGkgPSBEciAqIHRhYmxlRGkgKyBEaSAqIHRhYmxlRHI7XG5cbiAgICAgICAgLy8gUHJlLUZpbmFsIHZhbHVlc1xuICAgICAgICB2YXIgVDByID0gTUFyICsgTUNyO1xuICAgICAgICB2YXIgVDBpID0gTUFpICsgTUNpO1xuICAgICAgICB2YXIgVDFyID0gTUFyIC0gTUNyO1xuICAgICAgICB2YXIgVDFpID0gTUFpIC0gTUNpO1xuICAgICAgICB2YXIgVDJyID0gTUJyICsgTURyO1xuICAgICAgICB2YXIgVDJpID0gTUJpICsgTURpO1xuICAgICAgICB2YXIgVDNyID0gaW52ICogKE1CciAtIE1Ecik7XG4gICAgICAgIHZhciBUM2kgPSBpbnYgKiAoTUJpIC0gTURpKTtcblxuICAgICAgICAvLyBGaW5hbCB2YWx1ZXNcbiAgICAgICAgdmFyIEZBciA9IFQwciArIFQycjtcbiAgICAgICAgdmFyIEZBaSA9IFQwaSArIFQyaTtcblxuICAgICAgICB2YXIgRkJyID0gVDFyICsgVDNpO1xuICAgICAgICB2YXIgRkJpID0gVDFpIC0gVDNyO1xuXG4gICAgICAgIG91dFtBXSA9IEZBcjtcbiAgICAgICAgb3V0W0EgKyAxXSA9IEZBaTtcbiAgICAgICAgb3V0W0JdID0gRkJyO1xuICAgICAgICBvdXRbQiArIDFdID0gRkJpO1xuXG4gICAgICAgIC8vIE91dHB1dCBmaW5hbCBtaWRkbGUgcG9pbnRcbiAgICAgICAgaWYgKGkgPT09IDApIHtcbiAgICAgICAgICB2YXIgRkNyID0gVDByIC0gVDJyO1xuICAgICAgICAgIHZhciBGQ2kgPSBUMGkgLSBUMmk7XG4gICAgICAgICAgb3V0W0NdID0gRkNyO1xuICAgICAgICAgIG91dFtDICsgMV0gPSBGQ2k7XG4gICAgICAgICAgY29udGludWU7XG4gICAgICAgIH1cblxuICAgICAgICAvLyBEbyBub3Qgb3ZlcndyaXRlIG91cnNlbHZlc1xuICAgICAgICBpZiAoaSA9PT0gaHF1YXJ0ZXJMZW4pXG4gICAgICAgICAgY29udGludWU7XG5cbiAgICAgICAgLy8gSW4gdGhlIGZsaXBwZWQgY2FzZTpcbiAgICAgICAgLy8gTUFpID0gLU1BaVxuICAgICAgICAvLyBNQnI9LU1CaSwgTUJpPS1NQnJcbiAgICAgICAgLy8gTUNyPS1NQ3JcbiAgICAgICAgLy8gTURyPU1EaSwgTURpPU1EclxuICAgICAgICB2YXIgU1QwciA9IFQxcjtcbiAgICAgICAgdmFyIFNUMGkgPSAtVDFpO1xuICAgICAgICB2YXIgU1QxciA9IFQwcjtcbiAgICAgICAgdmFyIFNUMWkgPSAtVDBpO1xuICAgICAgICB2YXIgU1QyciA9IC1pbnYgKiBUM2k7XG4gICAgICAgIHZhciBTVDJpID0gLWludiAqIFQzcjtcbiAgICAgICAgdmFyIFNUM3IgPSAtaW52ICogVDJpO1xuICAgICAgICB2YXIgU1QzaSA9IC1pbnYgKiBUMnI7XG5cbiAgICAgICAgdmFyIFNGQXIgPSBTVDByICsgU1QycjtcbiAgICAgICAgdmFyIFNGQWkgPSBTVDBpICsgU1QyaTtcblxuICAgICAgICB2YXIgU0ZCciA9IFNUMXIgKyBTVDNpO1xuICAgICAgICB2YXIgU0ZCaSA9IFNUMWkgLSBTVDNyO1xuXG4gICAgICAgIHZhciBTQSA9IG91dE9mZiArIHF1YXJ0ZXJMZW4gLSBpO1xuICAgICAgICB2YXIgU0IgPSBvdXRPZmYgKyBoYWxmTGVuIC0gaTtcblxuICAgICAgICBvdXRbU0FdID0gU0ZBcjtcbiAgICAgICAgb3V0W1NBICsgMV0gPSBTRkFpO1xuICAgICAgICBvdXRbU0JdID0gU0ZCcjtcbiAgICAgICAgb3V0W1NCICsgMV0gPSBTRkJpO1xuICAgICAgfVxuICAgIH1cbiAgfVxufTtcblxuLy8gcmFkaXgtMiBpbXBsZW1lbnRhdGlvblxuLy9cbi8vIE5PVEU6IE9ubHkgY2FsbGVkIGZvciBsZW49NFxuRkZULnByb3RvdHlwZS5fc2luZ2xlUmVhbFRyYW5zZm9ybTIgPSBmdW5jdGlvbiBfc2luZ2xlUmVhbFRyYW5zZm9ybTIob3V0T2ZmLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgb2ZmLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgc3RlcCkge1xuICBjb25zdCBvdXQgPSB0aGlzLl9vdXQ7XG4gIGNvbnN0IGRhdGEgPSB0aGlzLl9kYXRhO1xuXG4gIGNvbnN0IGV2ZW5SID0gZGF0YVtvZmZdO1xuICBjb25zdCBvZGRSID0gZGF0YVtvZmYgKyBzdGVwXTtcblxuICBjb25zdCBsZWZ0UiA9IGV2ZW5SICsgb2RkUjtcbiAgY29uc3QgcmlnaHRSID0gZXZlblIgLSBvZGRSO1xuXG4gIG91dFtvdXRPZmZdID0gbGVmdFI7XG4gIG91dFtvdXRPZmYgKyAxXSA9IDA7XG4gIG91dFtvdXRPZmYgKyAyXSA9IHJpZ2h0UjtcbiAgb3V0W291dE9mZiArIDNdID0gMDtcbn07XG5cbi8vIHJhZGl4LTRcbi8vXG4vLyBOT1RFOiBPbmx5IGNhbGxlZCBmb3IgbGVuPThcbkZGVC5wcm90b3R5cGUuX3NpbmdsZVJlYWxUcmFuc2Zvcm00ID0gZnVuY3Rpb24gX3NpbmdsZVJlYWxUcmFuc2Zvcm00KG91dE9mZixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9mZixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHN0ZXApIHtcbiAgY29uc3Qgb3V0ID0gdGhpcy5fb3V0O1xuICBjb25zdCBkYXRhID0gdGhpcy5fZGF0YTtcbiAgY29uc3QgaW52ID0gdGhpcy5faW52ID8gLTEgOiAxO1xuICBjb25zdCBzdGVwMiA9IHN0ZXAgKiAyO1xuICBjb25zdCBzdGVwMyA9IHN0ZXAgKiAzO1xuXG4gIC8vIE9yaWdpbmFsIHZhbHVlc1xuICBjb25zdCBBciA9IGRhdGFbb2ZmXTtcbiAgY29uc3QgQnIgPSBkYXRhW29mZiArIHN0ZXBdO1xuICBjb25zdCBDciA9IGRhdGFbb2ZmICsgc3RlcDJdO1xuICBjb25zdCBEciA9IGRhdGFbb2ZmICsgc3RlcDNdO1xuXG4gIC8vIFByZS1GaW5hbCB2YWx1ZXNcbiAgY29uc3QgVDByID0gQXIgKyBDcjtcbiAgY29uc3QgVDFyID0gQXIgLSBDcjtcbiAgY29uc3QgVDJyID0gQnIgKyBEcjtcbiAgY29uc3QgVDNyID0gaW52ICogKEJyIC0gRHIpO1xuXG4gIC8vIEZpbmFsIHZhbHVlc1xuICBjb25zdCBGQXIgPSBUMHIgKyBUMnI7XG5cbiAgY29uc3QgRkJyID0gVDFyO1xuICBjb25zdCBGQmkgPSAtVDNyO1xuXG4gIGNvbnN0IEZDciA9IFQwciAtIFQycjtcblxuICBjb25zdCBGRHIgPSBUMXI7XG4gIGNvbnN0IEZEaSA9IFQzcjtcblxuICBvdXRbb3V0T2ZmXSA9IEZBcjtcbiAgb3V0W291dE9mZiArIDFdID0gMDtcbiAgb3V0W291dE9mZiArIDJdID0gRkJyO1xuICBvdXRbb3V0T2ZmICsgM10gPSBGQmk7XG4gIG91dFtvdXRPZmYgKyA0XSA9IEZDcjtcbiAgb3V0W291dE9mZiArIDVdID0gMDtcbiAgb3V0W291dE9mZiArIDZdID0gRkRyO1xuICBvdXRbb3V0T2ZmICsgN10gPSBGRGk7XG59O1xuIiwiY29uc3QgZDNfcGVha3MgPSByZXF1aXJlKFwiZDMtcGVha3NcIilcblxuZnVuY3Rpb24gZGJhbXAoZGIpIHtcbiAgICByZXR1cm4gTWF0aC5wb3coMTAsIGRiICogMC4wNSlcbn1cblxuY2xhc3MgTVNEIHtcblxuICAgIGNvbnN0cnVjdG9yKG51bUJpbnMsIGhpc3RvcnlTaXplKXtcbiAgICAgICAgdGhpcy5tc2QgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmYgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeSA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyAqIGhpc3RvcnlTaXplKTtcbiAgICAgICAgdGhpcy5sYXN0TWFnbml0dWRlcyA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubGFzdE1hZ0RpZmYgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmTm9ybWFsaXplID0gMSAvIGhpc3RvcnlTaXplO1xuICAgIH1cblxuICAgIGFkZEZvckJpbihiaW5JbmRleCwgbWFnRGIsIHNtb290aGluZykge1xuICAgICAgICB0aGlzLm1hZ0RpZmZbYmluSW5kZXhdID0gbWFnRGIgLSB0aGlzLmxhc3RNYWduaXR1ZGVzW2JpbkluZGV4XTtcbiAgICAgICAgY29uc3QgbWFnRGlmZkRpZmYgPSBNYXRoLnBvdyh0aGlzLm1hZ0RpZmZbYmluSW5kZXhdIC0gdGhpcy5sYXN0TWFnRGlmZltiaW5JbmRleF0sMik7XG4gICAgICAgIHRoaXMubXNkW2JpbkluZGV4XSArPSAoMS1zbW9vdGhpbmcpICogbWFnRGlmZkRpZmYgKyBzbW9vdGhpbmcgKiB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeVtiaW5JbmRleF0gO1xuICAgICAgICB0aGlzLm1hZ0RpZmZEaWZmSGlzdG9yeVtiaW5JbmRleF0gPSBtYWdEaWZmRGlmZjtcblxuICAgICAgICB0aGlzLmxhc3RNYWdEaWZmW2JpbkluZGV4XSA9IHRoaXMubWFnRGlmZltiaW5JbmRleF07XG4gICAgICAgIHRoaXMubGFzdE1hZ25pdHVkZXNbYmluSW5kZXhdID0gdGhpcy5tYWdEYltiaW5JbmRleF07XG4gICAgfVxuXG59XG5cbmNsYXNzIE1hZ25pdHVkZXNIaXN0b3J5IHtcblxuICAgIGNvbnN0cnVjdG9yKG51bUJpbnMsIGhpc3RvcnlTaXplLCBtaW5QZWFrVGhyPS00MCkge1xuICAgICAgICB0aGlzLm51bUJpbnMgPSBudW1CaW5zO1xuICAgICAgICB0aGlzLnJOdW1CaW5zID0gMS4wIC8gbnVtQmlucztcbiAgICAgICAgdGhpcy5oaXN0b3J5U2l6ZSA9IGhpc3RvcnlTaXplO1xuICAgICAgICB0aGlzLm1hZ0RiID0gbmV3IEZsb2F0MzJBcnJheShudW1CaW5zKTtcbiAgICAgICAgdGhpcy5tYWdTY2FsZSA9IDIgLyB0aGlzLm51bUJpbnM7XG5cbiAgICAgICAgdGhpcy5tc2QgPSBuZXcgTVNEKHRoaXMubnVtQmlucywgaGlzdG9yeVNpemUpO1xuXG4gICAgICAgIHRoaXMubWF4RGIgPSAtMTgwO1xuICAgICAgICB0aGlzLmF2ZXJhZ2VEYiA9IDA7XG4gICAgICAgIHRoaXMuYXZnVGhyID0gMC41O1xuXG4gICAgICAgIHRoaXMubWF4UGVha3MgPSAxMDtcbiAgICAgICAgdGhpcy5uYlBlYWtzID0gMDtcbiAgICAgICAgdGhpcy5taW5QZWFrVGhyID0gbWluUGVha1RocjtcbiAgICAgICAgdGhpcy5wZWFrVGhyID0gMDtcbiAgICAgICAgdGhpcy5wZWFrSW5kZXhlcyA9IG5ldyBJbnQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLnBlYWtIaXN0b3J5U2l6ZSA9IDEwMDtcbiAgICAgICAgdGhpcy5wZWFrSGlzdG9yeSA9IG5ldyBGbG9hdDMyQXJyYXkobnVtQmlucyAqIHRoaXMucGVha0hpc3RvcnlTaXplKTtcbiAgICAgICAgdGhpcy5wZWFrUGVyc2lzdGVuY2UgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLm1pblBlYWtQZXJzaXN0ZW5jZSA9IHRoaXMucGVha0hpc3RvcnlTaXplICogMC43NTtcbiAgICAgICAgdGhpcy5wZWFrV2lkdGggPSAxMTtcblxuICAgICAgICB0aGlzLmZiSW5kZXhlcyA9IG5ldyBJbnQzMkFycmF5KHRoaXMubnVtQmlucyk7XG4gICAgICAgIHRoaXMubmJGYiA9IDA7XG4gICAgICAgIHRoaXMuZmJIaXN0b3J5ID0gbmV3IEludDMyQXJyYXkoMjUpO1xuICAgICAgICB0aGlzLmZiSGlzdG9yeVBvcyA9IDA7XG5cbiAgICAgICAgdGhpcy5tYWdSZWR1Y3Rpb25zID0gbmV3IEZsb2F0MzJBcnJheSh0aGlzLm51bUJpbnMpO1xuICAgICAgICB0aGlzLm1hZ1JlZHVjdGlvbnNBbXAgPSBuZXcgRmxvYXQzMkFycmF5KHRoaXMubnVtQmlucyk7XG4gICAgICAgIHRoaXMubWFnUmVkdWN0aW9uc1VuaXREYiA9IDAuMTtcbiAgICAgICAgdGhpcy5tYWdSZWR1Y3Rpb25zVW5pdEFtcCA9IGRiYW1wKHRoaXMubWFnUmVkdWN0aW9uc1VuaXREYik7XG5cbiAgICAgICAgY29uc3Qgcmlja2VyID0gZDNfcGVha3Mucmlja2VyO1xuICAgICAgICAgIHRoaXMuZmluZFBlYWtzRm4gPSBkM19wZWFrcy5maW5kUGVha3MoKVxuICAgICAgICAgICAgLy8ua2VybmVsKHJpY2tlcilcbiAgICAgICAgICAgIC5nYXBUaHJlc2hvbGQoMTApXG4gICAgICAgICAgICAubWluU05SKDEpXG4gICAgICAgICAgICAud2lkdGhzKFsxLDIsM10pO1xuXG4gICAgICAgIGNvbnNvbGUubG9nKFwiRkZUIGNvbnN0cnVjdGVkXCIpXG4gICAgfVxuXG4gICAgc2hpZnRIaXN0b3J5KCkge1xuICAgICAgICAvL2NvbnN0IGxhc3RNZW1CbG9jayA9ICB0aGlzLm51bUJpbnMgKiAodGhpcy5oaXN0b3J5U2l6ZSAtIDEpO1xuICAgICAgICBjb25zdCBsYXN0UGVha01lbUJsb2NrID0gdGhpcy5udW1CaW5zICogKHRoaXMucGVha0hpc3RvcnlTaXplIC0gMSk7XG4gICAgICAgIGZvcihsZXQgYmluID0gMDsgYmluIDwgdGhpcy5udW1CaW5zOyArK2Jpbikge1xuICAgICAgICAgICAgLy8gdGhpcy5tc2QubXNkW2Jpbl0gLT0gdGhpcy5tc2QubWFnRGlmZkRpZmZIaXN0b3J5W2xhc3RNZW1CbG9jayArIGJpbl07XG4gICAgICAgICAgICB0aGlzLnBlYWtQZXJzaXN0ZW5jZVtiaW5dIC09IHRoaXMucGVha0hpc3RvcnlbbGFzdFBlYWtNZW1CbG9jayArIGJpbl07ICBcbiAgICAgICAgICAgIC8vaWYodGhpcy5tYWdSZWR1Y3Rpb25zW2Jpbl08MCkgdGhpcy5tYWdSZWR1Y3Rpb25zW2Jpbl0gKz0gdGhpcy5tYWdSZWR1Y3Rpb25zVW5pdERiXG4gICAgICAgICAgICAvL2lmKHRoaXMubWFnUmVkdWN0aW9uc0FtcFtiaW5dPDEpIHRoaXMubWFnUmVkdWN0aW9uc1tiaW5dICo9IHRoaXMubWFnUmVkdWN0aW9uc1VuaXRBbXBcbiAgICAgICAgfVxuXG4gICAgICAgIHRoaXMubWF4RGIgPSAtMTgwO1xuICAgICAgICB0aGlzLmF2ZXJhZ2VEYiA9IDA7XG4gICAgICAgIHRoaXMubmJQZWFrcyA9IDA7XG4gICAgICAgIC8vIHRoaXMubXNkLm1hZ0RpZmZEaWZmSGlzdG9yeS5jb3B5V2l0aGluKHRoaXMubnVtQmlucywwKTtcbiAgICAgICAgdGhpcy5wZWFrSGlzdG9yeS5jb3B5V2l0aGluKHRoaXMubnVtQmlucywwKTtcbiAgICB9XG5cbiAgICBhbmFseXNlU3BlY3RydW0oYnVmZmVyLCBtYXgpIHtcbiAgICAgICAgdGhpcy5zaGlmdEhpc3RvcnkoKTtcbiAgICAgICAgYnVmZmVyLmZvckVhY2goKG1hZywgYmluSW5kZXgpPT50aGlzLmFkZEZvckJpbihiaW5JbmRleCwgbWFnKSk7XG4gICAgICAgIHRoaXMuYXZlcmFnZURiID0gdGhpcy5hdmVyYWdlRGIgKiB0aGlzLnJOdW1CaW5zO1xuICAgICAgICB0aGlzLnBlYWtUaHIgPSB0aGlzLmF2ZXJhZ2VEYiArIChtYXggLSB0aGlzLmF2ZXJhZ2VEYikgKiB0aGlzLmF2Z1RocjtcbiAgICAgICAgaWYodGhpcy5wZWFrVGhyIDwgdGhpcy5taW5QZWFrVGhyKSB0aGlzLnBlYWtUaHIgPSB0aGlzLm1pblBlYWtUaHI7XG5cbiAgICAgICAgdGhpcy5maW5kUGVha3MoKTtcbiAgICAgICAgdGhpcy5maW5kRmVlZGJhY2tDYW5kaWRhdGVzKCk7XG4gICAgfVxuXG4gICAgYWRkRm9yQmluKGJpbkluZGV4LCBtYWdEYikge1xuICAgICAgICB0aGlzLm1hZ0RiW2JpbkluZGV4XSA9IG1hZ0RiO1xuICAgICAgICAvLyB0aGlzLm1zZC5hZGRGb3JCaW4oYmluSW5kZXgsIG1hZ0RiLCB0aGlzLnNtb290aGluZylcbiAgICAgICAgaWYobWFnRGIgPiB0aGlzLm1heERiKSB0aGlzLm1heERiID0gbWFnRGI7XG4gICAgICAgIHRoaXMuYXZlcmFnZURiICs9IG1hZ0RiO1xuICAgICAgICB0aGlzLnBlYWtIaXN0b3J5W2JpbkluZGV4XSA9IDA7XG4gICAgICAgIC8vdGhpcy5jaGVja1BlYWsoYmluSW5kZXgpO1xuICAgICAgICB0aGlzLmNvcnJlY3RNYWduaXR1ZGUoYmluSW5kZXgpXG4gICAgfVxuXG4gICAgZmluZFBlYWtzKCkge1xuICAgICAgICB2YXIgcGVha3MgPSB0aGlzLmZpbmRQZWFrc0ZuKHRoaXMubWFnRGIpO1xuICAgICAgICBwZWFrcy5mb3JFYWNoKHA9PnRoaXMuYWRkUGVhayhwLmluZGV4KSlcbiAgICAgICAgLy8gY29uc29sZS5sb2cocGVha3MpO1xuICAgIH0gICBcbiAgICBcbiAgICAvKmZpbmRQZWFrc19jdXN0b20oKXtcbiAgICAgICAgdmFyIG1hZ3MgPSB0aGlzLm1hZ0RiO1xuICAgICAgICBjb25zdCB0aHIgPSB0aGlzLnBlYWtUaHI7XG4gICAgICAgIGxldCBibG9ja01heCA9IC1JbmZpbml0eTtcbiAgICAgICAgbGV0IGJsb2NrTWF4QmluID0gLTE7XG4gICAgICAgIGxldCBsX2JpbiA9IDA7IGxldCByX2JpbiA9IDA7IFxuXG4gICAgICAgIGNvbnN0IGNoZWNrQmxvY2sgPSAoKT0+e1xuICAgICAgICAgICAgcl9iaW4gPSBsX2JpbjtcbiAgICAgICAgICAgIGJsb2NrTWF4ID0gLUluZmluaXR5O1xuICAgICAgICAgICAgYmxvY2tNYXhCaW4gPSAtMTtcbiAgICAgICAgICAgIHdoaWxlKHJfYmluIDwgbF9iaW4gKyB0aGlzLnBlYWtXaWR0aCAqIDIpIHtcbiAgICAgICAgICAgICAgICBpZihtYWdzW3JfYmluXSA+PSB0aHIgJiYgbWFnc1tyX2Jpbl0gPiBibG9ja01heCkgeyBibG9ja01heCA9IG1hZ3Nbcl9iaW5dOyBibG9ja01heEJpbiA9IHJfYmluIH1cbiAgICAgICAgICAgICAgICByX2JpbisrO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgY2hlY2tCbG9jaygpO1xuICAgICAgICBsZXQgYmluID0gcl9iaW4gLSB0aGlzLnBlYWtXaWR0aDtcbiAgICAgICAgaWYoYmxvY2tNYXhCaW4gPT0gYmluKSB0aGlzLmFkZFBlYWsoYmluKTtcbiAgICAgICAgcl9iaW4rKzsgbF9iaW4rKzsgYmluKys7XG5cbiAgICAgICAgd2hpbGUocl9iaW4gPCB0aGlzLm51bUJpbnMpe1xuICAgICAgICAgICAgaWYobWFnc1tyX2Jpbl0gPj0gdGhyICYmIG1hZ3Nbcl9iaW5dID4gYmxvY2tNYXgpIHsgXG4gICAgICAgICAgICAgICAgYmxvY2tNYXggPSBtYWdzW3JfYmluXTsgYmxvY2tNYXhCaW4gPSByX2JpblxuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBpZihibG9ja01heEJpbiA+IDAgJiYgYmxvY2tNYXhCaW4gPCBsX2JpbikgY2hlY2tCbG9jaygpO1xuICAgICAgICAgICAgICAgIGlmKGJsb2NrTWF4QmluID09IGJpbikgdGhpcy5hZGRQZWFrKGJpbik7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByX2JpbisrOyBsX2JpbisrOyBiaW4rKztcbiAgICAgICAgfVxuICAgIH0qL1xuXG4gICAgYWRkUGVhayhiaW5JbmRleCl7XG4gICAgICAgIHRoaXMucGVha0luZGV4ZXNbdGhpcy5uYlBlYWtzXSA9IGJpbkluZGV4O1xuICAgICAgICB0aGlzLnBlYWtIaXN0b3J5W2JpbkluZGV4XSA9IDE7XG4gICAgICAgIHRoaXMucGVha1BlcnNpc3RlbmNlW2JpbkluZGV4XSsrO1xuICAgICAgICB0aGlzLm5iUGVha3MrKztcbiAgICB9XG4gICAgXG4gICAgZmluZEZlZWRiYWNrQ2FuZGlkYXRlcygpIHtcbiAgICAgICAgdGhpcy5uYkZiID0gMDtcbiAgICAgICAgbGV0IG1pbk1zZEJpbiA9IC0xO1xuICAgICAgICBsZXQgbWluTXNkID0gSW5maW5pdHk7XG5cbiAgICAgICAgLy9jb25zb2xlLmxvZyh0aGlzLnBlYWtIaXN0b3J5KVxuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHRoaXMubmJQZWFrczsgKytpKSB7XG4gICAgICAgICAgICBjb25zdCBiaW4gPSB0aGlzLnBlYWtJbmRleGVzW2ldO1xuICAgICAgICAgICAgaWYodGhpcy5wZWFrUGVyc2lzdGVuY2VbYmluXSA+PSB0aGlzLm1pblBlYWtQZXJzaXN0ZW5jZSkge1xuICAgICAgICAgICAgLyppZih0aGlzLm1zZFtiaW5dIDwgbWluTXNkICYmIHRoaXMuc21vb3RoZWRNYWdEYltiaW4tMV0gPCB0aGlzLnNtb290aGVkTWFnRGJbYmluXSAmJiB0aGlzLnNtb290aGVkTWFnRGJbYmluKzFdIDwgdGhpcy5zbW9vdGhlZE1hZ0RiW2Jpbl0pIHtcbiAgICAgICAgICAgICAgICBtaW5Nc2QgPSB0aGlzLm1zZFtiaW5dO1xuICAgICAgICAgICAgICAgIG1pbk1zZEJpbiA9IGJpbjtcbiAgICAgICAgICAgIH0qL1xuICAgICAgICAgICAgdGhpcy5mYkluZGV4ZXNbdGhpcy5uYkZiKytdID0gYmluO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgLy90aGlzLmNvcnJlY3RNYWduaXR1ZGUoYmluKVxuXG4gICAgICAgIH1cblxuICAgICAgICBcbiAgICAgICAgLyppZihtaW5Nc2RCaW4gPj0gMCAmJiBtaW5Nc2QgPD0gMC4wMSkge1xuICAgICAgICAgICAgdGhpcy5mYkhpc3RvcnlbdGhpcy5mYkhpc3RvcnlQb3MrK10gPSBtaW5Nc2RCaW47XG4gICAgICAgICAgICBpZih0aGlzLmZiSGlzdG9yeS5yZWR1Y2UoKHQsdik9PnY9PW1pbk1zZEJpbj90KzE6dCwwKSA+PSB0aGlzLmZiSGlzdG9yeS5sZW5ndGggKiAwLjUpIHtcbiAgICAgICAgICAgICAgICB0aGlzLmZiSW5kZXhlc1t0aGlzLm5iRmIrK10gPSBtaW5Nc2RCaW47XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0aGlzLmZiSGlzdG9yeVt0aGlzLmZiSGlzdG9yeVBvcysrXSA9IC0xO1xuICAgICAgICB9XG4gICAgICAgIHRoaXMuZmJIaXN0b3J5LmZpbHRlcigodixpLGEpPT52ID49IDAgJiYgYS5pbmRleE9mKHYpPT09aSkuZm9yRWFjaChiaW49PntcbiAgICAgICAgICAgIHRoaXMuZmJJbmRleGVzW3RoaXMubmJGYl0gPSBiaW47XG4gICAgICAgICAgICB0aGlzLm5iRmIrKztcbiAgICAgICAgfSlcbiAgICAgICAgaWYodGhpcy5mYkhpc3RvcnlQb3MgPiB0aGlzLmZiSGlzdG9yeS5sZW5ndGgpIHRoaXMuZmJIaXN0b3J5UG9zID0gMDtcbiAgICAgICAgKi9cbiAgICB9XG5cbiAgICBjb3JyZWN0TWFnbml0dWRlKGJpbikge1xuICAgICAgICAvL3RoaXMubWFnUmVkdWN0aW9uc1tiaW5dICs9ICh0aGlzLnBlYWtUaHIgLSB0aGlzLm1hZ0RiW2Jpbl0pIC8gMlxuICAgICAgICBjb25zdCBjb3JyZWN0aW9uID0gKHRoaXMucGVha1RociAtIHRoaXMubWFnRGJbYmluXSk7XG4gICAgICAgIGlmKGNvcnJlY3Rpb24gPD0gMCApe1xuICAgICAgICAgICAgdGhpcy5tYWdSZWR1Y3Rpb25zW2Jpbl0gPSBjb3JyZWN0aW9uO1xuICAgICAgICAgICAgdGhpcy5tYWdSZWR1Y3Rpb25zQW1wW2Jpbl0gPSBkYmFtcCh0aGlzLm1hZ1JlZHVjdGlvbnNbYmluXSlcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHRoaXMubWFnUmVkdWN0aW9uc1tiaW5dID0gMDtcbiAgICAgICAgICAgIHRoaXMubWFnUmVkdWN0aW9uc0FtcFtiaW5dID0gMTtcbiAgICAgICAgfVxuICAgIH1cblxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHsgTVNELCBNYWduaXR1ZGVzSGlzdG9yeSB9IiwiY2xhc3MgQXV0b0dhaW4ge1xuICAgIGNvbnN0cnVjdG9yKGF1ZGlvQ29udGV4dCkge1xuICAgICAgICB0aGlzLmdhaW5JbmNyZW1lbnQgPSAwLjAwMTtcbiAgICAgICAgdGhpcy5nYWluRGVjcmVtZW50ID0gMC4wMDE7XG4gICAgICAgIHRoaXMubWluVm9sdW1lID0gLTEwO1xuICAgICAgICB0aGlzLm1heFZvbHVtZSA9IC0xMDtcbiAgICAgICAgdGhpcy5tYXhHYWluID0gMWU1O1xuICAgICAgICB0aGlzLm1pbkdhaW4gPSAtMWU1O1xuXG4gICAgICAgIHRoaXMuZ2FpbiA9IGF1ZGlvQ29udGV4dC5jcmVhdGVHYWluKCk7XG4gICAgICAgIHRoaXMubGltaXRlciA9IGF1ZGlvQ29udGV4dC5jcmVhdGVEeW5hbWljc0NvbXByZXNzb3IoKVxuICAgICAgICB0aGlzLmxpbWl0ZXIudGhyZXNob2xkLnZhbHVlID0gLTQwXG4gICAgICAgIHRoaXMubGltaXRlci5yYXRpby52YWx1ZSA9IDIwXG4gICAgICAgIHRoaXMubGltaXRlci5hdHRhY2sudmFsdWUgPSAwLjAwMVxuICAgICAgICB0aGlzLmxpbWl0ZXIucmVsZWFzZS52YWx1ZSA9IDAuMDFcblxuICAgICAgICB0aGlzLm5vdyA9ICgpPT5hdWRpb0NvbnRleHQuY3VycmVudFRpbWVcbiAgICB9XG5cbiAgICB1cGRhdGVHYWluKG1heERiKSB7XG4gICAgICAgIGlmKG1heERiIDw9IC0xODApIHJldHVyblxuICAgICAgICBjb25zdCBjdXJyZW50R2FpbiA9IDIwICogTWF0aC5sb2cxMCh0aGlzLmdhaW4uZ2Fpbi52YWx1ZSk7XG4gICAgICAgIHRoaXMuZ2Fpbi5nYWluLmNhbmNlbEFuZEhvbGRBdFRpbWUodGhpcy5ub3coKSlcblxuICAgICAgICBsZXQgbmV4dEdhaW4gPSBjdXJyZW50R2FpbjtcbiAgICAgICAgaWYobWF4RGIgPiB0aGlzLm1heFZvbHVtZSkge1xuICAgICAgICAgICAgbmV4dEdhaW4gKz0gKHRoaXMubWluVm9sdW1lIC0gbWF4RGIpICogdGhpcy5nYWluRGVjcmVtZW50XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBuZXh0R2FpbiArPSAodGhpcy5tYXhWb2x1bWUgLSBtYXhEYikgKiB0aGlzLmdhaW5JbmNyZW1lbnRcbiAgICAgICAgfVxuICAgICAgICBpZihuZXh0R2FpbiAhPSBjdXJyZW50R2Fpbikge1xuICAgICAgICAgICAgaWYodGhpcy5pc0dhaW5WYWxpZChuZXh0R2FpbikpIHtcbiAgICAgICAgICAgICAgICBuZXh0R2FpbiA9IE1hdGgucG93KDEwLCBuZXh0R2FpbiowLjA1KTtcbiAgICAgICAgICAgICAgICB0aGlzLmdhaW4uZ2Fpbi5saW5lYXJSYW1wVG9WYWx1ZUF0VGltZShuZXh0R2FpbiwgdGhpcy5ub3coKSlcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cblxuICAgIGlzR2FpblZhbGlkKGdhaW4pIHtcbiAgICAgICAgcmV0dXJuIGdhaW4gPiB0aGlzLm1pbkdhaW4gJiYgZ2FpbiA8IHRoaXMubWF4R2FpbjtcbiAgICB9XG5cbiAgICBnZXRDdXJyZW50VmFsdWUoKSB7IHJldHVybiB0aGlzLmdhaW4uZ2Fpbi52YWx1ZSB9XG4gICAgY29ubmVjdChkZXN0aW5hdGlvbikgeyByZXR1cm4gdGhpcy5nYWluLmNvbm5lY3QodGhpcy5saW1pdGVyKS5jb25uZWN0KGRlc3RpbmF0aW9uKSB9XG59XG5cbm1vZHVsZS5leHBvcnRzID0geyBBdXRvR2FpbiB9IiwiY29uc3QgRkZUID0gcmVxdWlyZSgnZmZ0LmpzJyk7XG5cbmNsYXNzIEZGVENvbnZGaWx0ZXIge1xuICAgIGNvbnN0cnVjdG9yKGF1ZGlvQ29udGV4dCwgbnVtQmlucykge1xuICAgICAgICB0aGlzLm51bUJpbnMgPSBudW1CaW5zO1xuICAgICAgICB0aGlzLmlyQnVmZmVyID0gYXVkaW9Db250ZXh0LmNyZWF0ZUJ1ZmZlcigyLCBudW1CaW5zLCBhdWRpb0NvbnRleHQuc2FtcGxlUmF0ZSlcbiAgICAgICAgdGhpcy5pckJ1ZmZlci5nZXRDaGFubmVsRGF0YSgwKVswXSA9ICgxKVxuXG4gICAgICAgIHRoaXMuZmZ0ID0gbmV3IEZGVChudW1CaW5zICogMik7XG4gICAgICAgIHRoaXMuZnJlcUNvbXBsZXhCdWZmZXIgPSB0aGlzLmZmdC5jcmVhdGVDb21wbGV4QXJyYXkoKTtcbiAgICAgICAgdGhpcy50aW1lQ29tcGxleEJ1ZmZlciA9IHRoaXMuZmZ0LmNyZWF0ZUNvbXBsZXhBcnJheSgpO1xuXG4gICAgICAgIHRoaXMubm9kZSA9IGF1ZGlvQ29udGV4dC5jcmVhdGVDb252b2x2ZXIoKTtcbiAgICAgICAgdGhpcy5ub2RlLm5vcm1hbGl6ZSA9IGZhbHNlXG4gICAgICAgIHRoaXMubm9kZS5idWZmZXIgPSB0aGlzLmlyQnVmZmVyO1xuICAgIH1cblxuICAgIHVwZGF0ZUtlcm5lbChyZWR1Y3Rpb25TcGVjdHJ1bSkge1xuICAgICAgICB0aGlzLmZyZXFDb21wbGV4QnVmZmVyLmZpbGwoMCk7XG4gICAgICAgIGZvcihsZXQgaSA9IDAsIGogPSAwOyBpIDwgdGhpcy5udW1CaW5zOyArK2ksIGorPTIpIHtcbiAgICAgICAgICAgIHRoaXMuZnJlcUNvbXBsZXhCdWZmZXJbal0gPSByZWR1Y3Rpb25TcGVjdHJ1bVtpXTtcbiAgICAgICAgICAgIHRoaXMuZnJlcUNvbXBsZXhCdWZmZXJbaisxXSA9IHJlZHVjdGlvblNwZWN0cnVtW2ldO1xuICAgICAgICB9O1xuICAgICAgICB0aGlzLmZmdC5jb21wbGV0ZVNwZWN0cnVtKHRoaXMuZnJlcUNvbXBsZXhCdWZmZXIpXG4gICAgICAgIHRoaXMuZmZ0LmludmVyc2VUcmFuc2Zvcm0odGhpcy50aW1lQ29tcGxleEJ1ZmZlciwgdGhpcy5mcmVxQ29tcGxleEJ1ZmZlcik7XG4gICAgICAgIHRoaXMuZmZ0LmZyb21Db21wbGV4QXJyYXkodGhpcy50aW1lQ29tcGxleEJ1ZmZlciwgdGhpcy5pckJ1ZmZlci5nZXRDaGFubmVsRGF0YSgwKSk7XG4gICAgICAgIHRoaXMubm9kZS5idWZmZXIgPSB0aGlzLmlyQnVmZmVyXG4gICAgfVxuXG4gICAgY29ubmVjdChkZXN0aW5hdGlvbikgeyByZXR1cm4gdGhpcy5ub2RlLmNvbm5lY3QoZGVzdGluYXRpb24pIH1cbn1cblxubW9kdWxlLmV4cG9ydHMgPSB7IEZGVENvbnZGaWx0ZXIgfSIsIm1vZHVsZS5leHBvcnRzID0ge1xuICAgIC4uLnJlcXVpcmUoXCIuL2FuYWwuanNcIiksXG4gICAgLi4ucmVxdWlyZShcIi4vYXV0b0dhaW4uanNcIiksXG4gICAgLi4ucmVxdWlyZShcIi4vZmZ0Q29udkZpbHRlci5qc1wiKSxcbiAgICAuLi5yZXF1aXJlKFwiLi9ub3RjaEZpbHRlckJhbmsuanNcIiksXG59IiwiY29uc3Qge2NhbGN1bGF0ZVFzLCBjaHJvbWF0aWNGaWx0ZXJ9ID0gcmVxdWlyZSgnLi4vc2NhbGVzLmpzJylcblxuY2xhc3MgTm90Y2hGaWx0ZXJCYW5rIHtcbiAgICBjb25zdHJ1Y3RvcihhdWRpb0NvbnRleHQsIHNjYWxlLCBxKSB7XG4gICAgICAgIHRoaXMuYXVkaW9Db250ZXh0ID0gYXVkaW9Db250ZXh0XG4gICAgICAgIGlmKCFzY2FsZSkge1xuICAgICAgICAgICAgY29uc3Qgc2NhbGVPYmogPSBjaHJvbWF0aWNGaWx0ZXIoKTtcbiAgICAgICAgICAgIHNjYWxlID0gc2NhbGVPYmouc2NhbGU7XG4gICAgICAgICAgICBxID0gc2NhbGVPYmoucXM7XG4gICAgICAgIH1cbiAgICAgICAgdGhpcy5zY2FsZSA9IHNjYWxlO1xuICAgICAgICB0aGlzLmZpbHRlcnMgPSBbXTtcbiAgICAgICAgdGhpcy5pbk1ldGVycyA9IFtdO1xuICAgICAgICB0aGlzLmdhaW5zID0gW107XG4gICAgICAgIHRoaXMub3V0TWV0ZXJzID0gW107XG5cbiAgICAgICAgdGhpcy5xID0gcSB8fCBjYWxjdWxhdGVRcyh0aGlzLnNjYWxlKTsgLy9jb25zb2xlLmxvZyhxKTtcblxuICAgICAgICB0aGlzLm1ha2VGaWx0ZXJzRnJvbVNjYWxlKHRoaXMuc2NhbGUsIHRoaXMucSk7XG4gICAgfVxuXG4gICAgbWFrZUZpbHRlcnNGcm9tU2NhbGUoc2NhbGUsIHEpIHtcbiAgICAgICAgdGhpcy5zY2FsZS5mb3JFYWNoKChmcmVxLGkpPT57XG4gICAgICAgICAgY29uc3QgZmlsdGVyID0gdGhpcy5hdWRpb0NvbnRleHQuY3JlYXRlQmlxdWFkRmlsdGVyKCk7XG4gICAgICAgICAgZmlsdGVyLnR5cGUgPSBcIm5vdGNoXCI7IGZpbHRlci5mcmVxdWVuY3kuc2V0VmFsdWVBdFRpbWUoZnJlcSwgdGhpcy5hdWRpb0NvbnRleHQuY3VycmVudFRpbWUpO1xuICAgICAgICAgIGZpbHRlci5nYWluLnZhbHVlID0gMTtcbiAgICAgICAgICAvLyBUT0RPOiBuZWVkIHRvIGNhbGN1bGF0ZSBRIG1vcmUgYWNjdXJhdGVseVxuICAgICAgICAgIGZpbHRlci5RLnZhbHVlID0gdHlwZW9mIHEgPT0gJ29iamVjdCcgPyBxW2ldIDogcSA7XG4gICAgICAgICAgdGhpcy5maWx0ZXJzLnB1c2goZmlsdGVyKTtcbiAgICAgICAgfSlcblxuICAgICAgICB0aGlzLmZpbHRlcnMuZm9yRWFjaCgoZixuKSA9PiB7IGlmKG48dGhpcy5maWx0ZXJzLmxlbmd0aC0xKWYuY29ubmVjdCh0aGlzLmZpbHRlcnNbbisxXSkgfSk7XG4gICAgfVxuXG4gICAgZ2V0RnJlcXVlbmN5UmVzcG9uc2UoZnJlcXMsIHJlc3BvbnNlKSB7XG4gICAgICAgIGNvbnN0IHRtcCA9IG5ldyBGbG9hdDMyQXJyYXkoZnJlcXMubGVuZ3RoKTtcbiAgICAgICAgY29uc3QgXyA9IG5ldyBGbG9hdDMyQXJyYXkoZnJlcXMubGVuZ3RoKTtcbiAgICAgICAgdGhpcy5maWx0ZXJzLmZvckVhY2goZj0+e1xuICAgICAgICAgICAgZi5nZXRGcmVxdWVuY3lSZXNwb25zZShmcmVxcyx0bXAsXylcbiAgICAgICAgICAgIHRtcC5mb3JFYWNoKCh0LG4pPT5yZXNwb25zZVtuXSo9dClcbiAgICAgICAgfSlcbiAgICB9XG5cbiAgICBjb25uZWN0SW5wdXQoaW5wdXQpIHtcbiAgICAgICAgaW5wdXQuY29ubmVjdCh0aGlzLmZpbHRlcnNbMF0pXG4gICAgfVxuICAgIGNvbm5lY3Qob3V0cHV0KSB7XG4gICAgICAgIHRoaXMuZmlsdGVyc1t0aGlzLmZpbHRlcnMubGVuZ3RoLTFdLmNvbm5lY3Qob3V0cHV0KVxuICAgICAgICByZXR1cm4gb3V0cHV0XG4gICAgfVxuXG59XG5cbm1vZHVsZS5leHBvcnRzID0geyBOb3RjaEZpbHRlckJhbmsgfSIsImNvbnN0IFNxdWlkYmFja0dyYXBoID0gcmVxdWlyZSgnLi9ncmFwaC5qcycpXG5jb25zdCBHcmFwaCA9IG5ldyBTcXVpZGJhY2tHcmFwaCgpO1xuY29uc3QgeyBcbiAgICBBdXRvR2FpbiwgXG4gICAgRkZUQ29udkZpbHRlcixcbiAgICBNYWduaXR1ZGVzSGlzdG9yeVxufSA9IHJlcXVpcmUoXCIuL2F1ZGlvLXByb2Nlc3NvcnMvaW5kZXguanNcIilcblxuY2xhc3MgU3F1aWRiYWNrRkZUUHJvY2VzcyB7XG5cbiAgICBjb25zdHJ1Y3RvcihhdWRpb0NvbnRleHQsIGZmdFNpemUgPSAxMDI0KSB7XG4gICAgICAgIHRoaXMuYXVkaW9Db250ZXh0ID0gYXVkaW9Db250ZXh0O1xuICAgIH1cblxuICAgIGFzeW5jIHN0YXJ0KCkge1xuICAgICAgICB0aGlzLmluaXRGRlROb2RlKCk7XG4gICAgICAgIHRoaXMuaW5pdEJ1ZmZlcnMoKTtcbiAgICAgICAgdGhpcy5pbml0Q29udigpO1xuICAgICAgICB0aGlzLmluaXRFeHRyZW1lc0ZpbHRlcnMoKTtcbiAgICAgICAgdGhpcy5pbml0R2FpbigpO1xuICAgICAgICBhd2FpdCB0aGlzLmluaXRJbnB1dCgpO1xuICAgICAgICBpZih0aGlzLmlucHV0KSB0aGlzLmNvbm5lY3RBbGwoKTtcbiAgICAgICAgcmVxdWVzdEFuaW1hdGlvbkZyYW1lKCgpPT50aGlzLnVwZGF0ZSgpKTtcbiAgICB9XG5cbiAgICBpbml0R2FpbigpIHtcbiAgICAgICAgdGhpcy5hdXRvR2FpbiA9IG5ldyBBdXRvR2Fpbih0aGlzLmF1ZGlvQ29udGV4dClcbiAgICB9XG5cbiAgICBpbml0QnVmZmVycygpIHtcbiAgICAgICAgY29uc3QgbnVtQmlucyA9IHRoaXMuZmZ0Tm9kZS5mcmVxdWVuY3lCaW5Db3VudDtcbiAgICAgICAgdGhpcy5udW1CaW5zID0gbnVtQmlucztcbiAgICAgICAgdGhpcy5iaW5Ub0ZyZXEgPSB0aGlzLmF1ZGlvQ29udGV4dC5zYW1wbGVSYXRlIC8gMiAvIG51bUJpbnM7XG4gICAgICAgIHRoaXMuZmZ0QnVmZmVyID0gbmV3IEZsb2F0MzJBcnJheShudW1CaW5zKTtcbiAgICAgICAgdGhpcy5vdXRGZnRCdWZmZXIgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLmZiQnVmZmVyID0gbmV3IEludDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubnVtRmIgPSAwO1xuXG4gICAgICAgIHRoaXMubWF4RGIgPSAtMTgwO1xuXG4gICAgICAgIHRoaXMuYW5hbCA9IG5ldyBNYWduaXR1ZGVzSGlzdG9yeShudW1CaW5zLCAxMCwgLTYwKTtcbiAgICB9XG5cbiAgICBpbml0Q29udigpIHtcbiAgICAgICAgdGhpcy5jb252ID0gbmV3IEZGVENvbnZGaWx0ZXIodGhpcy5hdWRpb0NvbnRleHQsIHRoaXMubnVtQmlucylcbiAgICB9XG5cbiAgICBpbml0RXh0cmVtZXNGaWx0ZXJzKCkge1xuICAgICAgICBjb25zdCBoaXBhc3MgPSB0aGlzLmF1ZGlvQ29udGV4dC5jcmVhdGVCaXF1YWRGaWx0ZXIoKTtcbiAgICAgICAgaGlwYXNzLnR5cGUgPSBcImhpZ2hwYXNzXCI7IGhpcGFzcy5mcmVxdWVuY3kuc2V0VmFsdWVBdFRpbWUodGhpcy5iaW5Ub0ZyZXEgKiAyLCB0aGlzLm5vdygpKTtcbiAgICAgICAgaGlwYXNzLlEuc2V0VmFsdWVBdFRpbWUoMC43MDcsIHRoaXMuYXVkaW9Db250ZXh0LmN1cnJlbnRUaW1lKTtcblxuICAgICAgICBjb25zdCBsb3Bhc3MgPSB0aGlzLmF1ZGlvQ29udGV4dC5jcmVhdGVCaXF1YWRGaWx0ZXIoKTtcbiAgICAgICAgbG9wYXNzLnR5cGUgPSBcImxvd3Bhc3NcIjsgbG9wYXNzLmZyZXF1ZW5jeS5zZXRWYWx1ZUF0VGltZSgxNTAwMCwgdGhpcy5ub3coKSk7XG4gICAgICAgIGxvcGFzcy5RLnNldFZhbHVlQXRUaW1lKDAuNzA3LCB0aGlzLmF1ZGlvQ29udGV4dC5jdXJyZW50VGltZSk7XG5cbiAgICAgICAgdGhpcy5oaXBhc3MgPSBoaXBhc3M7IHRoaXMubG9wYXNzID0gbG9wYXNzO1xuICAgIH1cblxuICAgIGluaXRGRlROb2RlKCkge1xuICAgICAgICB0aGlzLmZmdE5vZGUgPSB0aGlzLmF1ZGlvQ29udGV4dC5jcmVhdGVBbmFseXNlcigpO1xuICAgICAgICB0aGlzLmZmdE5vZGUuZmZ0U2l6ZSA9IDUxMjtcbiAgICAgICAgdGhpcy5mZnROb2RlLnNtb290aGluZ1RpbWVDb25zdGFudCA9IDAuOTtcblxuICAgICAgICB0aGlzLm91dEZmdE5vZGUgPSB0aGlzLmF1ZGlvQ29udGV4dC5jcmVhdGVBbmFseXNlcigpO1xuICAgICAgICB0aGlzLm91dEZmdE5vZGUuZmZ0U2l6ZSA9IHRoaXMuZmZ0Tm9kZS5mZnRTaXplO1xuICAgICAgICB0aGlzLm91dEZmdE5vZGUuc21vb3RoaW5nVGltZUNvbnN0YW50ID0gdGhpcy5mZnROb2RlLnNtb290aGluZ1RpbWVDb25zdGFudDtcbiAgICB9XG5cbiAgICBhc3luYyBpbml0SW5wdXQoKSB7XG4gICAgICAgIHRyeSB7XG4gICAgICAgICAgIGNvbnN0IGNvbnN0cmFpbnRzID0ge2F1ZGlvOiB7IG5vaXNlU3VwcHJlc3Npb246IGZhbHNlLCBlY2hvQ2FuY2VsbGF0aW9uOiBmYWxzZSwgYXV0b0dhaW5Db250cm9sOiBmYWxzZSB9fTtcbiAgICAgICAgICAgY29uc3Qgc3RyZWFtID0gYXdhaXQgbmF2aWdhdG9yLm1lZGlhRGV2aWNlcy5nZXRVc2VyTWVkaWEoY29uc3RyYWludHMpO1xuICAgICAgICAgICB0aGlzLmlucHV0ID0gdGhpcy5hdWRpb0NvbnRleHQuY3JlYXRlTWVkaWFTdHJlYW1Tb3VyY2Uoc3RyZWFtKTtcbiAgICAgICAgICAgY29uc29sZS5sb2coXCJbSW5wdXRdIEdvdCBhY2Nlc3MgdG8gaW5wdXQgZGV2aWNlXCIpXG4gICAgICAgICB9IGNhdGNoKGVycikge1xuICAgICAgICAgICAgY29uc29sZS5lcnJvcihcIltJbnB1dF0gQ2FuJ3QgYWNjZXNzIHVzZXIgaW5wdXQgZGV2aWNlXCIpXG4gICAgICAgICB9XG4gICAgfVxuXG4gICAgY29ubmVjdEFsbCgpIHtcbiAgICAgICAgY29uc3Qgc2lnbmFsID0gdGhpcy5pbnB1dFxuICAgICAgICAuY29ubmVjdCh0aGlzLmhpcGFzcylcbiAgICAgICAgLmNvbm5lY3QodGhpcy5sb3Bhc3MpXG5cbiAgICAgICAgc2lnbmFsLmNvbm5lY3QodGhpcy5mZnROb2RlKVxuXG4gICAgICAgIHNpZ25hbC5jb25uZWN0KHRoaXMuY29udi5ub2RlKS5jb25uZWN0KHRoaXMuYXV0b0dhaW4uZ2FpbikuY29ubmVjdCh0aGlzLmF1dG9HYWluLmxpbWl0ZXIpLmNvbm5lY3QodGhpcy5hdWRpb0NvbnRleHQuZGVzdGluYXRpb24pXG4gICAgICAgIHRoaXMuY29udi5jb25uZWN0KHRoaXMub3V0RmZ0Tm9kZSk7XG4gICAgfVxuXG4gICAgdXBkYXRlKCl7XG4gICAgICAgIHJlcXVlc3RBbmltYXRpb25GcmFtZSgoKT0+dGhpcy51cGRhdGUoKSlcbiAgICAgICAgdGhpcy51cGRhdGVTcGVjdHJ1bSgpO1xuICAgICAgICB0aGlzLmF1dG9HYWluLnVwZGF0ZUdhaW4odGhpcy5hbmFsLm1heERiKTtcbiAgICAgICAgdGhpcy5jb252LnVwZGF0ZUtlcm5lbCh0aGlzLmFuYWwubWFnUmVkdWN0aW9uc0FtcCk7XG4gICAgICAgIHRoaXMuZHJhd1NwZWN0cnVtKCk7XG4gICAgfVxuXG4gICAgdXBkYXRlU3BlY3RydW0oKSB7XG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5nZXRGbG9hdEZyZXF1ZW5jeURhdGEodGhpcy5mZnRCdWZmZXIpO1xuICAgICAgICB0aGlzLm91dEZmdE5vZGUuZ2V0RmxvYXRGcmVxdWVuY3lEYXRhKHRoaXMub3V0RmZ0QnVmZmVyKTtcbiAgICAgICAgdGhpcy5tYXhEYiA9IHRoaXMuZmZ0Tm9kZS5tYXhEZWNpYmVscztcbiAgICAgICAgdGhpcy5hbmFsLmFuYWx5c2VTcGVjdHJ1bSh0aGlzLmZmdEJ1ZmZlciwgdGhpcy5tYXhEYik7XG4gICAgfVxuXG4gICAgLy8gZ3JhcGhzXG4gICAgZHJhd1NwZWN0cnVtKCkge1xuICAgICAgICAvL2NvbnN0IHBpdGNoSSA9IHRoaXMuYW5hbC5mYkluZGV4ZXNbMF07XG4gICAgICAgIGNvbnN0IHBpdGNoSSA9IHRoaXMuZmZ0QnVmZmVyLmluZGV4T2YodGhpcy5hbmFsLm1heERiKTtcbiAgICAgICAgY29uc3QgYmdDb2xvciA9IHRoaXMuYW5hbC5uYlBlYWtzID4gMCA/IEdyYXBoLnBpdGNoQW1wVG9IU0xBKHBpdGNoSSAqIHRoaXMuYmluVG9GcmVxLCB0aGlzLmZmdEJ1ZmZlcltwaXRjaEldKSA6ICdyZ2JhKDEwMCwxMDAsMTAwLDAuMyknO1xuICAgICAgICAvL0dyYXBoLmRyYXdCZyhcImNhbnZhcyNzcGVjdHJ1bVwiLCAncmdiYSgwLDAsMCwwLjUpJylcbiAgICAgICAgR3JhcGguZHJhd0JnKFwiY2FudmFzI3NwZWN0cnVtXCIsIGJnQ29sb3IpXG4gICAgICAgIEdyYXBoLmRyYXdHYWluKFwiY2FudmFzI3NwZWN0cnVtXCIsIFwicmVkXCIsIHRoaXMuYXV0b0dhaW4uZ2V0Q3VycmVudFZhbHVlKCksIC0xMDAsIDgwKVxuXG4gICAgICAgIC8vR3JhcGguZHJhd1Ntb290aExvZ0xpbmUoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5mZnRCdWZmZXIsIC0xMDAsIDApO1xuICAgICAgICBHcmFwaC5kcmF3TG9nTGluZShcImNhbnZhcyNzcGVjdHJ1bVwiLCB0aGlzLmZmdEJ1ZmZlciwgLTEwMCwgMCk7XG5cbiAgICAgICAgLy8gR3JhcGguZHJhd1Ntb290aExvZ0xpbmUoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5vdXRGZnRCdWZmZXIsIC0xMDAsIDAsIFwicmdiYSgwLDAsMCwwLjEpXCIpO1xuICAgICAgICBcbiAgICAgICAgLy9HcmFwaC5kcmF3U21vb3RoTG9nTGluZUludmVydGVkKFwiY2FudmFzI3NwZWN0cnVtXCIsIHRoaXMuYW5hbC5tYWdSZWR1Y3Rpb25zLCAtMTAwLCAwLCBcInJnYmEoMCwwLDAsMC41KVwiKTtcbiAgICAgICAgR3JhcGguZHJhd0xvZ0xpbmVJbnZlcnRlZChcImNhbnZhcyNzcGVjdHJ1bVwiLCB0aGlzLmFuYWwubWFnUmVkdWN0aW9ucywgLTEwMCwgMCwgXCJyZ2JhKDAsMCwwLDAuNSlcIik7XG5cbiAgICAgICAgLypcbiAgICAgICAgR3JhcGguZHJhd0hMaW5lKFwiY2FudmFzI3NwZWN0cnVtXCIsICB0aGlzLmFuYWwuYXZlcmFnZURiLCAtMTAwLCAwLCBcImdyZWVuXCIpO1xuICAgICAgICAvL0dyYXBoLmRyYXdITGluZShcImNhbnZhcyNzcGVjdHJ1bVwiLCB0aGlzLmFuYWwubWF4RGIsIC0xMDAsIDAsIFwicmVkXCIpO1xuICAgICAgICBHcmFwaC5kcmF3SExpbmUoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5hbmFsLnBlYWtUaHIsIC0xMDAsIDAsIFwiYmx1ZVwiKTtcbiAgICAgICAgKi9cblxuICAgICAgICBHcmFwaC5kcmF3VmVydGljYWxzKFwiY2FudmFzI3NwZWN0cnVtXCIsIHRoaXMuYW5hbC5uYlBlYWtzLCB0aGlzLmFuYWwucGVha0luZGV4ZXMsIC04MCwgMCwgXCJyZ2JhKDI1NSwyNTUsMjU1LDAuMDEpXCIpXG4gICAgICAgIEdyYXBoLmRyYXdWZXJ0aWNhbHMoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5udW1GYiwgdGhpcy5mYkJ1ZmZlciwgLTgwLCAwLCBcInJnYmEoMjU1LDAsMCwwLjUpXCIpXG5cbiAgICAgICAgLy9HcmFwaC5kcmF3QmcoXCJjYW52YXMjc2xvcGVcIiwgXCJ3aGl0ZVwiKVxuICAgICAgICAvL0dyYXBoLmRyYXdMaW5lKFwiY2FudmFzI3Nsb3BlXCIsIHRoaXMuY29udi5pckJ1ZmZlci5nZXRDaGFubmVsRGF0YSgwKSwgLTEsIDEsIFwicmdiYSgyNTUsMCwwLDEpXCIpO1xuICAgIH1cbn1cblxubW9kdWxlLmV4cG9ydHMgPSBTcXVpZGJhY2tGRlRQcm9jZXNzIiwiY29uc3QgU3BsaW5lID0gcmVxdWlyZShcImN1YmljLXNwbGluZVwiKVxuXG5jbGFzcyBTcXVpZGJhY2tHcmFwaCB7XG5cbiAgICBkcmF3Qmcoc2VsZWN0b3IsIGNvbG9yKSB7XG4gICAgICAgIGNvbnN0IGNhbnZhcyA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3Ioc2VsZWN0b3IpO1xuICAgICAgICBjb25zdCBjYW52YXNDdHggPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuICAgICAgICBjYW52YXNDdHguYmVnaW5QYXRoKCk7XG4gICAgICAgIGNhbnZhc0N0eC5maWxsU3R5bGUgPSBjb2xvcjtcbiAgICAgICAgY2FudmFzQ3R4LmZpbGxSZWN0KDAsIDAsIGNhbnZhcy53aWR0aCwgY2FudmFzLmhlaWdodCk7XG4gICAgICAgIGNhbnZhc0N0eC5jbG9zZVBhdGgoKTtcbiAgICB9XG5cbiAgICBkcmF3R2FpbihzZWxlY3RvciwgY29sb3IsIGFtcD0wLCBtaW49MCwgbWF4PTEpIHtcbiAgICAgICAgY29uc3QgY2FudmFzID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihzZWxlY3Rvcik7XG4gICAgICAgIGNvbnN0IGNhbnZhc0N0eCA9IGNhbnZhcy5nZXRDb250ZXh0KFwiMmRcIik7XG4gICAgICAgIGNvbnN0IGdhaW5IZWlnaHQgPSAoMjAgKiBNYXRoLmxvZzEwKGFtcCkgLSBtaW4pIC8gKG1heC1taW4pO1xuICAgICAgICAvL2NvbnNvbGUubG9nKGFtcCwgZ2FpbkhlaWdodCwgbWluLCBtYXgpXG4gICAgICAgIGNhbnZhc0N0eC5iZWdpblBhdGgoKTtcbiAgICAgICAgY2FudmFzQ3R4LmZpbGxTdHlsZSA9IGNvbG9yO1xuICAgICAgICBjYW52YXNDdHguZmlsbFJlY3QoY2FudmFzLndpZHRoLTIsIGNhbnZhcy5oZWlnaHQsIGNhbnZhcy53aWR0aCwgLWNhbnZhcy5oZWlnaHQgKiBnYWluSGVpZ2h0KTtcbiAgICAgICAgY2FudmFzQ3R4LmNsb3NlUGF0aCgpO1xuICAgIH1cblxuICAgIHBpdGNoQW1wVG9IU0xBKHBpdGNoLGFtcD0xLGFscGhhPTAuMSl7XG4gICAgICAgIGxldCBvY3RhdmVzID0gTWF0aC5sb2cyKHBpdGNoLzQ0MCk7XG4gICAgICAgIGxldCBwYyA9IDEyICogb2N0YXZlcyAlIDEyO1xuICAgICAgICBsZXQgaCA9IHBjLzEyKjM2MDtcbiAgICAgICAgbGV0IGwgPSAxMCAqIChvY3RhdmVzICsgNS41KTtcbiAgICAgICAgbGV0IHMgPSAyNTUgKyBhbXA7XG4gICAgICAgIGxldCBjb2RlID0gYGhzbGEoJHtofSwke3M8MD8wOnM+MTAwPzEwMDpzfSUsJHtsPDU/NTpsPjEwMD8xMDA6bH0lLCR7YWxwaGF9KWA7XG4gICAgICAgIC8vY29uc29sZS5sb2cocGl0Y2gscGMsb2N0YXZlcyxzLGgpO1xuICAgICAgICByZXR1cm4gY29kZTtcbiAgICB9XG5cbiAgICBkcmF3TGluZShzZWxlY3RvciwgZGF0YSwgbWluLCBtYXgsIGNvbG9yPVwicmdiYSgyNTUsMjU1LDI1NSwxKVwiKSB7XG4gICAgICAgIGNvbnN0IGNhbnZhcyA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3Ioc2VsZWN0b3IpO1xuICAgICAgICBjb25zdCBjYW52YXNDdHggPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuICAgICAgICBjb25zdCBiYXJXaWR0aCA9IChjYW52YXMud2lkdGggLyBkYXRhLmxlbmd0aCk7XG4gICAgICAgIGxldCBiYXJIZWlnaHQ7XG4gICAgICAgIGxldCB4ID0gMDtcbiAgICAgICAgY2FudmFzQ3R4LmJlZ2luUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHgubW92ZVRvKHgsIGNhbnZhcy5oZWlnaHQgLSAoZGF0YVswXSAtIChtaW4pICkgLyAobWF4LW1pbikgKiBjYW52YXMuaGVpZ2h0KTtcblxuICAgICAgICBmb3IobGV0IGkgPSAxOyBpIDwgZGF0YS5sZW5ndGg7ICsraSl7XG4gICAgICAgICAgICBiYXJIZWlnaHQgPSAoZGF0YVtpXSAtIChtaW4pICkgLyAobWF4LW1pbikgKiBjYW52YXMuaGVpZ2h0O1xuICAgICAgICAgICAgY2FudmFzQ3R4LmxpbmVUbyh4ICsgYmFyV2lkdGgsIGNhbnZhcy5oZWlnaHQtYmFySGVpZ2h0KTtcbiAgICAgICAgICAgIHggKz0gYmFyV2lkdGg7XG4gICAgICAgIH1cbiAgICAgICAgY2FudmFzQ3R4LnN0cm9rZVN0eWxlID0gY29sb3I7XG4gICAgICAgIGNhbnZhc0N0eC5zdHJva2VXaWR0aCA9IDAuMTtcbiAgICAgICAgY2FudmFzQ3R4LnN0cm9rZSgpO1xuICAgIH1cblxuICAgIGRyYXdMb2dMaW5lKHNlbGVjdG9yLCBkYXRhLCBtaW4sIG1heCwgY29sb3I9XCJyZ2JhKDI1NSwyNTUsMjU1LDAuNSlcIikge1xuICAgICAgICBjb25zdCBjYW52YXMgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKHNlbGVjdG9yKTtcbiAgICAgICAgY29uc3QgY2FudmFzQ3R4ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcbiAgICAgICAgY29uc3QgYmFyV2lkdGggPSAoY2FudmFzLndpZHRoIC8gKGRhdGEubGVuZ3RoIC0gMSkpO1xuICAgICAgICBsZXQgYmFySGVpZ2h0O1xuICAgICAgICBjb25zdCBsb2dJVG9YID0gY2FudmFzLndpZHRoIC8gTWF0aC5sb2coZGF0YS5sZW5ndGgpXG4gICAgICAgIGNhbnZhc0N0eC5iZWdpblBhdGgoKTtcbiAgICAgICAgY2FudmFzQ3R4Lm1vdmVUbygwLCBjYW52YXMuaGVpZ2h0KTtcbiAgICAgICAgY2FudmFzQ3R4LmxpbmVUbygwLCBjYW52YXMuaGVpZ2h0IC0gKGRhdGFbMF0gLSAobWluKSApIC8gKG1heC1taW4pICogY2FudmFzLmhlaWdodCk7XG5cbiAgICAgICAgZm9yKGxldCBpID0gMTsgaSA8IGRhdGEubGVuZ3RoOyArK2kpe1xuICAgICAgICAgICAgYmFySGVpZ2h0ID0gKGRhdGFbaV0gLSBtaW4pIC8gKG1heC1taW4pICogY2FudmFzLmhlaWdodDtcbiAgICAgICAgICAgIGNvbnN0IHggPSAoTWF0aC5sb2coaSkgKyBNYXRoLmxvZyhpKzEpKSAqIGxvZ0lUb1ggLyAyO1xuICAgICAgICAgICAgY2FudmFzQ3R4LmxpbmVUbyh4LCBjYW52YXMuaGVpZ2h0LWJhckhlaWdodCk7XG4gICAgICAgIH1cbiAgICAgICAgY2FudmFzQ3R4LmxpbmVUbyhjYW52YXMud2lkdGgsIGNhbnZhcy5oZWlnaHQpXG4gICAgICAgIGNhbnZhc0N0eC5jbG9zZVBhdGgoKTtcbiAgICAgICAgY2FudmFzQ3R4LmZpbGxTdHlsZSA9IGNvbG9yXG4gICAgICAgIGNhbnZhc0N0eC5maWxsKCk7XG4gICAgfVxuXG4gICAgZHJhd1Ntb290aExvZ0xpbmUoc2VsZWN0b3IsIGRhdGEsIG1pbiwgbWF4LCBjb2xvcj1cInJnYmEoMjU1LDI1NSwyNTUsMC41KVwiLCB1cHNhbXBsZT0zLCBzdHJva2UgPSBmYWxzZSkge1xuICAgICAgICBjb25zdCBjYW52YXMgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKHNlbGVjdG9yKTtcbiAgICAgICAgY29uc3QgY2FudmFzQ3R4ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcbiAgICAgICAgY29uc3QgYmFyV2lkdGggPSAoY2FudmFzLndpZHRoIC8gKGRhdGEubGVuZ3RoIC0gMSkpO1xuICAgICAgICBsZXQgYmFySGVpZ2h0O1xuICAgICAgICBjb25zdCBsb2dJVG9YID0gY2FudmFzLndpZHRoIC8gTWF0aC5sb2coZGF0YS5sZW5ndGgpO1xuICAgICAgICBjb25zdCB4cyA9IFsuLi5kYXRhLmtleXMoKV07XG4gICAgICAgIGRhdGEgPSBuZXcgU3BsaW5lKHhzLCBkYXRhKTtcbiAgICAgICAgY2FudmFzQ3R4LmJlZ2luUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHgubW92ZVRvKDAsIGNhbnZhcy5oZWlnaHQpO1xuICAgICAgICBjb25zdCBpbmMgPSAxL3Vwc2FtcGxlO1xuICAgICAgICBmb3IobGV0IGkgPSAwOyBpIDwgKGRhdGEueXMubGVuZ3RoKTsgaSs9aW5jKXtcbiAgICAgICAgICAgIGNvbnN0IHZhbCA9IGRhdGEuYXQoaSk7IC8vIGRhdGFbaV1cbi8vY29uc29sZS5sb2coaSlcbiAgICAgICAgICAgIGJhckhlaWdodCA9ICggdmFsIC0gbWluKSAvIChtYXgtbWluKSAqIGNhbnZhcy5oZWlnaHQ7XG4gICAgICAgICAgICBjb25zdCB4ID0gKE1hdGgubG9nKGkpICsgTWF0aC5sb2coaSsxKSkgKiBsb2dJVG9YIC8gMjtcbiAgICAgICAgICAgIGNhbnZhc0N0eC5saW5lVG8oeCwgY2FudmFzLmhlaWdodC1iYXJIZWlnaHQpO1xuICAgICAgICB9XG4gICAgICAgIGNhbnZhc0N0eC5saW5lVG8oY2FudmFzLndpZHRoLCBjYW52YXMuaGVpZ2h0KVxuICAgICAgICBjYW52YXNDdHguY2xvc2VQYXRoKCk7XG4gICAgICAgIGlmKHN0cm9rZSkge1xuICAgICAgICAgICAgY2FudmFzQ3R4LnN0cm9rZVN0eWxlID0gY29sb3JcbiAgICAgICAgICAgIGNhbnZhc0N0eC5zdHJva2UoKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGNhbnZhc0N0eC5maWxsU3R5bGUgPSBjb2xvclxuICAgICAgICAgICAgY2FudmFzQ3R4LmZpbGwoKTtcbiAgICAgICAgfVxuXG4gICAgfVxuXG4gICAgZHJhd0xvZ0xpbmVJbnZlcnRlZChzZWxlY3RvciwgZGF0YSwgbWluLCBtYXgsIGNvbG9yPVwicmdiYSgyNTUsMjU1LDI1NSwwLjUpXCIpIHtcbiAgICAgICAgY29uc3QgY2FudmFzID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihzZWxlY3Rvcik7XG4gICAgICAgIGNvbnN0IGNhbnZhc0N0eCA9IGNhbnZhcy5nZXRDb250ZXh0KFwiMmRcIik7XG4gICAgICAgIGNvbnN0IGJhcldpZHRoID0gKGNhbnZhcy53aWR0aCAvIChkYXRhLmxlbmd0aCAtIDEpKTtcbiAgICAgICAgbGV0IGJhckhlaWdodDtcbiAgICAgICAgY29uc3QgbG9nSVRvWCA9IGNhbnZhcy53aWR0aCAvIE1hdGgubG9nKGRhdGEubGVuZ3RoIC0gMSlcbiAgICAgICAgY2FudmFzQ3R4LmJlZ2luUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHgubW92ZVRvKDAsIDApO1xuICAgICAgICBjYW52YXNDdHgubGluZVRvKDAsIGNhbnZhcy5oZWlnaHQtKGRhdGFbMF0gLSAobWluKSApIC8gKG1heC1taW4pICogY2FudmFzLmhlaWdodCk7XG5cbiAgICAgICAgZm9yKGxldCBpID0gMTsgaSA8IGRhdGEubGVuZ3RoOyArK2kpe1xuICAgICAgICAgICAgYmFySGVpZ2h0ID0gKGRhdGFbaV0gLSBtaW4pIC8gKG1heC1taW4pICogY2FudmFzLmhlaWdodDtcbiAgICAgICAgICAgIGNvbnN0IHggPSAoTWF0aC5sb2coaSkgKyBNYXRoLmxvZyhpKzEpKSAqIGxvZ0lUb1ggLyAyO1xuICAgICAgICAgICAgY2FudmFzQ3R4LmxpbmVUbyh4LCBjYW52YXMuaGVpZ2h0LWJhckhlaWdodCk7XG4gICAgICAgIH1cbiAgICAgICAgY2FudmFzQ3R4LmxpbmVUbyhjYW52YXMud2lkdGgsIDApXG4gICAgICAgIGNhbnZhc0N0eC5jbG9zZVBhdGgoKTtcbiAgICAgICAgY2FudmFzQ3R4LmZpbGxTdHlsZSA9IGNvbG9yXG4gICAgICAgIGNhbnZhc0N0eC5maWxsKCk7XG4gICAgfVxuXG4gICAgZHJhd1Ntb290aExvZ0xpbmVJbnZlcnRlZChzZWxlY3RvciwgZGF0YSwgbWluLCBtYXgsIGNvbG9yPVwicmdiYSgyNTUsMjU1LDI1NSwwLjUpXCIsIHVwc2FtcGxlID0gMykge1xuICAgICAgICBjb25zdCBjYW52YXMgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKHNlbGVjdG9yKTtcbiAgICAgICAgY29uc3QgY2FudmFzQ3R4ID0gY2FudmFzLmdldENvbnRleHQoXCIyZFwiKTtcbiAgICAgICAgY29uc3QgYmFyV2lkdGggPSAoY2FudmFzLndpZHRoIC8gKGRhdGEubGVuZ3RoIC0gMSkpO1xuICAgICAgICBsZXQgYmFySGVpZ2h0O1xuICAgICAgICBjb25zdCBsb2dJVG9YID0gY2FudmFzLndpZHRoIC8gTWF0aC5sb2coZGF0YS5sZW5ndGggLSAxKVxuICAgICAgICBjb25zdCB4cyA9IFsuLi5kYXRhLmtleXMoKV07XG4gICAgICAgIGRhdGEgPSBuZXcgU3BsaW5lKHhzLCBkYXRhKTtcbiAgICAgICAgY2FudmFzQ3R4LmJlZ2luUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHgubW92ZVRvKDAuNSwgMCk7XG4gICAgICAgIGNvbnN0IGluYyA9IDEvdXBzYW1wbGU7XG4gICAgICAgIGZvcihsZXQgaSA9IDA7IGkgPCAoZGF0YS55cy5sZW5ndGgpOyBpKz1pbmMpe1xuICAgICAgICAgICAgY29uc3QgdmFsID0gZGF0YS5hdChpKTsgLy8gZGF0YVtpXVxuICAgICAgICAgICAgYmFySGVpZ2h0ID0gKHZhbCAtIG1pbikgLyAobWF4LW1pbikgKiBjYW52YXMuaGVpZ2h0O1xuICAgICAgICAgICAgY29uc3QgeCA9IChNYXRoLmxvZyhpKSArIE1hdGgubG9nKGkrMSkpICogbG9nSVRvWCAvIDI7XG4gICAgICAgICAgICBjYW52YXNDdHgubGluZVRvKHgsIGNhbnZhcy5oZWlnaHQtYmFySGVpZ2h0KTtcbiAgICAgICAgfVxuICAgICAgICBjYW52YXNDdHgubGluZVRvKGNhbnZhcy53aWR0aCwgMClcbiAgICAgICAgY2FudmFzQ3R4LmNsb3NlUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHguZmlsbFN0eWxlID0gY29sb3JcbiAgICAgICAgY2FudmFzQ3R4LmZpbGwoKTtcbiAgICB9XG5cbiAgICBkcmF3VmVydGljYWxzKHNlbGVjdG9yLCBudW1Qb2ludHMsIGRhdGEsIG1pbiwgbWF4LCBjb2xvcj1cInJnYmEoMjU1LDI1NSwyNTUsMC41KVwiKSB7XG4gICAgICAgIGNvbnN0IGNhbnZhcyA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3Ioc2VsZWN0b3IpO1xuICAgICAgICBjb25zdCBjYW52YXNDdHggPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuICAgICAgICBjb25zdCBiYXJXaWR0aCA9IChjYW52YXMud2lkdGggLyAoZGF0YS5sZW5ndGggLSAxKSk7XG4gICAgICAgIGNvbnN0IGxvZ0lUb1ggPSBjYW52YXMud2lkdGggLyBNYXRoLmxvZyhkYXRhLmxlbmd0aCAtIDEpO1xuICAgICAgICBjYW52YXNDdHguZmlsbFN0eWxlID0gY29sb3I7XG4gICAgICAgIGNhbnZhc0N0eC5iZWdpblBhdGgoKTtcbiAgICAgICAgXG4gICAgICAgIGZvcihsZXQgbiA9IDA7IG4gPCBudW1Qb2ludHM7ICsrbil7XG4gICAgICAgICAgICBjb25zdCBpID0gZGF0YVtuXTtcbiAgICAgICAgICAgIC8vY29uc29sZS5sb2coaSk7XG4gICAgICAgICAgICBjb25zdCB4ID0gTWF0aC5sb2coaSkgKiBsb2dJVG9YO1xuICAgICAgICAgICAgY2FudmFzQ3R4LmZpbGxSZWN0KHgsIDAsIE1hdGgubG9nKGkrMSkgKiBsb2dJVG9YIC0geCwgY2FudmFzLmhlaWdodClcbiAgICAgICAgICAgIC8vY2FudmFzQ3R4Lm1vdmVUbyh4ICsgYmFyV2lkdGgqMC41LCBjYW52YXMuaGVpZ2h0KVxuICAgICAgICAgICAgLy9jYW52YXNDdHgubGluZVRvKHggKyBiYXJXaWR0aCowLjUsMCk7XG4gICAgICAgIH1cbiAgICAgICAgLy9jYW52YXNDdHguc3Ryb2tlU3R5bGUgPSBjb2xvclxuICAgICAgICAvL2NhbnZhc0N0eC5saW5lV2lkdGggPSAxO1xuICAgICAgICAvL2NhbnZhc0N0eC5zdHJva2UoKTtcbiAgICAgICAgY2FudmFzQ3R4LmZpbGwoKTtcbiAgICB9XG5cbiAgICBkcmF3SExpbmUoc2VsZWN0b3IsIGgsIG1pbiwgbWF4LCBjb2xvcj1cInJnYmEoMjU1LDI1NSwyNTUsMC41KVwiKSB7XG4gICAgICAgIGNvbnN0IGNhbnZhcyA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3Ioc2VsZWN0b3IpO1xuICAgICAgICBjb25zdCBjYW52YXNDdHggPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuICAgICAgICBjb25zdCB5ID0gY2FudmFzLmhlaWdodCAqICgxIC0gKGggLSBtaW4pIC8gKG1heC1taW4pKTtcbiAgICAgICAgY2FudmFzQ3R4LmJlZ2luUGF0aCgpO1xuICAgICAgICBjYW52YXNDdHgubW92ZVRvKDAsIHkpO1xuICAgICAgICBjYW52YXNDdHgubGluZVRvKGNhbnZhcy53aWR0aCwgeSk7XG4gICAgICAgIGNhbnZhc0N0eC5zdHJva2VTdHlsZSA9IGNvbG9yXG4gICAgICAgIGNhbnZhc0N0eC5saW5lV2lkdGggPSAwLjU7XG4gICAgICAgIGNhbnZhc0N0eC5zdHJva2UoKTtcbiAgICB9XG5cbn1cblxubW9kdWxlLmV4cG9ydHMgPSBTcXVpZGJhY2tHcmFwaCIsIlwidXNlIHN0cmljdFwiO1xuXG5jb25zdCBTcXVpZGJhY2tXb3JrbGV0UHJvY2VzcyA9IHJlcXVpcmUoJy4vd29ya2xldFByb2Nlc3MuanMnKVxuY29uc3QgU3F1aWRiYWNrRkZUUHJvY2VzcyA9IHJlcXVpcmUoJy4vZmZ0UHJvY2Vzcy5qcycpXG5cbmFzeW5jIGZ1bmN0aW9uIGluaXQoKSB7XG4gICAgY29uc3QgYXVkaW9Db250ZXh0ID0gbmV3IEF1ZGlvQ29udGV4dCgpXG4gICAgaWYgKGF1ZGlvQ29udGV4dC5hdWRpb1dvcmtsZXQgPT09IHVuZGVmaW5lZCkge1xuICAgICAgICBoYW5kbGVOb1dvcmtsZXQoKTtcbiAgICAgICAgcmV0dXJuO1xuICAgIH1cbiAgICBsZXQgcHJvY2VzcyA9IG5ldyBTcXVpZGJhY2tGRlRQcm9jZXNzKGF1ZGlvQ29udGV4dCwgNTEyKTtcbiAgICBhd2FpdCBwcm9jZXNzLnN0YXJ0KClcbn1cblxuZnVuY3Rpb24gaGFuZGxlTm9Xb3JrbGV0KCkge1xuICAgIGxldCAkbm9Xb3JrbGV0ID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihcIiNuby13b3JrbGV0XCIpO1xuICAgICRub1dvcmtsZXQuc3R5bGUuZGlzcGxheSA9ICdibG9jayc7XG4gICAgbGV0ICR0aW1lbGluZSA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoXCIudGltZWxpbmVcIik7XG4gICAgJHRpbWVsaW5lLnN0eWxlLmRpc3BsYXkgPSAnbm9uZSc7XG4gICAgbGV0ICRjb250cm9scyA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoXCIuY29udHJvbHNcIik7XG4gICAgJGNvbnRyb2xzLnN0eWxlLmRpc3BsYXkgPSAnbm9uZSc7XG59XG5cbmZ1bmN0aW9uIHJlc2l6ZShjYW52YXMpIHtcbiAgLy8gTG9va3VwIHRoZSBzaXplIHRoZSBicm93c2VyIGlzIGRpc3BsYXlpbmcgdGhlIGNhbnZhcy5cbiAgdmFyIGRpc3BsYXlXaWR0aCAgPSBjYW52YXMuY2xpZW50V2lkdGg7XG4gIHZhciBkaXNwbGF5SGVpZ2h0ID0gY2FudmFzLmNsaWVudEhlaWdodDtcbiBcbiAgLy8gQ2hlY2sgaWYgdGhlIGNhbnZhcyBpcyBub3QgdGhlIHNhbWUgc2l6ZS5cbiAgaWYgKGNhbnZhcy53aWR0aCAgIT0gZGlzcGxheVdpZHRoIHx8XG4gICAgICBjYW52YXMuaGVpZ2h0ICE9IGRpc3BsYXlIZWlnaHQpIHtcbiBcbiAgICAvLyBNYWtlIHRoZSBjYW52YXMgdGhlIHNhbWUgc2l6ZVxuICAgIGNhbnZhcy53aWR0aCAgPSBkaXNwbGF5V2lkdGg7XG4gICAgY2FudmFzLmhlaWdodCA9IGRpc3BsYXlIZWlnaHQ7XG4gIH1cbn1cblxuZnVuY3Rpb24gcmVzaXplQWxsQ2FudmFzKCkge1xuICAgIGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3JBbGwoXCJjYW52YXNcIikuZm9yRWFjaChjYW52YXM9PnJlc2l6ZShjYW52YXMpKVxufVxuXG5cbndpbmRvdy5hZGRFdmVudExpc3RlbmVyKCdsb2FkJywgKCk9PntcbiAgICBjb25zdCBidXR0b24gPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKFwiYnV0dG9uI3N0YXJ0XCIpXG4gICAgYnV0dG9uLmFkZEV2ZW50TGlzdGVuZXIoXCJjbGlja1wiLCAoKT0+e1xuICAgICAgICBpbml0KCk7XG4gICAgICAgIGJ1dHRvbi5jbGFzc0xpc3QuYWRkKFwic3RhcnRlZFwiKVxuICAgIH0pO1xuICAgIHJlc2l6ZUFsbENhbnZhcygpO1xufSk7XG5cbndpbmRvdy5hZGRFdmVudExpc3RlbmVyKCdyZXNpemUnLCAoKT0+e1xuICAgIHJlc2l6ZUFsbENhbnZhcygpXG59KTtcbiIsImNvbnN0IHBhcnRjaFNjYWxlNyA9IFszMCwzMC4zNzUsMzEuNSwzMiwzMy4zMzMzMywzMy43NSwzNC4yODU3MSwzNSwzNS41NTU1NiwzNiwzNy41LDM4LjU3MTQzLDM5LjM3NSw0MCw0MC41LDQyLDQyLjg1NzE0LDQ0LjQ0NDQ0LDQ1LDQ2LjY2NjY3LDQ3LjE0Mjg2LDQ4LDUwLDUwLjYyNSw1MS40Mjg1Nyw1Mi41LDUzLjMzMzMzLDU0LDU2LjI1LDU3LjE0Mjg2LDU5LjI1OTI2LDYwLDYwLjc1LDYzLDY0LDY2LjY2NjY3LDY3LjUsNjguNTcxNDMsNzAsNzEuMTExMTEsNzIsNzUsNzcuMTQyODYsNzguNzUsODAsODEsODQsODUuNzE0MjksODguODg4ODksOTAsOTMuMzMzMzMsOTQuMjg1NzEsOTYsMTAwLDEwMS4yNSwxMDIuODU3MSwxMDUsMTA2LjY2NjcsMTA4LDExMi41LDExNC4yODU3LDExOC41MTg1LDEyMCwxMjEuNSwxMjYsMTI4LDEzMy4zMzMzLDEzNSwxMzcuMTQyOSwxNDAsMTQyLjIyMjIsMTQ0LDE1MCwxNTQuMjg1NywxNTcuNSwxNjAsMTYyLDE2OCwxNzEuNDI4NiwxNzcuNzc3OCwxODAsMTg2LjY2NjcsMTg4LjU3MTQsMTkyLDIwMCwyMDIuNSwyMDUuNzE0MywyMTAsMjEzLjMzMzMsMjE2LDIyNSwyMjguNTcxNCwyMzcuMDM3LDI0MCwyNDMsMjUyLDI1NiwyNjYuNjY2NywyNzAsMjc0LjI4NTcsMjgwLDI4NC40NDQ0LDI4OCwzMDAsMzA4LjU3MTQsMzE1LDMyMCwzMjQsMzM2LDM0Mi44NTcxLDM1NS41NTU2LDM2MCwzNzMuMzMzMywzNzcuMTQyOSwzODQsNDAwLDQwNSw0MTEuNDI4Niw0MjAsNDI2LjY2NjcsNDMyLDQ1MCw0NTcuMTQyOSw0NzQuMDc0MSw0ODAsNDg2LDUwNCw1MTIsNTMzLjMzMzMsNTQwLDU0OC41NzE0LDU2MCw1NjguODg4OSw1NzYsNjAwLDYxNy4xNDI5LDYzMCw2NDAsNjQ4LDY3Miw2ODUuNzE0Myw3MTEuMTExMSw3MjAsNzQ2LjY2NjcsNzU0LjI4NTcsNzY4LDgwMCw4MTAsODIyLjg1NzEsODQwLDg1My4zMzMzLDg2NCw5MDAsOTE0LjI4NTcsOTQ4LjE0ODEsOTYwLDk3MiwxMDA4LDEwMjQsMTA2Ni42NjcsMTA4MCwxMDk3LjE0MywxMTIwLDExMzcuNzc4LDExNTIsMTIwMCwxMjM0LjI4NiwxMjYwLDEyODAsMTI5NiwxMzQ0LDEzNzEuNDI5LDE0MjIuMjIyLDE0NDAsMTQ5My4zMzMsMTUwOC41NzEsMTUzNiwxNjAwLDE2MjAsMTY0NS43MTQsMTY4MCwxNzA2LjY2NywxNzI4LDE4MDAsMTgyOC41NzEsMTg5Ni4yOTYsMTkyMCwxOTQ0LDIwMTYsMjA0OCwyMTMzLjMzMywyMTYwLDIxOTQuMjg2LDIyNDAsMjI3NS41NTYsMjMwNCwyNDAwLDI0NjguNTcxLDI1MjAsMjU2MCwyNTkyLDI2ODgsMjc0Mi44NTcsMjg0NC40NDQsMjg4MCwyOTg2LjY2NywzMDE3LjE0MywzMDcyLDMyMDAsMzI0MCwzMjkxLjQyOSwzMzYwLDM0MTMuMzMzLDM0NTYsMzYwMCwzNjU3LjE0MywzNzkyLjU5MywzODQwLDM4ODgsNDAzMiw0MDk2LDQyNjYuNjY3LDQzMjAsNDM4OC41NzEsNDQ4MCw0NTUxLjExMSw0NjA4LDQ4MDAsNDkzNy4xNDMsNTA0MCw1MTIwLDUxODQsNTM3Niw1NDg1LjcxNCw1Njg4Ljg4OSw1NzYwLDU5NzMuMzMzLDYwMzQuMjg2LDYxNDQsNjQwMCw2NDgwLDY1ODIuODU3LDY3MjAsNjgyNi42NjcsNjkxMiw3MjAwLDczMTQuMjg2LDc1ODUuMTg1LDc2ODAsNzc3Niw4MDY0LDgxOTIsODUzMy4zMzMsODY0MCw4Nzc3LjE0Myw4OTYwLDkxMDIuMjIyLDkyMTYsOTYwMCw5ODc0LjI4NiwxMDA4MCwxMDI0MCwxMDM2OCwxMDc1MiwxMDk3MS40MywxMTM3Ny43OCwxMTUyMCwxMTk0Ni42NywxMjA2OC41NywxMjI4OCwxMjgwMCwxMjk2MCwxMzE2NS43MSwxMzQ0MCwxMzY1My4zMywxMzgyNCwxNDQwMCwxNDYyOC41NywxNTE3MC4zNywxNTM2MCwxNTU1MiwxNjEyOCwxNjM4NCwxNzA2Ni42NywxNzI4MCwxNzU1NC4yOSwxNzkyMCwxODIwNC40NF07XG5cbmZ1bmN0aW9uIGNocm9tYXRpYyhtaW49MzAsbWF4PTIwMDAwLGRpdj0xMil7XG4gIHZhciBuX29jdGF2ZXMgPSBNYXRoLmxvZzIobWF4L21pbik7XG4gIHZhciBuX2ZyZXFzID0gTWF0aC5jZWlsKG5fb2N0YXZlcyAqIGRpdikrMTtcbiAgdmFyIHNjYWxlID0gW21pbl07XG4gIGZvcihsZXQgaT0xOyBpPG5fZnJlcXM7IGkrKyl7XG4gICAgc2NhbGUgPSBbLi4uc2NhbGUsIHNjYWxlW3NjYWxlLmxlbmd0aC0xXSAqIE1hdGgucG93KDIsIDEvZGl2KSBdO1xuICB9O1xuICByZXR1cm4gc2NhbGU7XG59XG5cbmZ1bmN0aW9uIGNocm9tYXRpY0ZpbHRlcihtaW49MzAsbWF4PTIwMDAwLGRpdj0xMil7XG4gICAgY29uc3Qgc2NhbGUgPSBjaHJvbWF0aWMobWluLCBtYXgsIGRpdik7XG4gICAgY29uc3QgcXMgPSBjYWxjdWxhdGVRcyhzY2FsZSk7XG4gICAgY29uc3QgYXZncyA9IHNjYWxlLm1hcCgoZnJlcSxpKT0+e1xuICAgICAgICBpZihpID09IDAgfHwgaSA9PSBzY2FsZS5sZW5ndGggLSAxKSByZXR1cm4gZnJlcVxuICAgICAgICByZXR1cm4gc2NhbGVbaS0xXSArIHNjYWxlW2krMV0gLyAyXG4gICAgfSlcblxuICAgIHJldHVybiB7c2NhbGUsIHFzfVxufVxuXG5mdW5jdGlvbiBiaW5zKG51bUJpbnMsIHNyPTQ0MTAwKSB7XG4gICAgY29uc3QgYmluVG9GcmVxID0gc3IgLyAyIC8gbnVtQmlucztcbiAgICB2YXIgc2NhbGUgPSBuZXcgQXJyYXkobnVtQmlucyk7XG4gICAgdmFyIHFzID0gbmV3IEFycmF5KG51bUJpbnMpO1xuICAgIGZvcihsZXQgaSA9IDA7IGkgPCBudW1CaW5zOyArK2kpIHtcbiAgICAgICAgc2NhbGVbaV0gPSBiaW5Ub0ZyZXEgKiBpO1xuICAgICAgICBxc1tpXSA9IGk7XG4gICAgfVxuICAgIHJldHVybiB7c2NhbGUsIHFzfVxuXG59XG5cbmZ1bmN0aW9uIGNhbGN1bGF0ZVFzKHNjYWxlLCBvdmVybGFwPTApe1xuICAgIHJldHVybiBzY2FsZS5tYXAoKGZyZXEsZnJlcV9pKT0+e1xuICAgICAgICBsZXQgZjAsIGYxLCBkZWx0YV9mO1xuICAgICAgICBpZihmcmVxX2k9PTApe1xuICAgICAgICAgICAgZjAgPSBmcmVxOyBmMSA9IHNjYWxlW2ZyZXFfaSsxXTtcbiAgICAgICAgfSBlbHNlIGlmKGZyZXFfaSA9PSAoc2NhbGUubGVuZ3RoLTEpKSB7XG4gICAgICAgICAgICBmMSA9IGZyZXE7IGYwID0gc2NhbGVbZnJlcV9pLTFdO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgZjAgPSBzY2FsZVtmcmVxX2ktMV07ICBmMSA9IHNjYWxlW2ZyZXFfaSsxXTtcbiAgICAgICAgICAgIC8qaWYob3ZlcmxhcCA+IDApe1xuICAgICAgICAgICAgICBkZWx0YV9mID0gTWF0aC5tYXgoZGVsdGFfZjAsZGVsdGFfZjEpICogMiAqIG92ZXJsYXA7XG4gICAgICAgICAgICB9ZWxzZXtcbiAgICAgICAgICAgICAgZGVsdGFfZiA9IE1hdGgubWluKGRlbHRhX2YwLGRlbHRhX2YxKSAqIDI7XG4gICAgICAgICAgICB9Ki9cbiAgICAgICAgfVxuICAgICAgICBkZWx0YV9mID0gZjEgLSBmMDtcbiAgICAgICAgZnJlcSA9ICAoZjErZjApIC8gMjtcbiAgICAgICAgcmV0dXJuIGZyZXEvZGVsdGFfZlxuICAgIH0pO1xufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHsgcGFydGNoU2NhbGU3LCBjaHJvbWF0aWMsIGNocm9tYXRpY0ZpbHRlciwgY2FsY3VsYXRlUXMsIGJpbnMgfTtcbiIsImNvbnN0IFNxdWlkYmFja0dyYXBoID0gcmVxdWlyZSgnLi9ncmFwaC5qcycpXG5jb25zdCB7IGNocm9tYXRpY0ZpbHRlciB9ID0gcmVxdWlyZSgnLi9zY2FsZXMuanMnKVxuY29uc3QgR3JhcGggPSBuZXcgU3F1aWRiYWNrR3JhcGgoKTtcblxuY2xhc3MgU3F1aWRiYWNrV29ya2xldFByb2Nlc3Mge1xuXG4gICAgY29uc3RydWN0b3IoYXVkaW9Db250ZXh0LCBmZnRTaXplID0gMTAyNCkge1xuICAgICAgICB0aGlzLmF1ZGlvQ29udGV4dCA9IGF1ZGlvQ29udGV4dDtcbiAgICAgICAgdGhpcy5mZnRTaXplID0gZmZ0U2l6ZTtcbiAgICAgICAgdGhpcy5nYWluSW5jcmVtZW50ID0gMC4wMDE7XG4gICAgICAgIHRoaXMuZ2FpbkRlY3JlbWVudCA9IDAuMDAxO1xuICAgICAgICB0aGlzLm1pblZvbHVtZSA9IC0xMDtcbiAgICAgICAgdGhpcy5tYXhWb2x1bWUgPSAtMTA7XG5cbiAgICAgICAgdGhpcy5pbml0QnVmZmVycygpO1xuICAgIH1cblxuICAgIGFzeW5jIHN0YXJ0KCkge1xuICAgICAgICBhd2FpdCB0aGlzLmF1ZGlvQ29udGV4dC5hdWRpb1dvcmtsZXQuYWRkTW9kdWxlKCdwaGFzZS12b2NvZGVyLmpzJyk7XG4gICAgICAgIHRoaXMuaW5pdEV4dHJlbWVzRmlsdGVycygpO1xuICAgICAgICB0aGlzLmluaXRGRlROb2RlKCk7XG4gICAgICAgIHRoaXMuaW5pdEdhaW4oKTtcbiAgICAgICAgLy90aGlzLmluaXRGaWx0ZXIoKTtcbiAgICAgICAgYXdhaXQgdGhpcy5pbml0SW5wdXQoKTtcbiAgICAgICAgaWYodGhpcy5pbnB1dCkgdGhpcy5jb25uZWN0QWxsKCk7XG4gICAgfVxuXG4gICAgaW5pdEZpbHRlcigpIHtcbiAgICAgICAgLypjb25zdCBzY2FsZSA9IGNocm9tYXRpY0ZpbHRlcih0aGlzLmJpblRvRnJlcSwgMjIwMDAsIDEyKVxuICAgICAgICB0aGlzLmZpbHRlciA9IG5ldyBOb3RjaEZpbHRlckJhbmsodGhpcy5hdWRpb0NvbnRleHQsIHNjYWxlLnNjYWxlLCAxMDApO1xuICAgICAgICB0aGlzLmZyZXFSZXNwID0gbmV3IEZsb2F0MzJBcnJheSh0aGlzLmZmdFNpemUvMisxKTtcbiAgICAgICAgdGhpcy5mcmVxUmVzcC5maWxsKDEpO1xuICAgICAgICB0aGlzLmZpbHRlci5nZXRGcmVxdWVuY3lSZXNwb25zZShmcmVxcywgdGhpcy5mcmVxUmVzcCk7Ki9cblxuICAgICAgICBjb25zdCBmcmVxcyA9IG5ldyBGbG9hdDMyQXJyYXkodGhpcy5mZnRTaXplLzIrMSk7XG4gICAgICAgIGZyZXFzLmZvckVhY2goKGYsbik9PmZyZXFzW25dPW4qdGhpcy5iaW5Ub0ZyZXEpO1xuICAgICAgICAvLyBjb25zb2xlLmxvZyhzY2FsZS5zY2FsZSlcbiAgICB9XG5cbiAgICBpbml0R2FpbigpIHtcbiAgICAgICAgdGhpcy5tYWluR2FpbiA9IHRoaXMuYXVkaW9Db250ZXh0LmNyZWF0ZUdhaW4oKTtcbiAgICAgICAgdGhpcy5saW1pdGVyID0gdGhpcy5hdWRpb0NvbnRleHQuY3JlYXRlRHluYW1pY3NDb21wcmVzc29yKClcbiAgICAgICAgdGhpcy5saW1pdGVyLnRocmVzaG9sZC52YWx1ZSA9IC00MFxuICAgICAgICB0aGlzLmxpbWl0ZXIucmF0aW8udmFsdWUgPSAyMFxuICAgICAgICB0aGlzLmxpbWl0ZXIuYXR0YWNrLnZhbHVlID0gMC4wMDFcbiAgICAgICAgdGhpcy5saW1pdGVyLnJlbGVhc2UudmFsdWUgPSAwLjAxXG4gICAgfVxuXG4gICAgaW5pdEJ1ZmZlcnMoKSB7XG4gICAgICAgIHRoaXMuYmluVG9GcmVxID0gdGhpcy5hdWRpb0NvbnRleHQuc2FtcGxlUmF0ZSAvIHRoaXMuZmZ0U2l6ZTtcbiAgICAgICAgY29uc3QgbnVtQmlucyA9IHRoaXMuZmZ0U2l6ZS8yICsgMTtcbiAgICAgICAgdGhpcy5mZnRCdWZmZXIgPSBuZXcgRmxvYXQzMkFycmF5KG51bUJpbnMpO1xuICAgICAgICB0aGlzLnBlYWtzQnVmZmVyID0gbmV3IEludDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubnVtUGVha3MgPSAwO1xuICAgICAgICB0aGlzLmZiQnVmZmVyID0gbmV3IEludDMyQXJyYXkobnVtQmlucyk7XG4gICAgICAgIHRoaXMubnVtRmIgPSAwO1xuICAgICAgICB0aGlzLnNsb3BlQnVmZmVyID0gbmV3IEZsb2F0MzJBcnJheShudW1CaW5zKTtcbiAgICAgICAgdGhpcy5maWx0ZXJCdWZmZXIgPSBuZXcgSW50MzJBcnJheShudW1CaW5zKTtcblxuICAgICAgICB0aGlzLm1heERiID0gLTE4MDtcbiAgICB9XG5cbiAgICBub3coKSB7XG4gICAgICAgIHJldHVybiB0aGlzLmF1ZGlvQ29udGV4dC5jdXJyZW50VGltZVxuICAgIH1cblxuICAgIGluaXRFeHRyZW1lc0ZpbHRlcnMoKSB7XG4gICAgICAgIGNvbnN0IGhpcGFzcyA9IHRoaXMuYXVkaW9Db250ZXh0LmNyZWF0ZUJpcXVhZEZpbHRlcigpO1xuICAgICAgICBoaXBhc3MudHlwZSA9IFwiaGlnaHBhc3NcIjsgaGlwYXNzLmZyZXF1ZW5jeS5zZXRWYWx1ZUF0VGltZSh0aGlzLmJpblRvRnJlcSwgdGhpcy5ub3coKSk7XG4gICAgICAgIGhpcGFzcy5RLnNldFZhbHVlQXRUaW1lKDAuNzA3LCB0aGlzLm5vdygpKTtcblxuICAgICAgICBjb25zdCBsb3Bhc3MgPSB0aGlzLmF1ZGlvQ29udGV4dC5jcmVhdGVCaXF1YWRGaWx0ZXIoKTtcbiAgICAgICAgbG9wYXNzLnR5cGUgPSBcImxvd3Bhc3NcIjsgbG9wYXNzLmZyZXF1ZW5jeS5zZXRWYWx1ZUF0VGltZSgxODAwMCwgdGhpcy5ub3coKSk7XG5cbiAgICAgICAgdGhpcy5oaXBhc3MgPSBoaXBhc3M7IHRoaXMubG9wYXNzID0gbG9wYXNzO1xuICAgIH1cblxuICAgIGluaXRGRlROb2RlKCkge1xuICAgICAgICB0aGlzLmZmdE5vZGUgPSBuZXcgQXVkaW9Xb3JrbGV0Tm9kZSh0aGlzLmF1ZGlvQ29udGV4dCwgJ3BoYXNlLXZvY29kZXItcHJvY2Vzc29yJywge3Byb2Nlc3Nvck9wdGlvbnM6IHtibG9ja1NpemU6IHRoaXMuZmZ0U2l6ZX19KTtcbiAgICAgICAgdGhpcy5pbml0V29ya2xldFBvcnQoKTtcbiAgICB9XG5cbiAgICBpbml0V29ya2xldFBvcnQoKSB7XG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5wb3J0Lm9ubWVzc2FnZSA9ICh7ZGF0YX0pID0+IHtcbiAgICAgICAgICAgIHN3aXRjaChkYXRhWzBdKSB7XG4gICAgICAgICAgICAgICAgY2FzZSBcInNwZWN0cnVtXCI6IHRoaXMuZmZ0QnVmZmVyLnNldChkYXRhWzFdKTsgYnJlYWs7XG4gICAgICAgICAgICAgICAgY2FzZSBcInNsb3BlXCI6IHRoaXMuc2xvcGVCdWZmZXIuc2V0KGRhdGFbMV0pOyBicmVhaztcbiAgICAgICAgICAgICAgICBjYXNlIFwiaGlzdG9cIjogdGhpcy5oaXN0b0J1ZmZlci5zZXQoZGF0YVsxXSk7IGJyZWFrO1xuICAgICAgICAgICAgICAgIGNhc2UgXCJwZWFrc1wiOlxuICAgICAgICAgICAgICAgICAgICB0aGlzLm51bVBlYWtzID0gZGF0YVsxXVxuICAgICAgICAgICAgICAgICAgICB0aGlzLnBlYWtzQnVmZmVyLnNldChkYXRhWzJdKVxuICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICBjYXNlIFwiZmJcIjpcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5udW1GYiA9IGRhdGFbMV1cbiAgICAgICAgICAgICAgICAgICAgdGhpcy5mYkJ1ZmZlci5zZXQoZGF0YVsyXSlcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5vbkZiRGF0YSgpO1xuICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICBjYXNlICdtYXhEYic6IHRoaXMubWF4RGIgPSBkYXRhWzFdOyBicmVhaztcbiAgICAgICAgICAgICAgICBjYXNlICdyZWR1Y3Rpb25zJzogdGhpcy5maWx0ZXJCdWZmZXIgPSBkYXRhWzFdOyBicmVhaztcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIC8vY29uc3Qgc3VtID0gTWF0aC5tYXgoLi4uZGF0YSk7XG4gICAgICAgICAgICAvLyBpZihkYXRhWzBdID4gMCkgY29uc29sZS5sb2coXCJGYiBjYW5kaWRhdGVzOlwiLCBkYXRhWzBdLCBkYXRhWzFdKVxuICAgICAgICAgICAgLy8gY29uc29sZS5sb2coZGF0YSlcbiAgICAgICAgfTtcbiAgICAgICAgcmVxdWVzdEFuaW1hdGlvbkZyYW1lKCgpPT50aGlzLnNob3dTcGVjdHJ1bSgpKTtcbiAgICAgICAgcmVxdWVzdEFuaW1hdGlvbkZyYW1lKCgpPT50aGlzLnNob3dTbG9wZSgpKTtcbiAgICB9XG5cbiAgICBjb25uZWN0QWxsKCkge1xuICAgICAgICB0aGlzLmlucHV0XG4gICAgICAgIC5jb25uZWN0KHRoaXMuaGlwYXNzKVxuICAgICAgICAuY29ubmVjdCh0aGlzLmxvcGFzcylcbiAgICAgICAgLmNvbm5lY3QodGhpcy5mZnROb2RlKVxuXG4gICAgICAgIC8vdGhpcy5maWx0ZXIuY29ubmVjdElucHV0KHRoaXMubG9wYXNzKVxuICAgICAgICAvL3RoaXMuZmlsdGVyXG4gICAgICAgIHRoaXMubG9wYXNzLmNvbm5lY3QodGhpcy5tYWluR2FpbikuY29ubmVjdCh0aGlzLmxpbWl0ZXIpLmNvbm5lY3QodGhpcy5hdWRpb0NvbnRleHQuZGVzdGluYXRpb24pXG5cbiAgICAgICAgLy90aGlzLmxvcGFzcy5jb25uZWN0KHRoaXMubWFpbkdhaW4pLmNvbm5lY3QodGhpcy5hdWRpb0NvbnRleHQuZGVzdGluYXRpb24pXG4gICAgfVxuXG4gICAgYXN5bmMgaW5pdElucHV0KCkge1xuICAgICAgICB0cnkge1xuICAgICAgICAgICBjb25zdCBjb25zdHJhaW50cyA9IHthdWRpbzogeyBub2lzZVN1cHByZXNzaW9uOiBmYWxzZSwgZWNob0NhbmNlbGxhdGlvbjogZmFsc2UsIGF1dG9HYWluQ29udHJvbDogZmFsc2UgfX07XG4gICAgICAgICAgIGNvbnN0IHN0cmVhbSA9IGF3YWl0IG5hdmlnYXRvci5tZWRpYURldmljZXMuZ2V0VXNlck1lZGlhKGNvbnN0cmFpbnRzKTtcbiAgICAgICAgICAgdGhpcy5pbnB1dCA9IHRoaXMuYXVkaW9Db250ZXh0LmNyZWF0ZU1lZGlhU3RyZWFtU291cmNlKHN0cmVhbSk7XG4gICAgICAgICAgIGNvbnNvbGUubG9nKFwiR09UIElOUFVUXCIpXG4gICAgICAgICB9IGNhdGNoKGVycikge1xuICAgICAgICAgICAvKiBoYW5kbGUgdGhlIGVycm9yICovXG4gICAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ0FOVCBHRVQgVVNFUiBJTlBVVFwiKVxuICAgICAgICAgfVxuICAgIH1cbiAgICBcblxuICAgIG9uRmJEYXRhKCkge1xuICAgICAgICB0aGlzLnVwZGF0ZUdhaW4oKVxuICAgIH1cblxuICAgIHVwZGF0ZUdhaW4oKSB7XG4gICAgICAgIGlmKHRoaXMubWF4RGIgPD0gLTE4MCkgcmV0dXJuXG4gICAgICAgIGNvbnN0IGdhaW4gPSB0aGlzLm1haW5HYWluLmdhaW47XG4gICAgICAgIGdhaW4uY2FuY2VsQW5kSG9sZEF0VGltZSh0aGlzLm5vdygpKVxuICAgICAgICBjb25zdCBjdXJyZW50R2FpbiA9IDIwICogTWF0aC5sb2cxMCh0aGlzLm1haW5HYWluLmdhaW4udmFsdWUpO1xuICAgICAgICBsZXQgbmV4dEdhaW4gPSBjdXJyZW50R2FpbjtcbiAgICAgICAgaWYodGhpcy5tYXhEYiA+IHRoaXMubWF4Vm9sdW1lKSB7XG4gICAgICAgICAgICBuZXh0R2FpbiArPSAodGhpcy5taW5Wb2x1bWUgLSB0aGlzLm1heERiKSAqIHRoaXMuZ2FpbkRlY3JlbWVudFxuICAgICAgICAgICAgLy9jb25zb2xlLmxvZyhcImdhaW4tLVwiLCB0aGlzLm1heERiLCBuZXh0R2FpbilcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIG5leHRHYWluICs9ICh0aGlzLm1heFZvbHVtZSAtIHRoaXMubWF4RGIpICogdGhpcy5nYWluSW5jcmVtZW50XG4gICAgICAgICAgICAvL2NvbnNvbGUubG9nKFwiZ2FpbisrXCIsdGhpcy5tYXhEYiwgbmV4dEdhaW4pXG4gICAgICAgIH1cbiAgICAgICAgLy9jb25zb2xlLmxvZyh0aGlzLm1heERiLG5leHRHYWluKVxuICAgICAgICBpZihuZXh0R2FpbiAhPSBjdXJyZW50R2Fpbikge1xuICAgICAgICAgICAgbmV4dEdhaW4gPSBNYXRoLnBvdygxMCwgbmV4dEdhaW4qMC4wNSk7XG4gICAgICAgICAgICBpZihpc0Zpbml0ZShuZXh0R2FpbikgJiYgbmV4dEdhaW4gPiAwKVxuICAgICAgICAgICAgICAgIGdhaW4ubGluZWFyUmFtcFRvVmFsdWVBdFRpbWUobmV4dEdhaW4sIHRoaXMubm93KCkpXG4gICAgICAgIH1cblxuICAgIH1cblxuICAgIC8vIGdyYXBoc1xuICAgIHNob3dTcGVjdHJ1bSgpe1xuICAgICAgICByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoKCk9PnRoaXMuc2hvd1NwZWN0cnVtKCkpXG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5wb3J0LnBvc3RNZXNzYWdlKFwiZ2V0U3BlY3RydW1cIik7XG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5wb3J0LnBvc3RNZXNzYWdlKFwiZ2V0UGVha3NcIik7XG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5wb3J0LnBvc3RNZXNzYWdlKFwiZ2V0RmJDYW5kaWRhdGVzXCIpO1xuICAgICAgICB0aGlzLmZmdE5vZGUucG9ydC5wb3N0TWVzc2FnZShcImdldE1heERiXCIpO1xuICAgICAgICB0aGlzLmZmdE5vZGUucG9ydC5wb3N0TWVzc2FnZShcImdldE1hZ1JlZHVjdGlvbnNcIik7XG5cbiAgICAgICAgY29uc3QgYmdDb2xvciA9IHRoaXMubnVtRmIgPiAwID8gR3JhcGgucGl0Y2hBbXBUb0hTTEEodGhpcy5mYkJ1ZmZlclswXSAqIHRoaXMuYmluVG9GcmVxLCB0aGlzLmZmdEJ1ZmZlclt0aGlzLmZiQnVmZmVyWzBdXSkgOiAncmdiYSgxMDAsMTAwLDEwMCwxKSc7XG5cblxuICAgICAgICAvL0dyYXBoLmRyYXdCZyhcImNhbnZhcyNzcGVjdHJ1bVwiLCAncmdiYSgwLDAsMCwwLjUpJylcbiAgICAgICAgR3JhcGguZHJhd0JnKFwiY2FudmFzI3NwZWN0cnVtXCIsIGJnQ29sb3IpXG4gICAgICAgIEdyYXBoLmRyYXdHYWluKFwiY2FudmFzI3NwZWN0cnVtXCIsIFwicmVkXCIsIHRoaXMubWFpbkdhaW4uZ2Fpbi52YWx1ZSwgLTEwMCwgODApXG5cbiAgICAgICAgR3JhcGguZHJhd0xvZ0xpbmUoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5mZnRCdWZmZXIsIC0xMDAsIDApO1xuICAgICAgICAvL0dyYXBoLmRyYXdMb2dMaW5lKFwiY2FudmFzI3NwZWN0cnVtXCIsIHRoaXMuZnJlcVJlc3AsIDAsIDEsIFwicmVkXCIpO1xuXG4gICAgICAgIEdyYXBoLmRyYXdMb2dMaW5lSW52ZXJ0ZWQoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5maWx0ZXJCdWZmZXIsIC0xMDAsIDAsIFwicmdiYSgwLDAsMCwwLjUpXCIpO1xuXG4gICAgICAgIC8vR3JhcGguZHJhd1ZlcnRpY2FscyhcImNhbnZhcyNzcGVjdHJ1bVwiLCB0aGlzLm51bVBlYWtzLCB0aGlzLnBlYWtzQnVmZmVyLCAtODAsIDAsIFwicmdiYSgyNTUsMjU1LDI1NSwxKVwiKVxuICAgICAgICAvL0dyYXBoLmRyYXdWZXJ0aWNhbHMoXCJjYW52YXMjc3BlY3RydW1cIiwgdGhpcy5udW1GYiwgdGhpcy5mYkJ1ZmZlciwgLTgwLCAwLCBcInJnYmEoMjU1LDAsMCwwLjUpXCIpXG4gICAgfVxuXG4gICAgc2hvd1Nsb3BlKG5vZGUpe1xuICAgICAgICByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoKCk9PnRoaXMuc2hvd1Nsb3BlKCkpXG4gICAgICAgIHRoaXMuZmZ0Tm9kZS5wb3J0LnBvc3RNZXNzYWdlKFwiZ2V0U2xvcGVcIik7XG4gICAgICAgIC8vY29uc29sZS5sb2coc2xvcGVCdWZmZXIpXG4gICAgICAgIEdyYXBoLmRyYXdCZyhcImNhbnZhcyNzbG9wZVwiLCAnd2hpdGUnKVxuICAgICAgICBHcmFwaC5kcmF3TG9nTGluZShcImNhbnZhcyNzbG9wZVwiLCB0aGlzLnNsb3BlQnVmZmVyLCAwLCA1MDAsIFwiYmxhY2tcIik7XG4gICAgfVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IFNxdWlkYmFja1dvcmtsZXRQcm9jZXNzIl19
