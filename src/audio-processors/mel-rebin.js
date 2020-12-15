/*===========================================================================*\
 * Experimental implementation of MFCC.
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Mel-scale and its related coefficients used in
 * human speech analysis.
\*===========================================================================*/

class MelRebin {

    constructor(sampleRate, fftSize, bankCount, lowFrequency, highFrequency) {
      this.constructMelFilterBank(fftSize, bankCount, lowFrequency, highFrequency, sampleRate);
    }

    constructMelFilterBank(fftSize, nFilters, lowF, highF, sampleRate) {

          const filters = [];
          const startBins = new Uint16Array(nFilters);
          const freqs = new Float32Array(nFilters);

        var lowM = hzToMels(lowF),
          highM = hzToMels(highF),
          deltaM = (highM - lowM) / (nFilters+1);

        // Construct equidistant Mel values between lowM and highM.
        for (var i = 0; i < nFilters; i++) {
            freqs[i] = melsToHz(lowM + (i * deltaM));
            //startBins[i] = Math.floor((fftSize+1) * freqs[i] / (sampleRate/2));
            startBins[i] = Math.floor((fftSize/sampleRate) * freqs[i]);
        }

        // Construct one cone filter per bin.
        // Filters end up looking similar to [... 0, 0, 0.33, 0.66, 1.0, 0.66, 0.33, 0, 0...]
        for (var i = 0; i < startBins.length; i++) {
            const newFilter = [];
            var filterRange = (i != startBins.length-1) ? startBins[i+1] - startBins[i] : startBins[i] - startBins[i-1];
            const centerBin = startBins[i];
            const leftBin = startBins[i] - filterRange;
            const rightBin = startBins[i] + filterRange
            for (let bin = leftBin; bin < rightBin; bin++) {
              let coeff;
              // Left edge of cone
              if (bin < centerBin) coeff = 1.0 - ((centerBin - bin) / filterRange);
              // Peak of cone
              else if (bin == centerBin) coeff = 1.0;
              // Right edge of cone
              else if (bin < rightBin) coeff = 1.0 - (bin-centerBin) / filterRange;
              newFilter.push(coeff)
            }
            filters[i] = newFilter;
        }

        // console.log("mel bins", startBins, freqs, filters)

        // Store for debugging.
        this.startBins = startBins;

        this.filters = filters
        this.freqs = freqs
    }

    filter(freqPowers) {
      var ret = [];

      filters.forEach( (filter, fIx) => {
        var tot = 0;
        freqPowers.forEach( (fp, pIx) => {
          tot += fp * filter[pIx];
        });
        ret[fIx] = tot;
      }); 
      return ret;
    }

    filterUint8SpectrumInPlace(magnitudes, results, minDb=0, maxDb=1) {
        // console.log(this.filter)
        let max = -Infinity; let nMax = -1;
        let min = Infinity;
        let avg = 0;
        const byteToFloat = (maxDb-minDb) / 255;

        //console.log(maxDb, minDb, byteToFloat)

        for(let i=0; i< magnitudes.length; i++)
            magnitudes[i] = Math.pow(10, magnitudes[i] * 0.05)

        this.filters.forEach((filter,filter_i)=>{
            const startBin = this.startBins[filter_i];
            let sum = filter.reduce((t,coeff,i)=>
                t + (magnitudes[startBin+i] * coeff)
            ,0);

            sum = 20 * Math.log10(sum)
            results[filter_i] = sum;
            if(sum > max) {
                max = sum; nMax = filter_i
            } else if(sum < min) min = sum
            avg += sum
        })

        // console.log(results)

        return [min, max, nMax, avg / results.length];
    }
}

function melsToHz(mels) {
  return 700 * (Math.exp(mels / 1127) - 1);
}

function hzToMels(hertz) {
  return 1127 * Math.log(1 + hertz/700);
}

module.exports = { MelRebin }

