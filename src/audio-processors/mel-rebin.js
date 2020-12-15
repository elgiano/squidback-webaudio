class MelRebin {
    constructor(sampleRate, fftSize, div, minFreq, maxFreq) {
        this.freqToBin = fftSize / sampleRate;
        this.filter = this.calcSpectralFilter(div, minFreq, maxFreq);
    }

    calcSpectralFilter(div, minFreq, maxFreq) {
        console.log(`[mel] calculating mel ${div}-tet from ${minFreq} to ${maxFreq} Hz`)
        let bins = [], filters = [];
        const numFilters = Math.floor(Math.log2(maxFreq / minFreq) * div);
        const freqs = new Float32Array(numFilters);
        const interval = Math.pow(2, 1/div);
        let freq = minFreq;
        for(let i = 0; i < numFilters; ++i) {
            freqs[i] = freq;
            bins[i] = Math.floor(freq * this.freqToBin)
            freq *= interval;
        }

        const finalBin = Math.floor(freq * this.freqToBin);
        let leftBin = Math.floor(minFreq/interval * this.freqToBin)
        for (let i = 0; i < numFilters; i++) {
            let centerBin = bins[i];
            let rightBin = (i < numFilters - 1) ? bins[i+1] : finalBin;
            filters[i] = [];
            for (let bin = leftBin, j=0; bin <= rightBin; bin++, j++) {
                if(bin < centerBin)
                    filters[i][j] = 1.0 - (centerBin - bin) / (centerBin - leftBin);
                else if(bin == centerBin)
                    filters[i][j] = 1.0
                else
                    filters[i][j] = 1.0 - (bin - centerBin) / (rightBin - centerBin);
            }
            bins[i] = leftBin;
            leftBin = centerBin;
        }

        return {freqs, bins, filters}
    }

    filterSpectrum(magnitudes) {
        // console.log(this.filter)
        return this.filter.filters.map((filter,filter_i)=>{
            const startBin = this.filter.bins[filter_i];
            return filter.reduce((t,coeff,i)=>
                t + (magnitudes[startBin+i] * coeff)
            ,0);
        })
    }

    filterUint8SpectrumInPlace(magnitudes, results, minDb=0, maxDb=1) {
        // console.log(this.filter)
        let max = -Infinity; let nMax = -1;
        let min = Infinity;
        let avg = 0;
        const byteToFloat = (maxDb-minDb) / 255;

        this.filter.filters.forEach((filter,filter_i)=>{
            const startBin = this.filter.bins[filter_i];
            const sum = minDb + filter.reduce((t,coeff,i)=>
                t + (magnitudes[startBin+i] * coeff)
            ,0) * byteToFloat;
                
            results[filter_i] = sum;
            if(sum > max) {
                max = sum; nMax = filter_i
            } else if(sum < min) min = sum
            avg += sum
        })
        return [min, max, nMax, avg / results.length];
    }

    melsToHz(mels) {  return 700 * (Math.exp(mels / 1127) - 1) }
    hzToMels(hertz) { return 1127 * Math.log(1 + hertz/700) }

}

module.exports = {MelRebin}