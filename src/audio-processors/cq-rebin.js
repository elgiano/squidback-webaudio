class CQRebin {
    constructor(sampleRate, fftSize, div, minFreq, maxFreq) {
        this.freqToBin = fftSize / sampleRate;
        this.filter = this.calcSpectralFilter(div, minFreq, maxFreq);
        this.freqs = this.filter.freqs
    }

    calcSpectralFilter(div, minFreq, maxFreq) {
        console.log(`[mel] calculating cq ${div}-tet from ${minFreq} to ${maxFreq} Hz`)
        const numFilters = Math.floor(Math.log2(maxFreq / minFreq) * div);
        const startBins = new Uint16Array(numFilters);
        const freqs = new Float32Array(numFilters);
        const interval = Math.pow(2, 1/div);
        let freq = minFreq;
        // calc frequencies and store center bins in startBins
        for(let i = 0; i < numFilters; ++i) {
            freqs[i] = freq;
            startBins[i] = Math.floor(freq * this.freqToBin)
            freq *= interval;
        }

        const lastBin = Math.floor(freq * this.freqToBin);
        let leftBin = Math.floor(minFreq/interval * this.freqToBin)
        console.log(freqs, startBins)
        const filters = [];
        for (let i = 0; i < numFilters; i++) {
            let centerBin = startBins[i];
            let rightBin = (i < numFilters - 1) ? startBins[i+1] : lastBin;
            const leftRange = (centerBin - leftBin)
            const rightRange = (rightBin - centerBin)
            const newFilter = [];
            // calc triangular coefficients
            console.log(i, leftBin, centerBin, rightBin,)
            for (let bin = leftBin, j=0; bin <= rightBin; bin++, j++) {
                let coeff;
                if(bin < centerBin)
                    coeff = 1.0 - (centerBin - bin) / leftRange;
                else if(bin == centerBin)
                    coeff = 1.0
                else
                    coeff = 1.0 - (bin - centerBin) / rightRange;
                newFilter.push(Math.pow(coeff,2))
            }
            // normalize coefficients to sum to 1
            const coeffSum = newFilter.reduce((s,v)=>s+v,0)
            
            for(let c=0; c < newFilter.length; c++) newFilter[c] /= Math.pow(newFilter.length, 0.8)

            startBins[i] = leftBin;
            leftBin = centerBin;
            filters.push(newFilter)
        }

        // console.log("bins after adjusting", startBins, filters)

        return {freqs, startBins, filters}
    }

    filterSpectrum(magnitudes) {
        // console.log(this.filter)
        return this.filter.filters.map((filter,filter_i)=>{
            const startBin = this.filter.startBins[filter_i];
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
            const startBin = this.filter.startBins[filter_i];
            const sum = minDb + filter.reduce((t,coeff,i)=>
                t + (magnitudes[startBin+i] * coeff)
            ,0) * byteToFloat;
                
            results[filter_i] = sum;
            if(sum > max) {
                max = sum; nMax = filter_i
            } else if(sum < min) min = sum
            avg += sum
        })

        //console.log(this.filter, results)

        return [min, max, nMax, avg / results.length];
    }

    melsToHz(mels) {  return 700 * (Math.exp(mels / 1127) - 1) }
    hzToMels(hertz) { return 1127 * Math.log(1 + hertz/700) }

}

module.exports = {CQRebin}