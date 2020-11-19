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