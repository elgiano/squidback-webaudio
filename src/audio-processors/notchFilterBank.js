class NotchFilterBank {
    constructor(audioContext, scale, q) {
        this.audioContext = audioContext
        this.scale = scale;
        this.q = q; 
        this.filters = [];

        //|| calculateQs(this.scale); //console.log(q);

        this.makeFiltersFromScale(this.scale, this.q);
    }

    makeFiltersFromScale() {
        this.scale.forEach((freq,i)=>{
          const filter = this.audioContext.createBiquadFilter();
          filter.type = "peaking"; filter.frequency.setValueAtTime(freq, this.audioContext.currentTime);
          filter.gain.value = 1;
          filter.Q.value = typeof this.q == 'object' ? this.q[i] : this.q ;
          this.filters.push(filter);
        })

        this.filters.forEach((f,n) => { if(n<this.filters.length-1)f.connect(this.filters[n+1]) });
    }

    getFrequencyResponse(freqs, response) {
        response.fill(1);
        const tmp = new Float32Array(freqs.length);
        const _ = new Float32Array(freqs.length);
        this.filters.forEach(f=>{
            f.getFrequencyResponse(freqs,tmp,_)
            tmp.forEach((t,n)=>response[n]*=t)
        })
    }

    forEachFreqResponse(freqs, action) {
        const tmp = new Float32Array(freqs.length);
        const _ = new Float32Array(freqs.length);
        this.filters.forEach(f=>{
            f.getFrequencyResponse(freqs,tmp,_)
            action(tmp)
        })
    }

    connectInput(input) {
        input.connect(this.filters[0])
        return this.filters[this.filters.length-1]
    }
    connect(output) {
        this.filters[this.filters.length-1].connect(output)
        return output
    }

    setGains(gains) {
        gains.forEach((g,n)=>{
            this.filters[n].gain.value = g
        })
    }

}

module.exports = { NotchFilterBank }