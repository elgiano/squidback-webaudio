class NotchFilterBank {
    constructor(audioContext, scale, q) {
        this.audioContext = audioContext
        this.scale = scale;
        this.q = q; 
        this.filters = [];
        //|| calculateQs(this.scale); //console.log(q);

        this.makeFiltersFromScale(this.scale, this.q);

        this.totalFreqResponse = new Float32Array(this.filters.length);
        this.singleFreqResponse = new Array(this.filters.length)
        for(let i = 0; i < this.filters.length; i++)
            this.singleFreqResponse[i] = new Float32Array(this.filters.length);
    }

    makeFiltersFromScale() {
        this.scale.forEach((freq,i)=>{
          const filter = this.audioContext.createBiquadFilter();
          const type = i == 0 ? "lowshelf" : i == freq.length - 1 ? "highshelf" : "peaking";
          filter.type = type; filter.frequency.setValueAtTime(freq, this.audioContext.currentTime);
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

    updateFreqResponses(freqs) {
        this.totalFreqResponse.fill(1);
        const _ = new Float32Array(freqs.length);
        this.filters.forEach((f,i)=>{
            f.getFrequencyResponse(freqs,this.singleFreqResponse[i],_)
            this.singleFreqResponse[i].forEach((t,n)=>this.totalFreqResponse[n]*=t)
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