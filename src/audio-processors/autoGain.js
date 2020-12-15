class AutoGain {
    constructor(audioContext) {
        this.gainIncrement = 0.001;
        this.gainDecrement = 0.001;
        this.minVolume = -40;
        this.maxVolume = -10;
        this.maxGain = 100;
        this.minGain = -100;

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
        this.gain.gain.cancelScheduledValues(this.now())

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