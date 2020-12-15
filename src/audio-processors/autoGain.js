class AutoGain {
    constructor(audioContext) {
        this.gainIncrement = 0.05;
        this.gainDecrement = 0.05;
        this.desiredLevel = -20;
        this.maxGain = 20;
        this.minGain = -20;

        this.gain = audioContext.createGain();
        this.limiter = audioContext.createDynamicsCompressor()
        this.limiter.threshold.value = -40
        this.limiter.ratio.value = 20
        this.limiter.attack.value = 0.001
        this.limiter.release.value = 0.01

        this.smoothing = 0.99;

        this.now = ()=>audioContext.currentTime
    }

    updateGain(maxDb) {
        if(maxDb <= -180) return
        const currentGain = this.getCurrentValueDb();
        this.gain.gain.cancelScheduledValues(this.now())
        let nextGain = currentGain;
        if(maxDb > this.desiredLevel) {
            nextGain += (this.desiredLevel - maxDb) * this.gainDecrement
        } else {
            nextGain += (this.desiredLevel - maxDb) * this.gainIncrement
        }

        console.log(maxDb, this.desiredLevel, currentGain, nextGain)
        nextGain = this.smoothing * currentGain + (1-this.smoothing) * nextGain

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
    getCurrentValueDb() { return 20 * Math.log10(this.gain.gain.value) }
    connect(destination) { return this.gain.connect(this.limiter).connect(destination) }
}

module.exports = { AutoGain }