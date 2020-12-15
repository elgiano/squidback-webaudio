class AutoGain {
    constructor(audioContext) {
        this.gainIncrement = 0.05 / 100;
        this.gainDecrement = 0.05 / 100;
        this.desiredLevel = -40;
        this.tolerance = 5;
        this.maxGain = 20;
        this.minGain = -20;

        this.gain = audioContext.createGain();
        this.limiter = audioContext.createDynamicsCompressor()
        this.limiter.threshold.value = -40
        this.limiter.ratio.value = 20
        this.limiter.attack.value = 0.001
        this.limiter.release.value = 0.01

        this.smoothing = 0//0.99;

        this.now = ()=>audioContext.currentTime
    }

    updateGain(maxDb) {
        if(maxDb <= -180) return
        const currentGain = this.getCurrentValueDb();
        this.gain.gain.cancelScheduledValues(this.now())
        let nextGain = currentGain;
        const diff = this.desiredLevel - maxDb
        if(Math.abs(diff) > this.tolerance)
            nextGain += diff * (diff > 0 ? this.gainIncrement : this.gainDecrement)

        // console.log(maxDb, this.desiredLevel, diff, nextGain)
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