const Spline = require("cubic-spline")
const smooth = require("array-smooth")

class SquidbackGraph {

    constructor(canvas){ 
        this.canvas = canvas;
        this.canvasCtx = canvas.getContext("2d");
        this.cache = {} 
    }

    setupFreqAxis(min, max, binToFreq) {
        this.minFreq = min; this.maxFreq = max;
        this.rLogFreqRange = 1 / Math.log(max/min)
        this.binToFreq = binToFreq;
    }
    freqToUniX(freq) {
        return Math.log(freq/this.minFreq) * this.rLogFreqRange 
    }

    fillLine(data, min, max, color="rgba(255,255,255,1)", inverted=false) {
        const canvas = this.canvas;
        const canvasCtx = this.canvasCtx;

        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height)

        let barHeight = (data[0] - min ) / (max-min) * canvas.height;
        const barWidth = (canvas.width / (data.length - 1));
        let x = 0;
        for(let i = 0; i < (data.length - 1); i++){
            barHeight = (data[i] - min ) / (max-min) * canvas.height;
            canvasCtx.lineTo(x, canvas.height - barHeight);
            x += barWidth;
        }
        canvasCtx.lineTo(canvas.width, canvas.height - barHeight)
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color;
        canvasCtx.fill();
    }

    drawFreqFFT(data, min, max, color="rgba(255,255,255,0.5)", cacheId, fill=true, inverted = false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height * (1-inverted));
        let init = false;

        let barHeight;
        for(let i = 1; i < data.length; ++i){
            const freq = i * this.binToFreq;
            if (freq < this.minFreq) continue
            if (freq > this.maxFreq) break

            barHeight = (data[i] - min) / (max-min) * canvas.height;
            if (!init) { canvasCtx.lineTo(0, canvas.height - barHeight); init = true };

            canvasCtx.lineTo(this.freqToUniX(freq) * canvas.width, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height * (1-inverted))
        canvasCtx.closePath();
        if(fill) {
            canvasCtx.fillStyle = color
            canvasCtx.fill();
        } else {
            canvasCtx.strokeStyle = color
            canvasCtx.stroke();
        }
    }

    drawSmoothCQDB(data, min, max, color="rgba(255,255,255,0.5)",  freqs, upsample = 3, inverted = false) {
        const canvas = this.canvas;
        const canvasCtx = this.canvasCtx;

        const barWidth = canvas.width / (data.length - 1) / upsample;
        let barHeight;
        data = new Spline([...data.keys()], data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, -1);
        let x = 0;
        let val =  20 * Math.log10(data.at(0));
        barHeight = (val - min) / (max-min) * canvas.height;
        canvasCtx.lineTo(0, canvas.height - barHeight);
        const inc = 1/upsample;
        for(let i = 0; i < data.ys.length-1; i+=inc){
            val = 20 * Math.log10(data.at(i));
            barHeight = (val - min) / (max-min) * canvas.height;
            x += barWidth;
            canvasCtx.lineTo(x, canvas.height-barHeight);
            //canvasCtx.fillRect(x, 0, 0.1, canvas.height);
        }
        canvasCtx.lineTo(canvas.width, canvas.height-barHeight);
        canvasCtx.lineTo(canvas.width, -1)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawGrayscaleVerticals(data, alpha=0.005) {
        const canvasCtx = this.canvasCtx;
        const barWidth = (this.canvas.width / (data.length - 1));
        canvasCtx.beginPath();
        let x = -0.5 * barWidth;
        for(let n = 0; n < data.length; ++n){
            let v = (Math.min(Math.max(data[n], -100), 100) + 100) * 0.005;
            v *= 255;
            canvasCtx.fillStyle = `rgba(${v},${v},${v},${alpha})`;
            canvasCtx.fillRect(x, 0, barWidth, this.canvas.height)
            x += barWidth;
        }
        canvasCtx.closePath();
    }

    drawBg(color) {
        this.canvasCtx.beginPath();
        this.canvasCtx.fillStyle = color;
        this.canvasCtx.fillRect(0, 0, this.canvas.width, this.canvas.height);
        this.canvasCtx.closePath();
    }

    drawGain(color, db=0, min=0, max=1) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx
        const gainHeight = db / (db < 0 ? (-min) : max) * canvas.height / 2;
        canvasCtx.fillStyle = "rgba(30%,30%,30%,1)";
        canvasCtx.fillRect(canvas.width-2, 0, 2, canvas.height);
        canvasCtx.fillStyle = color;
        canvasCtx.fillRect(canvas.width-2, canvas.height / 2, 2,  - gainHeight);
    }

    pitchAmpToHSLA(pitch,amp=255,alpha=0.1){
        let octaves = Math.log2(pitch/440);
        let pc = 12 * octaves % 12;
        let h = pc/12*360;
        let l = 10 * (octaves + 5.5);
        let s = 255 + amp;
        let code = `hsla(${h},${s<0?0:s>100?100:s}%,${l<5?5:l>100?100:l}%,${alpha})`;
        //console.log(pitch,pc,octaves,s,h);
        return code;
    }

    drawLineInverted(data, min, max, color="rgba(255,255,255,1)"){
        this.drawLine(data, min, max, color="rgba(255,255,255,1)", true);
    }
    drawLine(data, min, max, color="rgba(255,255,255,1)", inverted=false) {
        const barWidth = this.canvas.width / (data.length - 1);
        let barHeight;
        let x = 0;
        this.canvasCtx.beginPath();
        this.canvasCtx.moveTo(x, this.canvas.height - (data[0] - (min) ) / (max-min) * this.canvas.height);

        for(let i = 0; i < data.length; ++i){
            barHeight = (data[i] - (min) ) / (max-min) * this.canvas.height;
            this.canvasCtx.lineTo(x, this.canvas.height-barHeight);
            x += barWidth;
        }
        this.canvasCtx.strokeStyle = color;
        this.canvasCtx.strokeWidth = 0.1;
        this.canvasCtx.stroke();
    }

    drawLinVerticals(numPoints, data, color="rgba(255,255,255,0.005)") {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx
        const barWidth = (canvas.width / (data.length - 1));
        canvasCtx.beginPath();
        canvasCtx.fillStyle = color;
        for(let n = 0; n < numPoints; ++n){
            const i = data[n];
            canvasCtx.fillRect(barWidth*(i-0.25), 0, barWidth*0.5, canvas.height)
        }
        canvasCtx.closePath();
    }

    drawSmoothLine(data, min, max, color="rgba(255,255,255,1)", upsample=5, inverted=false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        data = new Spline([...data.keys()], data);

        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height)

        let barHeight = (data.at(0) - min ) / (max-min) * canvas.height;
        const barWidth = (canvas.width / (data.ys.length - 1) / upsample);
        let x = 0;
        const inc = 1/upsample;


        for(let i = 0; i < (data.ys.length - 1); i+=inc){
            barHeight = (data.at(i) - min ) / (max-min) * canvas.height;
            canvasCtx.lineTo(x, canvas.height - barHeight);
            x += barWidth;
        }
        canvasCtx.lineTo(canvas.width, barHeight)
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color;
        canvasCtx.fill();
    }

    drawSmoothLineDB(data, min, max, color="rgba(255,255,255,0.5)",  freqs, upsample = 3, inverted = false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        // xMin *= canvas.width; xMax *= canvas.width;
        const xMin = Math.log(freqs[0]/40) / Math.log(500) * canvas.width 
        const xMax = Math.log(freqs[freqs.length - 1]/20) / Math.log(1000) * canvas.width 

        const barWidth = (xMax-xMin) / (data.length * upsample - 1);
        let barHeight;
        const zeroLevel = 0;//canvas.height * (1-inverted) - (0 - min) / (max-min) * canvas.height * (1-inverted);
        data = new Spline([...data.keys()], data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, zeroLevel);
        canvasCtx.lineTo(xMin, zeroLevel);
        let x = xMin + barWidth/2;
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = 20 * Math.log10(data.at(i)); // data[i]
            barHeight = (val - min) / (max-min) * canvas.height;
            x += barWidth;
            canvasCtx.lineTo(x, canvas.height-barHeight);
            //canvasCtx.fillRect(x, 0, 0.1, canvas.height);
        }
        canvasCtx.lineTo(xMax, zeroLevel);
        canvasCtx.lineTo(canvas.width, zeroLevel)
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawSmoothLineInverted(data, min, max, color="rgba(255,255,255,0.5)",  xMin=0, xMax=1, upsample = 3) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        xMin *= canvas.width; xMax *= canvas.width;

        const barWidth = (xMax-xMin) / (data.length * upsample - 1);
        let barHeight;
        const zeroLevel = canvas.height - (0 - min) / (max-min) * canvas.height;
        data = new Spline([...data.keys()], data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, 0);
        canvasCtx.lineTo(0, zeroLevel);
        canvasCtx.lineTo(xMin, zeroLevel);
        let x = xMin;
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = data.at(i); // data[i]
            barHeight = (val - min) / (max-min) * canvas.height;
            x += barWidth;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(xMax, zeroLevel);
        canvasCtx.lineTo(canvas.width, zeroLevel)
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    precomputeXs(width, len, inc=1) {
        const logIToX = width / Math.log(/*len - 1*/500)
        let xs = [];
        for(let i = 0; i < len; i+=inc)
            xs[i] = (Math.log(i * 44100 / 1024 / 40) /*+ Math.log(i+1)) / 2*/ * logIToX);
        return xs;
    }

    getXs(width, len, inc=1, cacheId) {
        let xs;
        if(cacheId) {
            if(!this.cache[cacheId])
                this.cache[cacheId] = this.precomputeXs(width, len, inc)
            xs = this.cache[cacheId];
        } else {
            xs = this.precomputeXs(width, len, inc)
        }
        return xs;
    }

    clearCache() { this.cache = {} }

    drawLogLineInverted(data, min, max, color="rgba(255,255,255,0.5)", cacheId){
        this.drawLogLine(data, min, max, color, cacheId, true)
    }
    drawLogLine(data, min, max, color="rgba(255,255,255,0.5)", cacheId, inverted = false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        const xs = this.getXs(canvas.width, data.length, 1, cacheId);

        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height * (1-inverted));
        canvasCtx.lineTo(0, canvas.height - (data[0] - (min) ) / (max-min) * canvas.height);

        let barHeight;
        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - min) / (max-min) * canvas.height;
            canvasCtx.lineTo(xs[i], canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height * (1-inverted))
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }


    drawSmoothLogLineInverted(data, min, max, color="rgba(255,255,255,0.5)", upsample = 3, cacheId) {
        this.drawSmoothLogLine(data, min, max, color, upsample, cacheId, true)
    }
    drawSmoothLogLine(data, min, max, color="rgba(255,255,255,0.5)", upsample = 3, cacheId, inverted = false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        data = new Spline([...data.keys()], data);
        const inc = 1/upsample;
        const xs = this.getXs(canvas.width, data.ys.length, inc, cacheId)

        canvasCtx.beginPath();
        canvasCtx.moveTo(0.5, canvas.height * (1-inverted));
        let barHeight;
        for(let i = 0; i < data.ys.length; i+=inc){
            const val = data.at(i);
            barHeight = (val - min) / (max-min) * canvas.height;
            canvasCtx.lineTo(xs[i], canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height * (1-inverted))
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawAvgLogLineInverted(data, min, max, color="rgba(255,255,255,0.5)", smoothWindow = 3, cacheId) {
        this.drawAvgLogLine(data, min, max, color, smoothWindow, cacheId, true)
    }
    drawAvgLogLine(data, min, max, color="rgba(255,255,255,0.5)", smoothWindow = 3, cacheId, inverted = false) {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

        data = smooth(data, smoothWindow);
        const xs = this.getXs(canvas.width, data.length, 1, cacheId)

        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height * (1-inverted));
        canvasCtx.lineTo(0, canvas.height-(data[0] - min) / (max-min) * canvas.height);
        let barHeight;
        for(let i = 0; i < data.length; i++){
            const val = data[i];
            barHeight = (val - min) / (max-min) * canvas.height;
            canvasCtx.lineTo(xs[i], canvas.height-barHeight);
            //canvasCtx.fillRect(xs[i], 0, 0.1, canvas.height);

        }
        canvasCtx.lineTo(canvas.width, canvas.height * (1-inverted))
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawVerticals(numPoints, data, color="rgba(255,255,255,0.5)") {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx
        //const logIToX = canvas.width / Math.log(data.length - 1);
        const logIToX = canvas.width / Math.log(/*len - 1*/500)

        canvasCtx.fillStyle = color;
        canvasCtx.beginPath();
        const i2f =  44100 / 1024 / 40;
        
        for(let n = 0; n < numPoints; ++n){
            const i = data[n];
            //console.log(i);
            const x = Math.log(i * i2f) * logIToX;
            canvasCtx.fillRect(x, 0, Math.log((i+1)*i2f) * logIToX - x, canvas.height)
            //canvasCtx.moveTo(x + barWidth*0.5, canvas.height)
            //canvasCtx.lineTo(x + barWidth*0.5,0);
        }
        //canvasCtx.strokeStyle = color
        //canvasCtx.lineWidth = 1;
        //canvasCtx.stroke();
        canvasCtx.fill();
    }

    drawHLine(h, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = this.canvas; const canvasCtx = this.canvasCtx

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