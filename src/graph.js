const Spline = require("cubic-spline")

class SquidbackGraph {

    drawBg(selector, color) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        canvasCtx.beginPath();
        canvasCtx.fillStyle = color;
        canvasCtx.fillRect(0, 0, canvas.width, canvas.height);
        canvasCtx.closePath();
    }

    drawGain(selector, color, amp=0, min=0, max=1) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const gainHeight = (20 * Math.log10(amp) - min) / (max-min);
        //console.log(amp, gainHeight, min, max)
        canvasCtx.beginPath();
        canvasCtx.fillStyle = color;
        canvasCtx.fillRect(canvas.width-2, canvas.height, canvas.width, -canvas.height * gainHeight);
        canvasCtx.closePath();
    }

    pitchAmpToHSLA(pitch,amp=1,alpha=0.1){
        let octaves = Math.log2(pitch/440);
        let pc = 12 * octaves % 12;
        let h = pc/12*360;
        let l = 10 * (octaves + 5.5);
        let s = 255 + amp;
        let code = `hsla(${h},${s<0?0:s>100?100:s}%,${l<5?5:l>100?100:l}%,${alpha})`;
        //console.log(pitch,pc,octaves,s,h);
        return code;
    }

    drawLine(selector, data, min, max, color="rgba(255,255,255,1)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / data.length);
        let barHeight;
        let x = 0;
        canvasCtx.beginPath();
        canvasCtx.moveTo(x, canvas.height - (data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - (min) ) / (max-min) * canvas.height;
            canvasCtx.lineTo(x + barWidth, canvas.height-barHeight);
            x += barWidth;
        }
        canvasCtx.strokeStyle = color;
        canvasCtx.strokeWidth = 0.1;
        canvasCtx.stroke();
    }

    drawLogLine(selector, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length)
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height);
        canvasCtx.lineTo(0, canvas.height - (data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawSmoothLogLine(selector, data, min, max, color="rgba(255,255,255,0.5)", upsample=3, stroke = false) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length);
        const xs = [...data.keys()];
        data = new Spline(xs, data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, canvas.height);
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = data.at(i); // data[i]
//console.log(i)
            barHeight = ( val - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, canvas.height)
        canvasCtx.closePath();
        if(stroke) {
            canvasCtx.strokeStyle = color
            canvasCtx.stroke();
        } else {
            canvasCtx.fillStyle = color
            canvasCtx.fill();
        }

    }

    drawLogLineInverted(selector, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length - 1)
        canvasCtx.beginPath();
        canvasCtx.moveTo(0, 0);
        canvasCtx.lineTo(0, canvas.height-(data[0] - (min) ) / (max-min) * canvas.height);

        for(let i = 1; i < data.length; ++i){
            barHeight = (data[i] - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawSmoothLogLineInverted(selector, data, min, max, color="rgba(255,255,255,0.5)", upsample = 3) {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        let barHeight;
        const logIToX = canvas.width / Math.log(data.length - 1)
        const xs = [...data.keys()];
        data = new Spline(xs, data);
        canvasCtx.beginPath();
        canvasCtx.moveTo(0.5, 0);
        const inc = 1/upsample;
        for(let i = 0; i < (data.ys.length); i+=inc){
            const val = data.at(i); // data[i]
            barHeight = (val - min) / (max-min) * canvas.height;
            const x = (Math.log(i) + Math.log(i+1)) * logIToX / 2;
            canvasCtx.lineTo(x, canvas.height-barHeight);
        }
        canvasCtx.lineTo(canvas.width, 0)
        canvasCtx.closePath();
        canvasCtx.fillStyle = color
        canvasCtx.fill();
    }

    drawVerticals(selector, numPoints, data, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
        const barWidth = (canvas.width / (data.length - 1));
        const logIToX = canvas.width / Math.log(data.length - 1);
        canvasCtx.fillStyle = color;
        canvasCtx.beginPath();
        
        for(let n = 0; n < numPoints; ++n){
            const i = data[n];
            //console.log(i);
            const x = Math.log(i) * logIToX;
            canvasCtx.fillRect(x, 0, Math.log(i+1) * logIToX - x, canvas.height)
            //canvasCtx.moveTo(x + barWidth*0.5, canvas.height)
            //canvasCtx.lineTo(x + barWidth*0.5,0);
        }
        //canvasCtx.strokeStyle = color
        //canvasCtx.lineWidth = 1;
        //canvasCtx.stroke();
        canvasCtx.fill();
    }

    drawHLine(selector, h, min, max, color="rgba(255,255,255,0.5)") {
        const canvas = document.querySelector(selector);
        const canvasCtx = canvas.getContext("2d");
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