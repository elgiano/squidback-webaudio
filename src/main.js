"use strict";

const SquidbackWorkletProcess = require('./workletProcess.js')
const SquidbackFFTProcess = require('./fftProcess.js')

async function init() {
    const audioContext = new AudioContext()
    if (audioContext.audioWorklet === undefined) {
        handleNoWorklet();
        return;
    }
    let process = new SquidbackFFTProcess(audioContext, 512);
    await process.start()
}

function handleNoWorklet() {
    let $noWorklet = document.querySelector("#no-worklet");
    $noWorklet.style.display = 'block';
    let $timeline = document.querySelector(".timeline");
    $timeline.style.display = 'none';
    let $controls = document.querySelector(".controls");
    $controls.style.display = 'none';
}

function resize(canvas) {
  // Lookup the size the browser is displaying the canvas.
  var displayWidth  = canvas.clientWidth;
  var displayHeight = canvas.clientHeight;
 
  // Check if the canvas is not the same size.
  if (canvas.width  != displayWidth ||
      canvas.height != displayHeight) {
 
    // Make the canvas the same size
    canvas.width  = displayWidth;
    canvas.height = displayHeight;
  }
}

function resizeAllCanvas() {
    document.querySelectorAll("canvas").forEach(canvas=>resize(canvas))
}


window.addEventListener('load', ()=>{
    const button = document.querySelector("button#start")
    button.addEventListener("click", ()=>{
        init();
        button.classList.add("started")
    });
    resizeAllCanvas();
});

window.addEventListener('resize', ()=>{
    resizeAllCanvas()
});
