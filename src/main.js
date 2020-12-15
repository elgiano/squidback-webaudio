"use strict";

const SquidbackFilterBankProcess = require('./process/filterBankProcess.js')
const RemoteStream = require('./remoteStream.js')

async function init() {
    const audioContext = new AudioContext()
    let process = new SquidbackFilterBankProcess(audioContext, document.querySelector("canvas#spectrum"));
    await process.start(2048)

    let remote = new RemoteStream(audioContext);
    remote.connect(process.output, process.input)

    window.addEventListener('resize', ()=>{
        resizeAllCanvas();
        process.clearGraphCache()
    });

    // let process = new SquidbackFFTProcess(audioContext);
    // await process.start(512)
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
        document.body.requestFullscreen()
    });
    resizeAllCanvas();
});


