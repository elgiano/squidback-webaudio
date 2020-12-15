"use strict";

if ("serviceWorker" in navigator) {
  window.addEventListener("load", function() {
    navigator.serviceWorker
      .register("./serviceWorker.js")
      .then(res => console.log("[PWA] service worker registered"))
      .catch(err => console.log("[PWA] service worker not registered", err))
  })
}

const SquidbackFilterBankProcess = require('./process/filterBankProcess.js')
const RemoteStream = require('./remoteStream.js')
const { AudioContext } = require('standardized-audio-context')

async function init(online=false) {
    const audioContext = new AudioContext()
    let process = new SquidbackFilterBankProcess(audioContext, document.querySelector("canvas#spectrum"));
    await process.start(2048)

    if (online) {
        let remote = new RemoteStream(audioContext);
        remote.connect(process.output, process.input, "squidback", ()=>{
            updateRemoteGui(remote)
        }) 
    }

    window.addEventListener('resize', ()=>{
        resizeAllCanvas();
        process.clearGraphCache()
    });

    // let process = new SquidbackFFTProcess(audioContext);
    // await process.start(512)
}


function resizeAllCanvas() {
    document.querySelectorAll("canvas").forEach(canvas=>{
      const {clientWidth, clientHeight} = canvas;
     
      if (canvas.width  != clientWidth || canvas.height != clientHeight) {
        canvas.width  = clientWidth; canvas.height = clientHeight;
      }
    })
}

let remoteGuiInterval;
function updateRemoteGui(remote) {
    const container = document.querySelector("#remoteGui")
    const statusElement = document.querySelector("#remoteStatus")
    const counterElement = document.querySelector("#remotePeerCounter")

    if (remote.isConnected) { statusElement.innerText = "Online"; statusElement.classList.add("online") }
    else { statusElement.innerText = "Offline"; statusElement.classList.remove("online") }

    counterElement.innerText = `(${remote.peerCount})`

    container.classList.add("visible")
    if(remoteGuiInterval) clearTimeout(remoteGuiInterval)
    setTimeout(()=>container.classList.remove("visible"), 2500)
}

window.addEventListener('load', ()=>{
    const button = document.querySelector("button#start")
    button.addEventListener("click", async ()=>{
        await init(true);
        button.classList.add("started")
        document.body.requestFullscreen()
    });
    resizeAllCanvas();
})
