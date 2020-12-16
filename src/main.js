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

    if (online) startRemote(process)

    window.addEventListener('resize', ()=>{
        resizeAllCanvas();
        process.clearGraphCache()
    });

    // let process = new SquidbackFFTProcess(audioContext);
    // await process.start(512)
}

function startRemote(process) {
    let remote = new RemoteStream(process.audioContext);
    remote.connect(process.output, process.audioContext.destination, "squidback", ()=>{
        updateRemoteGui(remote)
    });
    updateRemoteGui(remote)
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

function tryWakeLock() {
    // The wake lock sentinel.
    let wakeLock = null;

    // Function that attempts to request a wake lock.
    const requestWakeLock = async () => {
      try {
        wakeLock = await navigator.wakeLock.request('screen');
        wakeLock.addEventListener('release', () => {
          console.log('Wake Lock was released');
        });
        console.log('Wake Lock is active');
      } catch (err) {
        console.error(`${err.name}, ${err.message}`);
      }
    };
}

window.addEventListener('load', ()=>{
    const button = document.querySelector("button#start")
    button.addEventListener("click", async ()=>{
        await init(true);
        button.classList.add("started")
        document.body.requestFullscreen()
    });
    resizeAllCanvas();
    tryWakeLock();
})
