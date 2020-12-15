const staticSquidback = "squidback-v1"
const assets = [
  "/",
  "/index.html",
  "/style.css",
  "/main.js",
]

self.addEventListener("install", installEvent => {
  installEvent.waitUntil(
    caches.open(staticSquidback).then(cache => {
      cache.addAll(assets)
    })
  )
})

self.addEventListener("fetch", fetchEvent => {
  fetchEvent.respondWith(
    /*caches.match(fetchEvent.request).then(res => {
      return res || fetch(fetchEvent.request)
    })*/
    fetch(fetchEvent.request).catch( () =>
        caches.match(event.request)
    )
  )
})


        
