function clear_status() {
    document.getElementById("status").innerHTML = "";
  }

  function update_status(message) {
    let status_area = document.getElementById("status");
    let current = status_area.innerHTML;
    status_area.innerHTML = (current + "\n" + message).trim();
  }

  export async function upload() {
    clear_status()
    update_status("Uploading file")

    let formData = new FormData();
    formData.append("file", document.getElementById("file").files[0]);

    let response = await fetch("/upload", {method: "POST", body: formData});
    let ws_id = await response.json();

    let ws = new WebSocket("ws://{{ public_url }}/ws/" + ws_id["ws_id"])
    ws.onmessage = function(event) {
        update_status(event.data)
    }
  }