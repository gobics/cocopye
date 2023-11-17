function update_status(message, state) {
    let status = document.getElementById("status");
    status.innerHTML = "";

    if (state == "progress") {
        let spinner = document.createElement("span");
        spinner.classList.add("loading", "loading-spinner", "loading-sm", "mr-4");
        status.appendChild(spinner);
        status.appendChild(document.createTextNode(message));
    }

    if (state == "result") {
        document.getElementById("completeness").innerHTML = message.completeness;
        document.getElementById("contamination").innerHTML = message.contamination;
        document.getElementById("stage").innerHTML = message.stage;
        document.getElementById("taxonomy").innerHTML = message.taxonomy;

        let div = document.createElement("div");
        div.classList.add("text-green-800");
        div.appendChild(document.createTextNode("Finished."));
        status.appendChild(div);

        let results = document.getElementById("results");
        results.classList.remove("opacity-30");
        results.classList.add("opacity-100");
        document.getElementById("upload-btn").removeAttribute("disabled");
    }

    if (state == "error") {
        let div = document.createElement("div");
        div.classList.add("text-red-700");
        div.appendChild(document.createTextNode(message));
        status.appendChild(div);
    }
}

async function upload() {
    if (document.getElementById("file").files[0] === undefined) {
        update_status("No file selected.", "error");
        return;
    }
    else {
        update_status("Uploading file", "progress");
        let results = document.getElementById("results");
        results.classList.remove("opacity-100");
        results.classList.add("opacity-30");
        document.getElementById("upload-btn").setAttribute("disabled", "true");
    }

    let formData = new FormData();
    formData.append("file", document.getElementById("file").files[0]);

    let response = await fetch("/upload", { method: "POST", body: formData });
    let ws_id = await response.json();

    let ws = new WebSocket("ws://" + window.location.host + "/ws/" + ws_id["ws_id"])
    ws.onmessage = function (event) {
        let received = JSON.parse(event.data);
        update_status(received.content, received.status)
    }
}

document.getElementById("upload-btn").addEventListener('click', upload);