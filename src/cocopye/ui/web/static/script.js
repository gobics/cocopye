// Copyright 2023 Niklas Birth
//
// This file is part of CoCoPyE.
//
// CoCoPyE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CoCoPyE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CoCoPyE. If not, see <https://www.gnu.org/licenses/>.

function update_status(message, state) {
    let status = document.getElementById("status");
    status.innerHTML = "";

    if (state === "progress") {
        let spinner = document.createElement("span");
        spinner.classList.add("loading", "loading-spinner", "loading-sm", "mr-4");
        status.appendChild(spinner);
        status.appendChild(document.createTextNode(message));
    }

    if (state === "result") {
        let error = document.getElementById("stage-0-error");
        error.style.display = "none";

        let method = "-";
        if(message.stage === 2) method = "markers";
        if(message.stage === 3) method = "markers, neural network";

        document.getElementById("completeness").innerHTML = message.completeness;
        document.getElementById("contamination").innerHTML = message.contamination;
        document.getElementById("stage").innerHTML = method;
        document.getElementById("taxonomy").innerHTML = message.taxonomy;

        let div = document.createElement("div");
        div.classList.add("text-green-800");
        div.appendChild(document.createTextNode("Finished."));
        status.appendChild(div);

        let results = document.getElementById("results");
        document.getElementById("upload-btn").removeAttribute("disabled");

        if(message.stage === 1) {
            error.style.display = "block";

            document.getElementById("completeness").innerHTML ="-";
            document.getElementById("contamination").innerHTML = "-";
            document.getElementById("stage").innerHTML = "-";
            document.getElementById("taxonomy").innerHTML ="-";

            return;
        }

        results.classList.remove("opacity-30");
        results.classList.add("opacity-100");
    }

    if (state === "error") {
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