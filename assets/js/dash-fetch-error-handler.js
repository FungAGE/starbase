/**
 * Shows a general error toast when Dash callback requests fail
 * (413, timeout, 5xx, network errors). Without this, failed requests fail silently.
 */
(function () {
    "use strict";
    var originalFetch = window.fetch;
    if (!originalFetch) return;

    function showToast(title, message) {
        var container = document.getElementById("dash-error-toast-container");
        if (!container) {
            container = document.createElement("div");
            container.id = "dash-error-toast-container";
            container.style.cssText =
                "position:fixed;bottom:20px;right:20px;z-index:10000;display:flex;flex-direction:column;gap:8px;max-width:360px;";
            document.body.appendChild(container);
        }
        var toast = document.createElement("div");
        toast.style.cssText =
            "background:#e03131;color:white;padding:12px 16px;border-radius:8px;box-shadow:0 4px 12px rgba(0,0,0,0.15);font-family:system-ui,sans-serif;";
        var div = document.createElement("div");
        div.textContent = title + " " + message;
        toast.appendChild(div);
        container.appendChild(toast);
        setTimeout(function () {
            toast.style.transition = "opacity 0.3s";
            toast.style.opacity = "0";
            setTimeout(function () {
                if (toast.parentNode) toast.parentNode.removeChild(toast);
            }, 300);
        }, 8000);
    }

    window.fetch = function (input, init) {
        var url = typeof input === "string" ? input : (input && input.url) || "";
        return originalFetch.apply(this, arguments).then(
            function (response) {
                if (
                    response &&
                    !response.ok &&
                    url.indexOf("_dash-update-component") !== -1
                ) {
                    showToast(
                        "Submission failed.",
                        "Something went wrong. Please try again or refresh the page."
                    );
                }
                return response;
            },
            function (error) {
                if (url.indexOf("_dash-update-component") !== -1) {
                    showToast(
                        "Submission failed.",
                        "Something went wrong. Please try again or refresh the page."
                    );
                }
                throw error;
            }
        );
    };
})();
