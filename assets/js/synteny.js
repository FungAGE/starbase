// Synteny visualization helper functions
// These functions support the ClusterMap.js integration in Dash

// SVG serialization function (adapted from clinker.js)
function serialiseSVG(svg) {
    /* Saves the figure to SVG in its current state.
     * Clones the provided SVG and sets the width/height of the clone to the
     * bounding box of the original SVG. Thus, downloaded figures will be sized
     * correctly.
     * This function returns a new Blob, which can then be downloaded.
     */
    let node = svg.node().cloneNode(true);

    // Remove all hidden elements
    d3.select(node)
        .selectAll(".hidden")
        .remove();

    // Draw filtered node in hidden <div> to calculate real bounding box
    let container = svg.node().parentNode;
    let div = d3.select(container)
        .append("div")
        .attr("visibility", "hidden");
    div.append(() => node);
    let bbox = div.select("g.clusterMapG").node().getBBox();
    div.remove();

    const xmlns = "http://www.w3.org/2000/xmlns/";
    const xlinkns = "http://www.w3.org/1999/xlink";
    const xhtml = "http://www.w3.org/1999/xhtml";
    const svgns = "http://www.w3.org/2000/svg";

    node.setAttribute("width", bbox.width);
    node.setAttribute("height", bbox.height);
    node.setAttributeNS(xmlns, "xmlns", svgns);
    node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
    node.setAttributeNS(xmlns, "xmlns:xhtml", xhtml);

    // Adjust x/y of <g> to account for axis/title position.
    // Replaces the transform attribute, so drag/zoom is ignored.
    d3.select(node)
        .select("g.clusterMapG")
        .attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`);

    const serializer = new window.XMLSerializer;
    const string = serializer.serializeToString(node);
    return new Blob([string], {type: "image/svg+xml"});
}

// Download function
function downloadFile(blob, filename) {
    /* Downloads a given blob to filename.
     * This function appends a new anchor to the document, which points to the
     * supplied blob. The anchor.click() method is called to trigger the download,
     * then the anchor is removed.
     */
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

// Initialize ClusterMap visualization
function initializeClusterMap(containerId, data, config = {}) {
    console.log("Initializing ClusterMap for container:", containerId);
    
    const container = d3.select(`#${containerId}`);
    if (container.empty()) {
        console.error("Container not found:", containerId);
        return null;
    }

    // Clear existing content
    container.selectAll("*").remove();

    // Default configuration
    const defaultConfig = {
        plot: {
            scaleFactor: 30,
            scaleGenes: true,
        },
        cluster: {
            spacing: 50,
            alignLabels: true,
        },
        gene: {
            label: {
                show: false,
            }
        },
        legend: {
            show: true,
        },
        link: {
            show: true,
            threshold: 0.3,
        }
    };

    // Merge configurations
    const finalConfig = Object.assign({}, defaultConfig, config);

    try {
        // Create ClusterMap instance
        const chart = ClusterMap.ClusterMap().config(finalConfig);
        
        // Apply the chart to the container
        container.datum(data).call(chart);
        
        console.log("ClusterMap initialized successfully");
        return chart;
    } catch (error) {
        console.error("Error initializing ClusterMap:", error);
        container.html(`<div style="color: red; padding: 20px;">Error initializing visualization: ${error.message}</div>`);
        return null;
    }
}

// Update ClusterMap configuration
function updateClusterMapConfig(containerId, config) {
    console.log("Updating ClusterMap config for:", containerId);
    
    const container = d3.select(`#${containerId}`);
    if (container.empty()) {
        console.error("Container not found:", containerId);
        return;
    }

    // Try to get existing chart instance (this might need adjustment based on how we store it)
    const chart = container.datum() ? ClusterMap.ClusterMap().config(config) : null;
    
    if (chart && container.datum()) {
        try {
            container.call(chart);
            console.log("ClusterMap config updated successfully");
        } catch (error) {
            console.error("Error updating ClusterMap config:", error);
        }
    }
}

// Export SVG from ClusterMap
function exportClusterMapSVG(containerId, filename = 'synteny.svg') {
    console.log("Exporting SVG from:", containerId);
    
    const svg = d3.select(`#${containerId} svg.clusterMap`);
    if (svg.empty()) {
        console.error("No SVG found in container:", containerId);
        return;
    }

    try {
        const blob = serialiseSVG(svg);
        downloadFile(blob, filename);
        console.log("SVG exported successfully");
    } catch (error) {
        console.error("Error exporting SVG:", error);
    }
}

// Reset ClusterMap view (zoom/pan)
function resetClusterMapView(containerId) {
    console.log("Resetting view for:", containerId);
    
    const svg = d3.select(`#${containerId} svg.clusterMap`);
    if (svg.empty()) {
        console.error("No SVG found in container:", containerId);
        return;
    }

    try {
        // Reset zoom and pan to default
        const zoom = d3.zoom();
        svg.call(zoom.transform, d3.zoomIdentity.translate(20, 50).scale(1.2));
        console.log("View reset successfully");
    } catch (error) {
        console.error("Error resetting view:", error);
    }
}

// Make functions available globally for Dash callbacks
window.syntenyViz = {
    initialize: initializeClusterMap,
    updateConfig: updateClusterMapConfig,
    exportSVG: exportClusterMapSVG,
    resetView: resetClusterMapView,
    serialiseSVG: serialiseSVG,
    downloadFile: downloadFile
};