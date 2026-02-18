// Enhanced Synteny visualization helper functions
// These functions support the ClusterMap.js integration in Dash with improved tooltips and interactivity

// SVG serialization function (adapted from clinker.js)
function serialiseSVG(svg) {
    /* Saves the figure to SVG in its current state.
     * Clones the provided SVG and sets the width/height of the clone to the
     * bounding box of the original SVG.
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

    // Adjust x/y of <g> to account for axis/title position
    d3.select(node)
        .select("g.clusterMapG")
        .attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`);

    const serializer = new window.XMLSerializer;
    const string = serializer.serializeToString(node);
    return new Blob([string], {type: "image/svg+xml"});
}

// Download function
function downloadFile(blob, filename) {
    /* Downloads a given blob to filename */
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

// Enhanced tooltip creation
function createGeneTooltip(gene) {
    /**
     * Create a rich HTML tooltip for a gene with all metadata
     */
    if (!gene || !gene.metadata) {
        return "No data available";
    }
    
    let html = '<div class="gene-tooltip-content">';
    
    // Header with gene name
    html += `<div class="tooltip-header">
        <strong style="font-size: 14px; color: #228be6;">
            ${gene.label || gene.name || 'Unknown Gene'}
        </strong>
    </div>`;
    
    // Basic information
    html += '<div class="tooltip-section">';
    html += `<div class="tooltip-row">
        <span class="tooltip-label">Position:</span>
        <span class="tooltip-value">${gene.start.toLocaleString()} - ${gene.end.toLocaleString()} (${(gene.end - gene.start).toLocaleString()} bp)</span>
    </div>`;
    html += `<div class="tooltip-row">
        <span class="tooltip-label">Strand:</span>
        <span class="tooltip-value">${gene.strand > 0 ? 'Forward (+)' : 'Reverse (-)'}</span>
    </div>`;
    
    // Metadata from GFF
    const metadata = gene.metadata;
    const importantFields = ['Target_ID', 'Alias', 'SeqID', 'type', 'source', 'score'];
    
    for (const field of importantFields) {
        if (metadata[field] && metadata[field] !== '.' && metadata[field] !== '') {
            html += `<div class="tooltip-row">
                <span class="tooltip-label">${field}:</span>
                <span class="tooltip-value">${metadata[field]}</span>
            </div>`;
        }
    }
    html += '</div>';
    
    // Additional attributes
    const additionalAttrs = Object.keys(metadata).filter(k => 
        !importantFields.includes(k) && 
        !['id'].includes(k) &&
        metadata[k] && 
        metadata[k] !== '.' && 
        metadata[k] !== ''
    );
    
    if (additionalAttrs.length > 0) {
        html += '<div class="tooltip-section">';
        html += '<div style="margin-top: 8px; margin-bottom: 4px; font-weight: bold; color: #868e96;">Additional Attributes:</div>';
        for (const attr of additionalAttrs) {
            html += `<div class="tooltip-row">
                <span class="tooltip-label">${attr}:</span>
                <span class="tooltip-value">${metadata[attr]}</span>
            </div>`;
        }
        html += '</div>';
    }
    
    html += '</div>';
    return html;
}

// Enhanced gene detail panel update
function updateGeneDetailPanel(gene, panelId = 'gene-details-panel') {
    /**
     * Update the gene details panel with comprehensive gene information
     */
    const panel = document.getElementById(panelId);
    if (!panel || !gene || !gene.metadata) {
        return;
    }
    
    let html = '<div class="gene-details-expanded">';
    
    // Header
    html += `<div style="margin-bottom: 16px;">
        <div style="font-size: 18px; font-weight: bold; color: #228be6; margin-bottom: 8px;">
            ${gene.label || gene.name || 'Unknown Gene'}
        </div>
        <div style="font-size: 12px; color: #868e96;">
            Click on genes in the visualization to explore their annotations
        </div>
    </div>`;
    
    // Create tabbed sections
    html += '<div class="detail-tabs">';
    
    // Basic Information Tab
    html += '<div class="detail-section">';
    html += '<div class="section-title">📍 Genomic Location</div>';
    html += '<table class="detail-table">';
    html += `<tr><td class="label">Coordinates:</td><td>${gene.start.toLocaleString()} - ${gene.end.toLocaleString()}</td></tr>`;
    html += `<tr><td class="label">Length:</td><td>${(gene.end - gene.start).toLocaleString()} bp</td></tr>`;
    html += `<tr><td class="label">Strand:</td><td>${gene.strand > 0 ? 'Forward (+)' : 'Reverse (-)'}</td></tr>`;
    if (gene.metadata.phase && gene.metadata.phase !== '.') {
        html += `<tr><td class="label">Phase:</td><td>${gene.metadata.phase}</td></tr>`;
    }
    html += '</table>';
    html += '</div>';
    
    // Annotation Information
    html += '<div class="detail-section">';
    html += '<div class="section-title">🏷️ Annotation</div>';
    html += '<table class="detail-table">';
    if (gene.metadata.type) {
        html += `<tr><td class="label">Feature Type:</td><td><span class="badge">${gene.metadata.type}</span></td></tr>`;
    }
    if (gene.metadata.source) {
        html += `<tr><td class="label">Source:</td><td>${gene.metadata.source}</td></tr>`;
    }
    if (gene.metadata.score && gene.metadata.score !== '.') {
        html += `<tr><td class="label">Score:</td><td>${gene.metadata.score}</td></tr>`;
    }
    if (gene.metadata.Target_ID) {
        html += `<tr><td class="label">Target ID:</td><td><code>${gene.metadata.Target_ID}</code></td></tr>`;
    }
    if (gene.metadata.Alias) {
        html += `<tr><td class="label">Alias:</td><td><code>${gene.metadata.Alias}</code></td></tr>`;
    }
    if (gene.metadata.SeqID) {
        html += `<tr><td class="label">Sequence ID:</td><td><code>${gene.metadata.SeqID}</code></td></tr>`;
    }
    html += '</table>';
    html += '</div>';
    
    // All Attributes
    const allAttrs = Object.keys(gene.metadata).filter(k => 
        !['type', 'source', 'score', 'Target_ID', 'Alias', 'SeqID', 'id', 'phase'].includes(k) &&
        gene.metadata[k] && 
        gene.metadata[k] !== '.'
    );
    
    if (allAttrs.length > 0) {
        html += '<div class="detail-section">';
        html += '<div class="section-title">📋 All Attributes</div>';
        html += '<table class="detail-table">';
        for (const attr of allAttrs.sort()) {
            html += `<tr><td class="label">${attr}:</td><td>${gene.metadata[attr]}</td></tr>`;
        }
        html += '</table>';
        html += '</div>';
    }
    
    html += '</div>'; // Close detail-tabs
    html += '</div>'; // Close gene-details-expanded
    
    panel.innerHTML = html;
    
    // Add some inline styles for the detail panel
    const style = document.createElement('style');
    style.textContent = `
        .gene-details-expanded {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            font-size: 13px;
        }
        .detail-section {
            margin-bottom: 20px;
            background: #f8f9fa;
            padding: 12px;
            border-radius: 6px;
        }
        .section-title {
            font-weight: bold;
            margin-bottom: 10px;
            color: #495057;
            font-size: 14px;
        }
        .detail-table {
            width: 100%;
            border-collapse: collapse;
        }
        .detail-table tr {
            border-bottom: 1px solid #dee2e6;
        }
        .detail-table td {
            padding: 6px 8px;
        }
        .detail-table td.label {
            font-weight: 600;
            color: #495057;
            width: 35%;
            vertical-align: top;
        }
        .detail-table code {
            background: #e9ecef;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Monaco', 'Menlo', monospace;
            font-size: 12px;
        }
        .badge {
            background: #228be6;
            color: white;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 600;
        }
    `;
    if (!document.getElementById('gene-detail-styles')) {
        style.id = 'gene-detail-styles';
        document.head.appendChild(style);
    }
}

// Initialize ClusterMap visualization with enhanced features
function initializeClusterMap(containerId, data, config = {}) {
    console.log("Initializing Enhanced ClusterMap for container:", containerId);
    
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
            },
            shape: {
                onClick: function(event, gene) {
                    // Enhanced click handler
                    console.log("Gene clicked:", gene);
                    updateGeneDetailPanel(gene);
                    
                    // Highlight selected gene
                    d3.selectAll('.genePolygon')
                        .style('stroke-width', '1px')
                        .style('opacity', 0.7);
                    d3.select(event.target)
                        .style('stroke-width', '3px')
                        .style('stroke', '#228be6')
                        .style('opacity', 1);
                }
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
        
        // Add enhanced tooltips
        addEnhancedTooltips(containerId);
        
        console.log("Enhanced ClusterMap initialized successfully");
        return chart;
    } catch (error) {
        console.error("Error initializing ClusterMap:", error);
        container.html(`<div style="color: red; padding: 20px;">Error initializing visualization: ${error.message}</div>`);
        return null;
    }
}

// Add enhanced tooltip functionality
function addEnhancedTooltips(containerId) {
    /**
     * Add interactive tooltips to all genes in the visualization
     */
    const container = d3.select(`#${containerId}`);
    
    // Create tooltip div if it doesn't exist
    let tooltip = d3.select('body').selectAll('.enhanced-gene-tooltip').data([null]);
    tooltip = tooltip.enter()
        .append('div')
        .attr('class', 'enhanced-gene-tooltip')
        .style('position', 'absolute')
        .style('background', 'rgba(255, 255, 255, 0.98)')
        .style('color', '#212529')
        .style('padding', '12px')
        .style('border-radius', '6px')
        .style('font-size', '12px')
        .style('pointer-events', 'none')
        .style('z-index', '10000')
        .style('opacity', 0)
        .style('box-shadow', '0 4px 12px rgba(0,0,0,0.15)')
        .style('border', '1px solid #dee2e6')
        .style('max-width', '400px')
        .merge(tooltip);
    
    // Add tooltip CSS
    const tooltipStyle = document.createElement('style');
    tooltipStyle.textContent = `
        .tooltip-header {
            padding-bottom: 8px;
            margin-bottom: 8px;
            border-bottom: 2px solid #228be6;
        }
        .tooltip-section {
            margin-top: 4px;
        }
        .tooltip-row {
            display: flex;
            padding: 3px 0;
            font-size: 11px;
        }
        .tooltip-label {
            font-weight: 600;
            color: #495057;
            min-width: 100px;
            margin-right: 8px;
        }
        .tooltip-value {
            color: #212529;
            word-break: break-word;
        }
    `;
    if (!document.getElementById('tooltip-styles')) {
        tooltipStyle.id = 'tooltip-styles';
        document.head.appendChild(tooltipStyle);
    }
    
    // Attach tooltip to gene polygons
    container.selectAll('.genePolygon')
        .on('mouseenter', function(event, d) {
            tooltip
                .html(createGeneTooltip(d))
                .style('left', (event.pageX + 15) + 'px')
                .style('top', (event.pageY - 10) + 'px')
                .transition()
                .duration(200)
                .style('opacity', 1);
            
            // Highlight gene on hover
            d3.select(this)
                .style('opacity', 1)
                .style('cursor', 'pointer');
        })
        .on('mousemove', function(event) {
            tooltip
                .style('left', (event.pageX + 15) + 'px')
                .style('top', (event.pageY - 10) + 'px');
        })
        .on('mouseleave', function() {
            tooltip
                .transition()
                .duration(200)
                .style('opacity', 0);
            
            d3.select(this)
                .style('opacity', 0.9);
        });
}

// Update ClusterMap configuration
function updateClusterMapConfig(containerId, config) {
    console.log("Updating ClusterMap config for:", containerId);
    
    const container = d3.select(`#${containerId}`);
    if (container.empty()) {
        console.error("Container not found:", containerId);
        return;
    }

    const chart = container.datum() ? ClusterMap.ClusterMap().config(config) : null;
    
    if (chart && container.datum()) {
        try {
            container.call(chart);
            addEnhancedTooltips(containerId);
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
        const zoom = d3.zoom();
        svg.call(zoom.transform, d3.zoomIdentity.translate(20, 50).scale(1.2));
        console.log("View reset successfully");
    } catch (error) {
        console.error("Error resetting view:", error);
    }
}

// Export gene data to CSV
function exportGenesToCSV(containerId, filename = 'genes.csv') {
    console.log("Exporting genes to CSV from:", containerId);
    
    const svg = d3.select(`#${containerId} svg.clusterMap`);
    if (svg.empty()) {
        console.error("No SVG found in container:", containerId);
        return;
    }
    
    try {
        const genes = [];
        svg.selectAll('.gene').each(function(d) {
            if (d && d.metadata) {
                genes.push({
                    name: d.label || d.name,
                    start: d.start,
                    end: d.end,
                    strand: d.strand > 0 ? '+' : '-',
                    length: d.end - d.start,
                    ...d.metadata
                });
            }
        });
        
        if (genes.length === 0) {
            console.warn("No genes found to export");
            return;
        }
        
        // Convert to CSV
        const headers = Object.keys(genes[0]);
        let csv = headers.join(',') + '\n';
        genes.forEach(gene => {
            csv += headers.map(h => {
                const val = gene[h];
                return typeof val === 'string' && val.includes(',') ? `"${val}"` : val;
            }).join(',') + '\n';
        });
        
        // Download
        const blob = new Blob([csv], { type: 'text/csv' });
        downloadFile(blob, filename);
        console.log(`Exported ${genes.length} genes to CSV`);
    } catch (error) {
        console.error("Error exporting to CSV:", error);
    }
}

// Make functions available globally for Dash callbacks
window.syntenyViz = {
    initialize: initializeClusterMap,
    updateConfig: updateClusterMapConfig,
    exportSVG: exportClusterMapSVG,
    exportCSV: exportGenesToCSV,
    resetView: resetClusterMapView,
    serialiseSVG: serialiseSVG,
    downloadFile: downloadFile,
    createGeneTooltip: createGeneTooltip,
    updateGeneDetailPanel: updateGeneDetailPanel,
    addEnhancedTooltips: addEnhancedTooltips
};
