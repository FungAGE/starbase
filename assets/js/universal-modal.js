// Universal Modal System for Starbase
// Provides consistent modal styling and behavior across all pages

class UniversalModal {
    constructor() {
        this.modal = null;
        this.isInitialized = false;
    }

    init() {
        if (this.isInitialized) return;

        // Create modal HTML
        const modalHTML = `
            <div id="universal-modal" style="
                display: none;
                position: fixed;
                z-index: 10000;
                left: 0;
                top: 0;
                width: 100%;
                height: 100%;
                background-color: rgba(0,0,0,0.5);
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
            ">
                <div style="
                    background-color: white;
                    margin: 2% auto;
                    padding: 0;
                    border-radius: 8px;
                    width: 90%;
                    max-width: 800px;
                    max-height: 90vh;
                    overflow-y: auto;
                    box-shadow: 0 10px 25px rgba(0,0,0,0.2);
                    position: relative;
                ">
                    <div style="
                        position: sticky;
                        top: 0;
                        background: white;
                        padding: 20px;
                        margin: -20px -20px 20px -20px;
                        border-bottom: 1px solid #dee2e6;
                        z-index: 1;
                        border-radius: 8px 8px 0 0;
                    ">
                        <span id="universal-modal-close" style="
                            color: #868e96;
                            float: right;
                            font-size: 28px;
                            font-weight: bold;
                            cursor: pointer;
                            transition: color 0.2s;
                            padding: 10px;
                            margin: -10px -10px -10px 10px;
                        ">&times;</span>
                        <h2 id="universal-modal-title" style="
                            margin: 0 0 8px 0;
                            color: #1a1a1a;
                            font-size: clamp(20px, 5vw, 28px);
                            font-weight: 700;
                            line-height: 1.2;
                            text-align: center;
                            padding: 16px 20px;
                            width: 100%;
                            box-sizing: border-box;
                        "></h2>
                    </div>
                    <div id="universal-modal-content" style="
                        padding: 0;
                        color: #212529;
                        line-height: 1.6;
                        width: 100%;
                    "></div>
                </div>
            </div>
        `;

        document.body.insertAdjacentHTML('beforeend', modalHTML);
        this.modal = document.getElementById('universal-modal');

        // Add event listeners
        document.getElementById('universal-modal-close').onclick = () => this.close();

        this.modal.onclick = (e) => {
            if (e.target === this.modal) {
                this.close();
            }
        };

        // Add CSS styles
        const style = document.createElement('style');
        style.textContent = `
            #universal-modal-close:hover {
                color: #000;
            }

            #universal-modal-content {
                max-width: 900px;
                margin: 0 auto;
            }

            #universal-modal-title {
                box-sizing: border-box;
                padding: 0 20px;
            }

            .modal-container {
                display: flex;
                flex-direction: column;
                gap: 20px;
                align-items: center;
                width: 100%;
                box-sizing: border-box;
                padding: 0 20px;
            }

            .modal-section {
                margin-bottom: 24px;
                background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
                border: 1px solid #e9ecef;
                border-radius: 12px;
                padding: 24px;
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
                transition: box-shadow 0.2s ease;
                width: 100%;
                box-sizing: border-box;
            }

            .modal-section:hover {
                box-shadow: 0 4px 16px rgba(0, 0, 0, 0.08);
            }

            .section-title {
                font-weight: 700;
                font-size: 20px;
                margin-bottom: 20px;
                color: #1a1a1a;
                text-align: center;
                position: relative;
                padding-bottom: 12px;
            }

            .section-title::after {
                content: '';
                position: absolute;
                bottom: 0;
                left: 50%;
                transform: translateX(-50%);
                width: 60px;
                height: 3px;
            }

            .section-content {
                color: #495057;
            }

            .modal-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
                gap: 20px;
                margin-bottom: 0;
            }

            .modal-row {
                display: flex;
                justify-content: center;
                align-items: center;
                gap: 16px;
                padding: 12px 16px;
                background: rgba(255, 255, 255, 0.7);
                border-radius: 8px;
                border: 1px solid #f1f3f5;
                transition: all 0.2s ease;
                min-height: 48px;
            }

            .modal-row:hover {
                background: rgba(255, 255, 255, 0.9);
                border-color: #dee2e6;
                transform: translateY(-1px);
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
            }

            .modal-label {
                font-weight: 600;
                color: #495057;
                font-size: 14px;
                letter-spacing: 0.5px;
                flex-shrink: 0;
                min-width: 120px;
            }

            .modal-value {
                color: #212529;
                font-weight: 500;
                text-align: center;
                flex: 1;
                word-break: break-word;
            }

            .modal-badge {
                display: inline-block;
                padding: 4px 8px;
                border-radius: 4px;
                font-size: 12px;
                font-weight: 500;
            }

            .badge-green {
                background-color: #e9faf0;
                color: #099268;
            }

            .badge-yellow {
                background-color: #fff8e1;
                color: #f59f00;
            }

            .badge-blue {
                background-color: #e7f5ff;
                color: #1c7ed6;
            }

            .badge-red {
                background-color: #fef2f2;
                color: #dc2626;
            }

            .modal-link {
                color: #228be6;
                text-decoration: none;
                font-weight: 500;
            }

            .modal-link:hover {
                text-decoration: underline;
            }

            .quality-tags {
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                justify-content: center;
                margin-top: 4px;
            }

            .alert {
                padding: 16px 20px;
                border-radius: 8px;
                margin-bottom: 20px;
                border: 1px solid;
                text-align: center;
                background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.04);
            }

            .alert-red {
                background: linear-gradient(135deg, #fef2f2 0%, #fef2f2 100%);
                color: #dc2626;
                border-color: #fecaca;
            }

            .alert-yellow {
                background: linear-gradient(135deg, #fffbeb 0%, #fffbeb 100%);
                color: #d97706;
                border-color: #fed7aa;
            }

            .alert-blue {
                background: linear-gradient(135deg, #f0f9ff 0%, #f0f9ff 100%);
                color: #0284c7;
                border-color: #bae6fd;
            }
        `;
        document.head.appendChild(style);

        this.isInitialized = true;
    }

    show(title, content) {
        if (!this.isInitialized) {
            this.init();
        }

        document.getElementById('universal-modal-title').innerHTML = title;
        document.getElementById('universal-modal-content').innerHTML = content;
        this.modal.style.display = 'block';

        // Prevent body scroll when modal is open
        document.body.style.overflow = 'hidden';
    }

    close() {
        this.modal.style.display = 'none';
        document.body.style.overflow = 'auto';
    }

    // Static method to render modal content from structured data
    static renderAccessionModal(data) {
        let html = '';

        // Handle error case
        if (data.error) {
            return `
                <div class="modal-container">
                    <div class="alert alert-red">
                        <div style="font-weight: 600; margin-bottom: 8px;">Error</div>
                        <div>${data.error}</div>
                    </div>
                </div>
            `;
        }

        // Wrapper for better centering
        html += '<div class="modal-container">';

        // Basic ship information section
        if (data.familyName || data.genomes_present) {
            html += `
                <div class="modal-section">
                    <div class="section-title">🚀 Starship Information</div>
                    <div class="section-content">
                        <div class="modal-grid">
                            ${data.familyName ? `
                                <div class="modal-row">
                                    <span class="modal-label">Starship Family:</span>
                                    <span class="modal-value">${data.familyName}</span>
                                </div>
                            ` : ''}
                            ${data.navis_name && ['Aristaeus', 'Arwing', 'Atlantia', 'Chrysoar', 'Defiant', 'Enterprise', 'Galactica', 'Hephaestus', 'Horizon', 'Mithridate', 'Osiris', 'Phoenix', 'Prometheus', 'Typhon', 'Voyager', 'Wallaby'].includes(data.navis_name) ? `
                                <div class="modal-row">
                                    <span class="modal-label">Starship Navis:</span>
                                    <span class="modal-value">${data.navis_name}</span>
                                </div>
                            ` : ''}
                            ${data.genomes_present ? `
                                <div class="modal-row">
                                    <span class="modal-label">Genomes Present:</span>
                                    <span class="modal-badge badge-blue">${data.genomes_present}</span>
                                </div>
                            ` : ''}

                            </div>
                    </div>
                </div>
            `;
        }

        // Curation status and quality tags
        if (data.curated_status || (data.quality_tags && data.quality_tags.length > 0)) {
            html += `
                <div class="modal-section">
                    <div class="section-title">Curation Status</div>
                    <div class="section-content">
                        <div class="modal-grid">
                            ${data.curated_status ? `
                                <div class="modal-row">
                                    <span class="modal-label">Status:</span>
                                    <span class="modal-badge badge-${data.curated_status === 'curated' ? 'green' : 'yellow'}">${data.curated_status}</span>
                                </div>
                            ` : ''}
                            ${data.quality_tags && data.quality_tags.length > 0 ? `
                                <div class="modal-row">
                                    <span class="modal-label">Quality:</span>
                                    <div class="quality-tags">
                                        ${data.quality_tags.map(tag => `<span class="modal-badge badge-red">${tag}</span>`).join('')}
                                    </div>
                                </div>
                            ` : ''}
                        </div>
                    </div>
                </div>
            `;
        }

        // Taxonomy section
        if (data.order || data.family || data.species_name) {
            html += `
                <div class="modal-section">
                    <div class="section-title">Taxonomy</div>
                    <div class="section-content">
                        <div class="modal-grid">
                            ${data.order ? `
                                <div class="modal-row">
                                    <span class="modal-label">Order:</span>
                                    <span class="modal-value">${data.order}</span>
                                </div>
                            ` : ''}
                            ${data.family ? `
                                <div class="modal-row">
                                    <span class="modal-label">Family:</span>
                                    <span class="modal-value">${data.family}</span>
                                </div>
                            ` : ''}
                            ${data.species_name ? `
                                <div class="modal-row">
                                    <span class="modal-label">Species:</span>
                                    <span class="modal-value">${data.species_name}</span>
                                </div>
                            ` : ''}
                            ${data.tax_id ? `
                                <div class="modal-row">
                                    <span class="modal-label">NCBI Taxonomy ID:</span>
                                    <span class="modal-value">
                                        ${data.tax_id !== 'N/A' ?
                                            `<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${data.tax_id}" target="_blank" class="modal-link">${data.tax_id}</a>` :
                                            'N/A'
                                        }
                                    </span>
                                </div>
                            ` : ''}
                        </div>
                    </div>
                </div>
            `;
        }

        // Genome details section
        if (data.assembly_accession || data.genome_source || data.contig_id || data.element_position || data.element_length) {
            html += `
                <div class="modal-section">
                    <div class="section-title">Genome Details</div>
                    <div class="section-content">
                        <div class="modal-grid">
                            ${data.assembly_accession ? `
                                <div class="modal-row">
                                    <span class="modal-label">Assembly Accession:</span>
                                    <span class="modal-value">${data.assembly_accession}</span>
                                </div>
                            ` : ''}
                            ${data.genome_source ? `
                                <div class="modal-row">
                                    <span class="modal-label">Genome Source:</span>
                                    <span class="modal-value">${data.genome_source}</span>
                                </div>
                            ` : ''}
                            ${data.contig_id ? `
                                <div class="modal-row">
                                    <span class="modal-label">Contig ID:</span>
                                    <span class="modal-value">${data.contig_id}</span>
                                </div>
                            ` : ''}
                            ${data.element_position ? `
                                <div class="modal-row">
                                    <span class="modal-label">Element Position:</span>
                                    <span class="modal-value">${data.element_position}</span>
                                </div>
                            ` : ''}
                            ${data.element_length ? `
                                <div class="modal-row">
                                    <span class="modal-label">Size:</span>
                                    <span class="modal-value">${data.element_length} bp</span>
                                </div>
                            ` : ''}
                        </div>
                    </div>
                </div>
            `;
        }

        // Close the modal container wrapper
        html += '</div>';

        return html;
    }
}

// Global instance - only create if not already exists
if (typeof window.universalModal === 'undefined') {
    window.universalModal = new UniversalModal();
}

// Global function for easy access - only create if not already exists
if (typeof window.showUniversalModal === 'undefined') {
    window.showUniversalModal = function(title, content) {
        window.universalModal.show(title, content);
    };
}

// Function to show accession modal with data fetching - only create if not already exists
if (typeof window.showAccessionModal === 'undefined') {
    window.showAccessionModal = async function(accessionId) {
        try {
            const response = await fetch(`/api/accession/accession_details/${accessionId}`);
            const data = await response.json();

            if (data.error) {
                window.showUniversalModal('Error', UniversalModal.renderAccessionModal(data));
                return;
            }

            const title = data.title || `Starship Accession: ${accessionId}`;
            const content = UniversalModal.renderAccessionModal(data);
            window.showUniversalModal(title, content);
        } catch (error) {
            console.error('Error fetching accession details:', error);
            window.showUniversalModal('Error', UniversalModal.renderAccessionModal({
                error: `Error loading details for ${accessionId}. Details: ${error.message}`
            }));
        }
    };
}
