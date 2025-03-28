
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        initBlaster: function(data) {
            if (data && data.blast_text) {
                try {
                    // Clear existing content
                    document.getElementById('blast-multiple-alignments').innerHTML = '';
                    document.getElementById('blast-alignments-table').innerHTML = '';
                    // document.getElementById('blast-single-alignment').innerHTML = '';
                    
                    // Initialize BlasterJS
                    var blasterjs = require("biojs-vis-blasterjs");
                    var instance = new blasterjs({
                        string: data.blast_text,
                        multipleAlignments: "blast-multiple-alignments",
                        alignmentsTable: "blast-alignments-table"
                        // singleAlignment: "blast-single-alignment"
                    });
                    
                    return window.dash_clientside.no_update;
                } catch (error) {
                    console.error('Error initializing BlasterJS:', error);
                    return error.toString();
                }
            }
            return window.dash_clientside.no_update;
        }
    }
});
