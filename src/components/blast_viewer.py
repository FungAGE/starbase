from dash import html
import dash_mantine_components as dmc

class BlastViewer(html.Div):
    def __init__(self, blast_results):
        super().__init__([
            dmc.Container([
                html.Div(id="blast-multiple-alignments"),
                html.Div(id="blast-alignments-table"),
                html.Div(id="blast-single-alignment"),
                # Include the JS directly in the component
                html.Script(src="/assets/js/blaster.min.js"),
                html.Script(src="/assets/js/html2canvas.min.js"),
                # Initialize blasterjs with the results
                html.Script(f"""
                    var alignments = `{blast_results}`;
                    var blasterjs = require("biojs-vis-blasterjs");
                    var instance = new blasterjs({{
                        string: alignments,
                        multipleAlignments: "blast-multiple-alignments",
                        alignmentsTable: "blast-alignments-table",
                        singleAlignment: "blast-single-alignment",
                    }});
                """)
            ])
        ])