if (!window.dash_clientside) {
    window.dash_clientside = {};
}

window.dash_clientside.handleCellClick = function(e) {
    if (!e) return null;
    if (e.colDef.field === 'accession_tag') {
        return JSON.stringify({
            accession: e.value,
            column: e.colDef.field,
            data: e.data
        });
    }
    return null;
}