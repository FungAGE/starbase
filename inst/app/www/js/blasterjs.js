require=(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
(function (global){
/** @preserve http://github.com/easeway/js-class */

// Class Definition using ECMA5 prototype chain

function inherit(dest, src, noParent) {
    while (src && src !== Object.prototype) {
        Object.getOwnPropertyNames(src).forEach(function (name) {
            if (name != '.class' && !dest.hasOwnProperty(name)) {
                var desc = Object.getOwnPropertyDescriptor(src, name);
                Object.defineProperty(dest, name, desc);
            }
        });
        if (noParent) {
            break;
        }
        src = src.__proto__;
    }
    return dest;
}

var Class = function (base, proto, options) {
    if (typeof(base) != 'function') {
        options = proto;
        proto = base;
        base = Object;
    }
    if (!proto) {
        proto = {};
    }
    if (!options) {
        options = {};
    }
    
    var meta = {
        name: options.name,
        base: base,
        implements: []
    }
    var classProto = Class.clone(proto);
    if (options.implements) {
        (Array.isArray(options.implements) ? options.implements : [options.implements])
            .forEach(function (implementedType) {
                if (typeof(implementedType) == 'function' && implementedType.prototype) {
                    meta.implements.push(implementedType);
                    Class.extend(classProto, implementedType.prototype);
                }
            });
    }
    classProto.__proto__ = base.prototype;
    var theClass = function () {
        if (typeof(this.constructor) == 'function') {
            this.constructor.apply(this, arguments);
        }
    };
    meta.type = theClass;
    theClass.prototype = classProto;
    Object.defineProperty(theClass, '.class.meta', { value: meta, enumerable: false, configurable: false, writable: false });
    Object.defineProperty(classProto, '.class', { value: theClass, enumerable: false, configurable: false, writable: false });
    if (options.statics) {
        Class.extend(theClass, options.statics);
    }
    return theClass;
};

Class.extend = inherit;

Class.clone = function (object) {
    return inherit({}, object);
};

function findType(meta, type) {
    while (meta) {
        if (meta.type.prototype === type.prototype) {
            return true;
        }
        for (var i in meta.implements) {
            var implType = meta.implements[i];
            var implMeta = implType['.class.meta'];
            if (implMeta) {
                if (findType(implMeta, type)) {
                    return true;
                }
            } else {
                for (var proto = implType.prototype; proto; proto = proto.__proto__) {
                    if (proto === type.prototype) {
                        return true;
                    }
                }
            }
        }
        meta = meta.base ? meta.base['.class.meta'] : undefined;
    }
    return false;
}

var Checker = Class({
    constructor: function (object) {
        this.object = object;
    },
    
    typeOf: function (type) {
        if (this.object instanceof type) {
            return true;
        }
        var meta = Class.typeInfo(this.object);
        return meta && findType(meta, type);
    }
});

// aliases
Checker.prototype.a = Checker.prototype.typeOf;
Checker.prototype.an = Checker.prototype.typeOf;

Class.is = function (object) {
    return new Checker(object);
};

Class.typeInfo = function (object) {
    var theClass = object.__proto__['.class'];
    return theClass ? theClass['.class.meta'] : undefined;
};

Class.VERSION = [0, 0, 2];

if (module) {
    module.exports = Class;
} else {
    global.Class = Class;   // for browser
}
}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],2:[function(require,module,exports){
(function(factory) {
  if(typeof exports === 'object') {
    factory(exports);
  } else {
    factory(this);
  }
}).call(this, function(root) { 

  var slice   = Array.prototype.slice,
      each    = Array.prototype.forEach;

  var extend = function(obj) {
    if(typeof obj !== 'object') throw obj + ' is not an object' ;

    var sources = slice.call(arguments, 1); 

    each.call(sources, function(source) {
      if(source) {
        for(var prop in source) {
          if(typeof source[prop] === 'object' && obj[prop]) {
            extend.call(obj, obj[prop], source[prop]);
          } else {
            obj[prop] = source[prop];
          }
        } 
      }
    });

    return obj;
  }

  root.extend = extend;
});

},{}],"biojs-vis-blasterjs":[function(require,module,exports){
/*
 * biojs-vis-blasterjs
 * https://github.com/sing-group/blasterjs
 *
 * Copyright (c) 2016 SING - Sistemas Informaticos de Nueva Generacion
 * Licensed under the MIT license.
 *
 *
 * BlasterJS
 *
 * @class
 * @extends Biojs
 *
 * @author <a href="mailto:aiblanco@uvigo.es">Aitor Blanco Miguez</a>
 * @version 0.1.1
 * @category 3
 *
 * @requires Bootstrap 3
 * @dependency <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous"></link>
 *
 * @requires <a href='html2canvas.js'>HTML2CANVAS</a>
 * @dependency <script type="text/javascript" src="lib/html2canvas.js"></script>
 *
 *
 * @param {Object} options An object with the options for BlasterJS component.
 *
 * @option {string} input
 *    Identifier of the INPUT tag where the BLAST output file should be selected.
 *
 * @option {string} multipleAlignments
 *    Identifier of the DIV tag where the multiple alignments should be displayed.
 *
 * @option {string} alignmentsTable
 *    Identifier of the DIV tag where the alignments table should be displayed.
 *
 * @option {string} singleAlignment
 *    Identifier of the DIV tag where the single should be displayed.
 * 
 * @example
 * var instance = new Biojs.blasterjs({
 *      input: "blastinput",
 *      multipleAlignments: "blast-multiple-alignments",
 *      alignmentsTable: "blast-alignments-table",
        singleAlignment: "blast-single-alignment"
 * });
 */


var blasterjs;
var Class = require('js-class');
module.exports = blasterjs = Class(
    /** @lends Biojs.blasterjs# */
    {
        constructor: function (options) {
		var self = this;

        require('js-extend').extend(this.opt, options);
            
        var css = 'table tbody tr td button:hover{ text-decoration: underline;}';
        var style = document.createElement('style');
        if (style.styleSheet) {
            style.styleSheet.cssText = css;
        } else {
            style.appendChild(document.createTextNode(css));
        }
        document.getElementsByTagName('head')[0].appendChild(style);
        document.getElementById(self.opt.input).addEventListener('change',  function(evt) { self._displayAlignments(evt, self) }, false);
            
        },
        
        /**
         * Default values for the options
         * @name Biojs.blasterjs-opt
         */
        opt : {
             input: "blastinput",
             multipleAlignments: "blast-multiple-alignments",
             alignmentsTable: "blast-alignments-table",
             singleAlignment: "blast-single-alignment"
        },
        
        /**
         * Private: Read and display BLAST alignments.
         * @ignore
         */
        _displayAlignments : function(evt, self) {
            if (window.File && window.FileReader && window.FileList && window.Blob) {
                var f = evt.target.files[0];
                if (f) {
                    var r = new FileReader();
                    r.onload = function(e) { 
                        var queryLenght = getQueryLenght(e.target.result); 
                        var alignments  = getAlignments(e.target.result);
                        createAlignmentComparison(alignments, queryLenght, true, self);
                        createAlignmentTable(alignments, self);
                        createSingleAlignment(alignments[0], self);
                    }
                    r.readAsText(f);
                } else { 
                    alert('Failed to load file');
                }
            } else {
              alert('The File APIs are not fully supported by your browser.');
            }  
        
            function BlastAlignment (description,
                                    length,
                                    score,
                                    eValue,
                                    identities,
                                    positives,
                                    gaps,
                                    queryStart,
                                    query,
                                    queryEnd,
                                    comparison,
                                    subjectStart,
                                    subject,
                                    subjectEnd){                
                this.description  = description;
                this.length       = length;
                this.score        = score;
                this.eValue       = eValue;
                this.identities   = identities;
                this.positives    = positives;
                this.gaps         = gaps;
                this.queryStart   = queryStart;
                this.query        = query;
                this.queryEnd     = queryEnd;
                this.comparison   = comparison;
                this.subjectStart = subjectStart;
                this.subject      = subject;
                this.subjectEnd   = subjectEnd;
            }

            String.prototype.startsWith = function(prefix) {
                return this.indexOf(prefix) === 0;
            }

            function getQueryLenght(content){
                var lines = content.split('\n');
                var length = 0;
                for (var i = 0; i < lines.length; i++){
                    if(lines[i].startsWith('Length=')){
                        length = lines[i].split('=')[1];
                        break;
                    }
                }
                return length;
            }

            function getAlignments(content){
                var lines = content.split('\n');
                return parseBlastText(content);  
            }

            function parseBlastText(content){
                var lines = content.split('\n');
                var alignments = [];
                for (var i = 0; i < lines.length; i++){
                    if(lines[i].startsWith('>')){
                        var line1 = lines[i].split('>')[1];
                        var line2 = "";
                        var currentLine = i;
                        while (line2 == ""){
                            currentLine = currentLine +1;
                            if(lines[currentLine].startsWith('Length=')){
                                line2 = lines[currentLine];
                            }else{
                                line1 += lines[currentLine];
                            }
                        }
                        var description  = line1;
                        var length       = line2.split('=')[1];

                        if(lines[currentLine+2].startsWith(' Features in this part of subject sequence:')){                
                            currentLine = currentLine+3;
                        }
                        var score        = lines[currentLine+2].split(',')[0].split(' ')[3];
                        var eValue       = lines[currentLine+2].split(',')[1].split(' ')[4];
                        var identities   = lines[currentLine+3].split(',')[0].split('(')[1].substr(0, lines[currentLine+3].split(',')[0].split('(')[1].length-2);
                        if(lines[0].startsWith('BLASTN')){
                            var positives = 'N/A';
                            var gaps    = lines[currentLine+3].split(',')[1].split('(')[1].substr(0, lines[currentLine+3].split(',')[1].split('(')[1].length-2);
                        }else{
                            var positives    = lines[currentLine+3].split(',')[1].split('(')[1].substr(0, lines[currentLine+3].split(',')[1].split('(')[1].length-2);
                            var gaps         = lines[currentLine+3].split(',')[2].split('(')[1].substr(0, lines[currentLine+3].split(',')[2].split('(')[1].length-2);
                        }
                        if(lines[currentLine+4].split(',')[0].split(' ')[1] == 'Frame' || lines[currentLine+4].startsWith(' Strand')){
                            currentLine = currentLine+1;   
                        }
                        var queryStart = lines[currentLine+5].substring(5).replace(/^\s+/g, '').split(' ')[0];
                        var query = lines[currentLine+5].substring(5).replace(/\s+/g, '').replace(/[0-9]/g, '');
                        var queryEnd = lines[currentLine+5].substring(5).replace(/^\s+/g, '').split(' ')[lines[currentLine+5].substring(5).replace(/^\s+/g, '').split(' ').length-1];
                        var comparison = lines[currentLine+6].replace(/^\s+/g, ''); 
                        var sbjctStart = lines[currentLine+7].substring(5).replace(/^\s+/g, '').split(' ')[0];
                        var sbjct = lines[currentLine+7].substring(5).replace(/\s+/g, '').replace(/[0-9]/g, '');
                        var sbjctEnd = lines[currentLine+7].substring(5).replace(/^\s+/g, '').split(' ')[lines[currentLine+7].substring(5).replace(/^\s+/g, '').split(' ').length-1];

                        currentLine = currentLine+9;
                        while (lines[currentLine].startsWith('Query')){
                            var nextQuery = lines[currentLine].substring(5).replace(/\s+/g, '').replace(/[0-9]/g, '');
                            query += nextQuery;
                            queryEnd = lines[currentLine].substring(5).replace(/^\s+/g, '').split(' ')[lines[currentLine].substring(5).replace(/^\s+/g, '').split(' ').length-1];
                            sbjct += lines[currentLine+2].substring(5).replace(/\s+/g, '').replace(/[0-9]/g, '');
                            sbjctEnd = lines[currentLine+2].substring(5).replace(/^\s+/g, '').split(' ')[lines[currentLine+2].substring(5).replace(/^\s+/g, '').split(' ').length-1];

                            var nextComparison = lines[currentLine+1].replace(/^\s+/g, '');
                            if(nextQuery.length > nextComparison.length){
                                var diference = nextQuery.length-nextComparison.length;
                                for(var j = 0; j < diference; j++){
                                    nextComparison = ' '+nextComparison;   
                                }
                            }
                            comparison += nextComparison;
                            currentLine = currentLine+4;
                        }   

                        var alignment = new BlastAlignment( description, length, score, eValue, identities, positives, gaps, queryStart, query, queryEnd, comparison, sbjctStart, sbjct, sbjctEnd );
                        alignments.push(alignment);
                    }
                }
                return alignments;
            }

            function getColor(colored, score){
                if(score<40){
                    return getDivColor(colored, 1);
                }else{
                    if(score>=40 || score<50){
                        return getDivColor(colored, 2);
                    }else{
                        if(score>=50 || score<80){
                            return getDivColor(colored, 3);
                        }else{
                            if(score>=80 || score<200){
                                return getDivColor(colored, 4);
                            }else{
                                return getDivColor(colored, 5);
                            }
                        }
                    }
                }
            }

            function getDivColorText(div){
                switch(div) {
                    case 1:
                        return '<40';
                        break;
                    case 2:
                        return '40-50';
                        break;
                    case 3:
                        return '50-80';
                        break;
                    case 4:
                        return '80-200';
                        break;
                    case 5:
                        return '>=200';
                        break;
                    default:
                        return '0';
                }
            }

            function getDivColor(colored, div){
                if(colored){
                    switch(div) {
                        case 1:
                            return '#5C6D7E';
                            break;
                        case 2:
                            return '#9B59B6';
                            break;
                        case 3:
                            return '#5CACE2';
                            break;
                        case 4:
                            return '#57D68D';
                            break;
                        case 5:
                            return '#C0392B';
                            break;
                        default:
                            return '#FFF';
                    }
                }else{
                    switch(div) {
                        case 1:
                            return '#BCBCBC';
                            break;
                        case 2:
                            return '#989898';
                            break;
                        case 3:
                            return '#747474';
                            break;
                        case 4:
                            return '#565656';
                            break;
                        case 5:
                            return '#343434';
                            break;
                        default:
                            return '#FFF';
                    }
                }
            }

            function createAlignmentDiv(color, width1, width2, alignment){ 
                var container = document.createElement('div');
                var div1      = document.createElement('div');
                var div2      = document.createElement('div');
                var divClear  = document.createElement('div');
                var a         = document.createElement('a');
                container.style.minHeight = '12px';
                div1.style.width           = width1+'px';
                div1.style.minHeight       = '4px';
                div1.style.float           = 'left';
                div2.style.width           = width2+'px';
                div2.style.minHeight       = '12px';
                div2.style.float           = 'left';
                div2.style.backgroundColor =  color;    
                div2.onmouseout            = function(){document.getElementById('defline').value=' Mouse over to show defline and scores, click to show alignments';};
                div2.onmouseover           = function(){document.getElementById('defline').value=' '+alignment.description+'. S='+alignment.score+' E='+alignment.eValue;};
                divClear.style.clear = 'both';
                a.href = '#'+alignment.description.split(' ')[0];  
                a.onclick=function(){ selectAlignment(alignment.description.split(' ')[0]); };
                a.appendChild(div2);
                container.appendChild(div1);
                container.appendChild(a);
                container.appendChild(divClear);
                return container;
            }

            function selectAlignment(alignment){
                var item = document.getElementById(alignment).parentElement.parentElement;
                var items = document.getElementsByClassName('alignment-table-description');
                var i;
                for (i = 0; i < items.length; i++) {
                    items[i].style.fontWeight = 'normal';
                    items[i].parentElement.parentElement.style.fontWeight = 'normal';
                }
                item.style.fontWeight = 'bold';
                document.getElementById(alignment).style.fontWeight = 'bold';
            }

            function createColorsDiv(colored){
                var container  = document.createElement('div');
                var divSpace   = document.createElement('div');
                var divClear   = document.createElement('div');
                container.style.color = '#EEE';
                divSpace.style.minWidth  = '50px';
                divSpace.style.minHeight = '10px';
                divSpace.style.float     = 'left';
                container.appendChild(divSpace);
                for(var i = 1; i <= 5; i++){
                    var div = document.createElement('div');
                    div.style.minWidth        = '100px';
                    div.style.textAlign       = 'center';
                    div.style.float           = 'left';        
                    div.innerHTML             = getDivColorText(i).bold();
                    div.style.backgroundColor = getDivColor(colored, i);
                    container.appendChild(div);
                }
                divClear.style.clear = 'both';
                container.appendChild(divClear);
                return container;
            }

            function createQueryDiv(colored){
                var container  = document.createElement('div');
                var divSpace   = document.createElement('div');
                var divQuery   = document.createElement('div');
                var divClear   = document.createElement('div');
                container.style.marginTop = '3px';
                container.style.color     = '#5C6D7E';
                container.style.fontSize  = '10px';
                divSpace.style.width = '50px';
                divSpace.innerHTML   = 'QUERY'.bold();
                divSpace.style.float = 'left';
                divQuery.style.width     = '500px';
                divQuery.style.height    = '10px';
                divQuery.style.float     = 'left';
                divQuery.style.marginTop = '2px';
                divClear.style.clear = 'both';
                if(colored){
                    divQuery.style.backgroundColor = '#C0392B'; 
                } else{
                    divQuery.style.backgroundColor = '#343434';         
                }
                container.appendChild(divSpace);
                container.appendChild(divQuery);
                container.appendChild(divClear);
                return container;
            }

            function createNumbersDiv(lenght){
                var container    = document.createElement('div');
                var divSpace     = document.createElement('div');
                var divNumbers   = document.createElement('div');
                var divClear     = document.createElement('div');
                container.style.marginBottom = '5px';
                container.style.fontSize     = '11px';
                divSpace.style.minWidth  = '50px';
                divSpace.style.minHeight = '10px';
                divSpace.style.float     = 'left';    
                divNumbers.style.float = 'left';
                divClear.style.clear = 'both';  
                divNumbers = divideDivNumbers(divNumbers, lenght);  
                container.appendChild(divSpace);
                container.appendChild(divNumbers);
                container.appendChild(divClear);
                return container;
            }

            function divideDivNumbers(container, lenght){
                var divClear = document.createElement('div');
                if(lenght > 4){
                    if(lenght % 5 == 0){
                        container = createDivisionsDivNumbers(container, 5, lenght/5, 100);
                    }else{
                        var pixels = 500/(5+((lenght%5)/5));
                        container = createDivisionsDivNumbers(container, 5, parseInt(lenght/5), parseInt(pixels));            
                        //Podemos hacerlo o no
                        var pxrest = parseInt(500-(pixels*5));
                        var div = document.createElement('div');
                        div.style.float = 'left';
                        div.style.width = pxrest+'px';
                        div.style.textAlign = 'right';
                        div.innerHTML = lenght;
                        container.appendChild(div);
                    }
                }else{
                    container = createDivisionsDivNumbers(container, lenght, 1, parseInt(500/lenght));
                }    
                divClear.style.clear = 'both'; 
                container.appendChild(divClear);
                return container;
            }

            function createDivisionsDivNumbers(container, divisions, size, pixels){
                for(var i = 0; i<divisions; i++){
                    if(i == 0){
                        var px2  = pixels/2;
                        var div1 = document.createElement('div');
                        var div2 = document.createElement('div');
                        div1.style.float     = 'left';
                        div1.style.width     = px2+'px';
                        div1.style.textAlign = 'left';
                        div1.innerHTML       = '0';
                        div2.style.float     = 'left';
                        div2.style.width     = px2+'px';
                        div2.style.textAlign = 'right';
                        div2.innerHTML       = size*(i+1);
                        container.appendChild(div1); 
                        container.appendChild(div2); 
                    }else{
                        var div = document.createElement('div');
                        div.style.float     = 'left';
                        div.style.width     = pixels+'px';
                        div.style.textAlign = 'right';
                        div.innerHTML       = size*(i+1);
                        container.appendChild(div);
                    }
                }
                return container;
            }


            function createHeader(container, colored, lenght){   
                var text    = document.createElement('div');
                var colors  = createColorsDiv(colored);
                var query   = createQueryDiv(colored);
                var numbers = createNumbersDiv(lenght);
                text.style.color         = '#5C6D7E';
                text.style.textAlign     = 'center';
                text.style.paddingBottom = '5px';
                text.innerHTML           = 'COLOR KEY FOR ALIGNMENT SCORES'.bold();
                container.appendChild(text);
                container.appendChild(colors);
                container.appendChild(query);
                container.appendChild(numbers);
                return container;
            }

            function createBody(alignments, queryLenght, container, colored){
                var alignmentContainer = document.createElement('div');
                alignmentContainer.style.paddingBottom = '10px';
                for(var i = 0; i < alignments.length; i++){
                    if(parseInt(alignments[i].queryStart)>parseInt(alignments[i].queryEnd)){
                        var queryStart = alignments[i].queryEnd;
                        var queryEnd   = alignments[i].queryStart;
                    }else{
                        var queryStart = alignments[i].queryStart;
                        var queryEnd   = alignments[i].queryEnd;             
                    }
                    var white     = parseInt(50+((500*(queryStart-1))/queryLenght));
                    var black     = parseInt(550-white-(500*(queryLenght-queryEnd)/queryLenght));
                    var alignment = createAlignmentDiv(getColor(colored, alignments[i].score), white, black, alignments[i]);
                    alignment.style.marginBottom = '4px';
                    alignmentContainer.appendChild(alignment);
                }
                container.appendChild(alignmentContainer);
                return container;
            }

            function createChangeColorButton(alignments, lenght, colored, self){
                var button = document.createElement('button');
                button.id                = 'changeColors';
                button.className         = 'btn';
                button.style.marginRight = '10px';
                button.style.marginTop   = '5px';
                if(colored == true){
                    var text = document.createTextNode('Click to change colors to grayscale');
                }else{
                    var text = document.createTextNode('Click to change colors to full colored');
                }
                button.appendChild(text);
                button.onclick=function(){ changeColors(alignments, lenght, button, colored, self); };
                return button;
            }

            function createDownloadButton(){
                var button = document.createElement('button');
                button.id              = 'downloadAlignments';
                button.className       = 'btn';
                button.style.marginTop = '5px';
                var text = document.createTextNode('Download as image'); 
                button.appendChild(text); 
                button.addEventListener('click', downloadAlignmentsImg);
                return button;
            }

            function changeColors(alignments, lenght, button, colored, self){   
                if(colored == true){
                    colored = false;
                }else{
                    colored = true;   
                }
                button.removeChild(button.childNodes[0]);
                var blastDiv = document.getElementById(self.opt.multipleAlignments);
                while (blastDiv.firstChild) {
                    blastDiv.removeChild(blastDiv.firstChild);
                }
                createAlignmentComparison(alignments, lenght, colored, self);
            }

            function downloadAlignmentsImg(){
                var buttons   = document.getElementById('blast-multiple-alignments-buttons');
                var input     = document.getElementById('defline');
                var container = document.getElementById('alignments-container');
                container.removeChild(buttons);
                container.removeChild(input);
                html2canvas(container, {
                  onrendered: function(canvas) {
                    document.body.appendChild(canvas);
                    var a = document.createElement('a');
                    document.body.appendChild(a);
                    a.href     = canvas.toDataURL('img/png');
                    a.download = 'alignments.png';
                    a.click();
                    document.body.removeChild(canvas);
                    document.body.removeChild(a);
                    container.appendChild(input);
                    container.appendChild(buttons);
                  }
                });   
            }

            function createAlignmentComparison(alignments, queryLenght, colored, self){
                var blastDiv  = document.getElementById(self.opt.multipleAlignments);
                while(blastDiv.hasChildNodes()){
                    blastDiv.removeChild(blastDiv.firstChild);	
                }
                var container        = document.createElement('div');
                var containerButtons = document.createElement('div');
                var input            = document.createElement('input');
                var buttonColors     = createChangeColorButton(alignments, queryLenght, colored, self);
                var buttonDownload   = createDownloadButton();
                blastDiv.style.paddingTop = '20px';
                input.id    = 'defline';
                input.name  = 'defline';
                input.type  = 'text';
                input.value = ' Mouse over to show defline and scores, click to show alignments';
                input.style.width   = '556px';
                input.style.padding = '1px';
                input.style.border  = 0;
                input.style.cursos  = 'auto';    
                container.id                    = 'alignments-container';
                container.style.border          = 'thin solid #DDD';
                container.style.margin          = '0 auto';
                container.style.padding         = '10px';
                container.style.maxWidth        = '580px';
                container.style.backgroundColor = '#FFF';
                container = createHeader(container, colored, queryLenght);
                container = createBody(alignments, queryLenght, container, colored);
                container.appendChild(input);    
                containerButtons.style.textAlign = 'right';
                containerButtons.id              = 'blast-multiple-alignments-buttons';
                blastDiv.style.minWidth        = '580px';
                containerButtons.appendChild(buttonColors); 
                containerButtons.appendChild(buttonDownload);   
                container.appendChild(containerButtons);
                blastDiv.appendChild(container);
            }

            function createTableHeader(){
                var table       = document.createElement('table');
                var thead       = document.createElement('thead');
                var theadTr     = document.createElement('tr');
                var theadDesc   = document.createElement('th');
                var theadScore  = document.createElement('th');
                var theadEval   = document.createElement('th');
                var theadIdent  = document.createElement('th');
                var theadPos    = document.createElement('th');
                var theadGaps   = document.createElement('th');
                table.className = "table table-striped";
                theadDesc.innerHTML  = 'Description'.bold();
                theadScore.innerHTML = 'Score'.bold();
                theadEval.innerHTML  = 'E value'.bold();
                theadIdent.innerHTML = 'Identities'.bold();
                theadPos.innerHTML   = 'Positives'.bold();
                theadGaps.innerHTML  = 'Gaps'.bold();
                theadTr.appendChild(theadDesc);
                theadTr.appendChild(theadScore);
                theadTr.appendChild(theadEval);
                theadTr.appendChild(theadIdent);
                theadTr.appendChild(theadPos);
                theadTr.appendChild(theadGaps);
                thead.appendChild(theadTr);
                table.appendChild(thead);
                return table;
            }

            function createTableButtons(alignments){
                var container = document.createElement('div');
                var butCSV    = document.createElement('button');
                var butImg    = document.createElement('button');
                container.style.textAlign = 'right';    
                butCSV.style.marginRight = '10px';
                butCSV.className         = 'btn';
                butCSV.innerHTML         = 'Download as csv';
                butCSV.onclick           = function(){ downloadTableCSV(alignments); };
                butImg.className = 'btn';
                butImg.innerHTML = 'Download as image';
                butImg.onclick   = function(){ downloadTableImg(); };
                container.appendChild(butCSV);
                container.appendChild(butImg);
                return container;
            }

            function createAlignmentTable(alignments, self){
                var tableDiv     = document.getElementById(self.opt.alignmentsTable);
                while(tableDiv.hasChildNodes()){
                    tableDiv.removeChild(tableDiv.firstChild);	
                }
                var imgContainer = document.createElement('div');
                var butContainer = createTableButtons(alignments);
                var table        = createTableHeader();
                var tbody        = document.createElement('tbody');
                tableDiv.style.paddingTop = '50px';
                imgContainer.style.backgroundColor = '#FFF';
                imgContainer.id                    = 'blast-alignments-table-img';
                for(var i = 0; i < alignments.length; i++){
                    var tr           = document.createElement('tr');
                    var tdDesc       = document.createElement('td');
                    var butDesc      = document.createElement('button');
                    butDesc.alignment = alignments[i];
                    butDesc.onclick   = function(){ location.href='#'+self.opt.singleAlignments;
                                                    createSingleAlignment(this.alignment, self);};
                    butDesc.id        = alignments[i].description.split(" ")[0];
                    butDesc.innerHTML = alignments[i].description;
                    butDesc.style.border = 0;
                    butDesc.style.padding = 0;
                    butDesc.style.display = 'inline';
                    butDesc.style.background = 'none';
                    butDesc.className = 'alignment-table-description';
                    tdDesc.appendChild(butDesc);
                    var tdScore       = document.createElement('td');
                    tdScore.innerHTML = alignments[i].score;
                    var tdEval       = document.createElement('td');
                    tdEval.innerHTML = alignments[i].eValue;
                    var tdIdent       = document.createElement('td');
                    tdIdent.innerHTML = alignments[i].identities+"%";
                    var tdPos       = document.createElement('td');
                    if(alignments[i].positives == "N/A"){
                        tdPos.innerHTML = alignments[i].positives;
                    }else{
                        tdPos.innerHTML = alignments[i].positives+"%";
                    }
                    var tdGaps       = document.createElement('td');
                    tdGaps.innerHTML = alignments[i].gaps+"%";
                    tr.appendChild(tdDesc);
                    tr.appendChild(tdScore);
                    tr.appendChild(tdEval);
                    tr.appendChild(tdIdent);
                    tr.appendChild(tdPos);
                    tr.appendChild(tdGaps);
                    tbody.appendChild(tr);
                }
                table.appendChild(tbody);
                imgContainer.appendChild(table);
                tableDiv.appendChild(imgContainer);
                tableDiv.appendChild(butContainer);
            }

            function downloadTableImg(){
                var items = document.getElementsByClassName('alignment-table-description');
                var i;
                for (i = 0; i < items.length; i++) {
                    items[i].style.fontWeight = 'normal';
                    items[i].parentElement.parentElement.style.fontWeight = 'normal';
                }
                var container = document.getElementById('blast-alignments-table-img');
                html2canvas(container, {
                  onrendered: function(canvas) {
                    document.body.appendChild(canvas);
                    var a = document.createElement('a');
                    document.body.appendChild(a);
                    a.href = canvas.toDataURL('img/png');
                    a.download = 'alignments-table.png';
                    a.click();
                    document.body.removeChild(canvas);
                    document.body.removeChild(a);
                  }
                });   
            }

            function downloadTableCSV(alignments){
                var csvContent = 'data:text/csv;charset=utf-8,';
                csvContent += 'Description; Score; eValue; Identities; Positives; Gaps\n';

                for(var i = 0; i < alignments.length; i++){
                    csvContent += alignments[i].description;
                    csvContent += '; ';
                    csvContent += alignments[i].score;
                    csvContent += '; ';
                    csvContent += alignments[i].eValue;
                    csvContent += '; ';
                    csvContent += alignments[i].identities;
                    csvContent += '; ';
                    csvContent += alignments[i].positives;
                    csvContent += '; ';
                    csvContent += alignments[i].gaps;
                    csvContent += '\n';
                }
                var encodedUri = encodeURI(csvContent);
                var link = document.createElement('a');
                link.setAttribute('href', encodedUri);
                link.setAttribute('download', 'alignments-table.csv');
                link.click(); 
            }

            function createSingleAlignment(alignment, self){
                var alignmentDiv = document.getElementById(self.opt.singleAlignment);
                while(alignmentDiv.hasChildNodes()){
                    alignmentDiv.removeChild(alignmentDiv.firstChild);	
                }

                var alignmentPre = document.createElement('pre');
                var alignmentContainer = document.createElement('div');
                var buttonDownload = document.createElement('button');
                var span = document.createElement('span');
                var queryDiv      = createSingleQueryDiv(alignment.query, alignment.queryStart, alignment.queryEnd);
                var comparisonDiv = createSingleComparisonDiv(alignment.comparison);
                var subjectDiv    = createSingleSubjectDiv(alignment.subject, alignment.subjectStart, alignment.subjectEnd); 

                alignmentPre.style.color         = '#2c3e50';
                alignmentPre.style.paddingTop    = '25px';
                alignmentPre.style.paddingBottom = '40px';
                alignmentPre.style.textAlign     = 'left';
                alignmentPre.style.fontFamily    = 'Helvetica,Arial,sans-serif';
                alignmentPre.id                  = 'blast-single-alignment-pre';
                alignmentContainer.style.margin     = '0 auto';
                alignmentContainer.style.display    = 'table';
                alignmentContainer.style.paddingTop = '30px';
                alignmentDiv.style.textAlign  = 'right';
                alignmentDiv.style.paddingTop = '50px';
                buttonDownload.className = 'btn';
                buttonDownload.innerHTML = 'Download as image';
                buttonDownload.onclick   = function(){ downloadSingleAlignmentImg(alignment); };
                span.innerHTML         = alignment.description;
                span.style.paddingLeft = '15px';
                span.style.fontWeight  = 'bold';
                span.style.fontSize    = '15px';
                span.style.fontFamily  =  'Helvetica,Arial,sans-serif';

                alignmentContainer.appendChild(queryDiv);
                alignmentContainer.appendChild(comparisonDiv);
                alignmentContainer.appendChild(subjectDiv);
                alignmentPre.appendChild(span);
                alignmentPre.appendChild(alignmentContainer);
                alignmentDiv.appendChild(alignmentPre);
                alignmentDiv.appendChild(buttonDownload);
            }

            function downloadSingleAlignmentImg(alignment){
                var container = document.getElementById('blast-single-alignment-pre');
                html2canvas(container, {
                  onrendered: function(canvas) {
                    document.body.appendChild(canvas);
                    var a = document.createElement('a');
                    document.body.appendChild(a);
                    a.href = canvas.toDataURL('img/png');
                    var tittle = alignment.description+'-alignment.png';
                    a.download = tittle;
                    a.click();
                    document.body.removeChild(canvas);
                    document.body.removeChild(a);
                  }
                });   
            }

            function createSingleQueryDiv(query, start, end){
                var textDiv  = document.createElement('div');
                var startDiv = document.createElement('div');
                var endDiv   = document.createElement('div');
                var queryDiv = document.createElement('div');
                textDiv.innerHTML         = 'Query'.bold();
                textDiv.style.display     = 'inline-block';
                textDiv.style.marginRight = '20px';
                textDiv.style.textAlign   = 'right';
                textDiv.style.width       = '55px';
                startDiv.innerHTML = String(start).bold();
                startDiv.style.display = 'inline-block';
                startDiv.style.marginRight = '20px';
                startDiv.style.width       = '25px';
                endDiv.innerHTML         = String(end).bold();
                endDiv.style.display     = 'inline-block';
                endDiv.style.marginLeft  = '20px';
                endDiv.style.marginRight = '70px';
                queryDiv.appendChild(textDiv);
                queryDiv.appendChild(startDiv);
                for(var i = 0; i < query.length; i++){
                    var div = document.createElement('div');
                    div.style.backgroundColor = getAminoColor(query.charAt(i));
                    div.innerHTML             = query.charAt(i).bold();
                    div.style.width           = '18px';
                    div.style.textAlign       = 'center';
                    div.style.display         = 'inline-block';
                    queryDiv.appendChild(div);
                }
                queryDiv.appendChild(endDiv);
                return queryDiv;
            }

            function createSingleComparisonDiv(comparison){
                var comparisonDiv = document.createElement('div');
                var spaceDiv      = document.createElement('div');
                spaceDiv.style.minWidth  = '120px';
                spaceDiv.style.minHeight = '1px';
                spaceDiv.style.display   = 'inline-block';
                comparisonDiv.appendChild(spaceDiv);
                for(var i = 0; i < comparison.length; i++){
                    var div = document.createElement('div');
                    div.style.backgroundColor = getAminoColor(comparison.charAt(i));
                    div.innerHTML             = comparison.charAt(i).bold();
                    div.style.width           = '18px';
                    div.style.textAlign       = 'center';
                    div.style.display         = 'inline-block';
                    comparisonDiv.appendChild(div);
                }
                return comparisonDiv;
            }

            function createSingleSubjectDiv(subject, start, end){
                var textDiv    = document.createElement('div');
                var startDiv   = document.createElement('div');
                var endDiv     = document.createElement('div');
                var subjectDiv = document.createElement('div');
                textDiv.innerHTML         = 'Subject'.bold();
                textDiv.style.display     = 'inline-block';
                textDiv.style.textAlign   = 'right';
                textDiv.style.marginRight = '20px';
                textDiv.style.width       = '55px';
                startDiv.style.width       = '25px';
                startDiv.innerHTML         = String(start).bold();
                startDiv.style.display     = 'inline-block';
                startDiv.style.marginRight = '20px';
                endDiv.innerHTML         = String(end).bold();
                endDiv.style.display     = 'inline-block';
                endDiv.style.marginLeft  = '20px';
                endDiv.style.marginRight = '70px';
                subjectDiv.appendChild(textDiv);
                subjectDiv.appendChild(startDiv);
                for(var i = 0; i < subject.length; i++){
                    var div = document.createElement('div');
                    div.style.backgroundColor = getAminoColor(subject.charAt(i));
                    div.innerHTML             = subject.charAt(i).bold();
                    div.style.width           = '18px';
                    div.style.textAlign       = 'center';
                    div.style.display         = 'inline-block';
                    subjectDiv.appendChild(div);
                }
                subjectDiv.appendChild(endDiv);
                return subjectDiv;
            }

            function getAminoColor(char){
                switch(char) {
                    case 'A':
                        return '#DBFA60';
                        break;
                    case 'C':
                        return '#F9FA60';
                        break;
                    case 'D':
                        return '#F9605F';
                        break;
                    case 'E':
                        return '#F9609C';
                        break;
                    case 'F':
                        return '#5FF99D';
                        break;
                    case 'G':
                        return '#F9BC5F';
                        break;
                    case 'H':
                        return '#609DF9';
                        break;
                    case 'I':
                        return '#99F95A';
                        break;
                    case 'K':
                        return '#A062FF';
                        break;
                    case 'L':
                        return '#7EF960';
                        break;
                    case 'M':
                        return '#63FF63';
                        break;
                    case 'N':
                        return '#D95DF9';
                        break;
                    case 'P':
                        return '#F9DA60';
                        break;
                    case 'Q':
                        return '#F955D8';
                        break;
                    case 'R':
                        return '#5360FB';
                        break;
                    case 'S':
                        return '#F97E60';
                        break;
                    case 'T':
                        return '#FFA563';
                        break;
                    case 'V':
                        return '#C0F86B';
                        break;
                    case 'W':
                        return '#FDD9F9';
                        break;
                    case 'Y':
                        return '#60F9DA';
                        break;
                    default:
                        return '#FFFFFF';
                }
            }
        }
});

},{"js-class":1,"js-extend":2}]},{},[]);
