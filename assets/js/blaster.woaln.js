require = (function e(t, n, l) {
    function i(s, a) {
        if (!n[s]) {
            if (!t[s]) {
                var p = "function" == typeof require && require;
                if (!a && p) return p(s, !0);
                if (r) return r(s, !0);
                var o = Error("Cannot find module '" + s + "'");
                throw ((o.code = "MODULE_NOT_FOUND"), o);
            }
            var d = (n[s] = { exports: {} });
            t[s][0].call(
                d.exports,
                function (e) {
                    return i(t[s][1][e] || e);
                },
                d,
                d.exports,
                e,
                t,
                n,
                l
            );
        }
        return n[s].exports;
    }
    for (var r = "function" == typeof require && require, s = 0; s < l.length; s++) i(l[s]);
    return i;
})(
    {
        1: [
            function (e, t, n) {
                (function (e) {
                    function n(e, t, n) {
                        for (
                            ;
                            t &&
                            t !== Object.prototype &&
                            (Object.getOwnPropertyNames(t).forEach(function (n) {
                                if (".class" != n && !e.hasOwnProperty(n)) {
                                    var l = Object.getOwnPropertyDescriptor(t, n);
                                    Object.defineProperty(e, n, l);
                                }
                            }),
                            !n);

                        )
                            t = t.__proto__;
                        return e;
                    }
                    var l = function (e, t, n) {
                        "function" != typeof e && ((n = t), (t = e), (e = Object)), (t = t || {});
                        var i = { name: (n = n || {}).name, base: e, implements: [] },
                            r = l.clone(t);
                        function s() {
                            "function" == typeof this.constructor && this.constructor.apply(this, arguments);
                        }
                        return (
                            n.implements &&
                                (Array.isArray(n.implements) ? n.implements : [n.implements]).forEach(function (e) {
                                    "function" == typeof e && e.prototype && (i.implements.push(e), l.extend(r, e.prototype));
                                }),
                            (r.__proto__ = e.prototype),
                            ((i.type = s).prototype = r),
                            Object.defineProperty(s, ".class.meta", { value: i, enumerable: !1, configurable: !1, writable: !1 }),
                            Object.defineProperty(r, ".class", { value: s, enumerable: !1, configurable: !1, writable: !1 }),
                            n.statics && l.extend(s, n.statics),
                            s
                        );
                    };
                    (l.extend = n),
                        (l.clone = function (e) {
                            return n({}, e);
                        });
                    var i = l({
                        constructor: function (e) {
                            this.object = e;
                        },
                        typeOf: function (e) {
                            if (this.object instanceof e) return !0;
                            var t = l.typeInfo(this.object);
                            return (
                                t &&
                                (function e(t, n) {
                                    for (; t; ) {
                                        if (t.type.prototype === n.prototype) return !0;
                                        for (var l in t.implements) {
                                            var i = t.implements[l],
                                                r = i[".class.meta"];
                                            if (r) {
                                                if (e(r, n)) return !0;
                                            } else for (var s = i.prototype; s; s = s.__proto__) if (s === n.prototype) return !0;
                                        }
                                        t = t.base ? t.base[".class.meta"] : void 0;
                                    }
                                    return !1;
                                })(t, e)
                            );
                        },
                    });
                    (i.prototype.a = i.prototype.typeOf),
                        (i.prototype.an = i.prototype.typeOf),
                        (l.is = function (e) {
                            return new i(e);
                        }),
                        (l.typeInfo = function (e) {
                            var t = e.__proto__[".class"];
                            return t ? t[".class.meta"] : void 0;
                        }),
                        (l.VERSION = [0, 0, 2]),
                        t ? (t.exports = l) : (e.Class = l);
                }.call(this, "undefined" != typeof global ? global : "undefined" != typeof self ? self : "undefined" != typeof window ? window : {}));
            },
            {},
        ],
        2: [
            function (e, t, n) {
                (function () {
                    var e = Array.prototype.slice,
                        t = Array.prototype.forEach,
                        n = function (l) {
                            if ("object" != typeof l) throw l + " is not an object";
                            var i = e.call(arguments, 1);
                            return (
                                t.call(i, function (e) {
                                    if (e) for (var t in e) "object" == typeof e[t] && l[t] ? n.call(l, l[t], e[t]) : (l[t] = e[t]);
                                }),
                                l
                            );
                        };
                    this.extend = n;
                }.call(this));
            },
            {},
        ],
        "biojs-vis-blasterjs": [
            function (e, t, n) {
                var l = e("js-class");
                t.exports = l({
                    constructor: function (t) {
                        var n = this;
                        e("js-extend").extend(this.opt, t);
                        var l = "table tbody tr td button:hover{ text-decoration: underline;}",
                            i = document.createElement("style");
                        if ((i.styleSheet ? (i.styleSheet.cssText = l) : i.appendChild(document.createTextNode(l)), document.getElementsByTagName("head")[0].appendChild(i), n.opt.string)) {
                            var r = { target: { files: [new Blob([n.opt.string], { type: "text/plain" })] } };
                            n._displayAlignments(r, n);
                        } else
                            document.getElementById(n.opt.input).addEventListener(
                                "change",
                                function (e) {
                                    n._displayAlignments(e, n);
                                },
                                !1
                            );
                    },
                    opt: { input: "blastinput", multipleAlignments: "blast-multiple-alignments", alignmentsTable: "blast-alignments-table", singleAlignment: "blast-single-alignment" },
                    _displayAlignments: function (e, t) {
                        if (window.File && window.FileReader && window.FileList && window.Blob) {
                            var n = e.target.files[0];
                            if (n) {
                                var l = new FileReader();
                                (l.onload = function (e) {
                                    try {
                                        var n,
                                            l = (n = e.target.result).split("\n")[2].startsWith("<BlastOutput>")
                                                ? (function (e) {
                                                      for (var t = [], n = e.split("\n"), l = 0, i = 0; i < n.length; i++) n[i].startsWith("<Iteration>") && l++;
                                                      if (1 == l) t.push(e);
                                                      else {
                                                          var r = 0,
                                                              s = "";
                                                          for (i = 0; i < n.length && !n[(r = i)].startsWith("<Iteration>"); i++) s = s + n[i] + "\n";
                                                          for (var a = 0; a < l; a++) {
                                                              var p = s + n[r] + "\n";
                                                              for (r++; void 0 !== n[r] && !n[r].startsWith("<Iteration>"); ) (p = p + n[r] + "\n"), r++;
                                                              t.push(p);
                                                          }
                                                      }
                                                      return t;
                                                  })(n)
                                                : (function (e) {
                                                      for (var t = [], n = e.split("\n"), l = 0, i = 0; i < n.length; i++) n[i].startsWith("Query=") && l++;
                                                      if (1 == l) t.push(e);
                                                      else {
                                                          var r = 0,
                                                              s = "";
                                                          for (i = 0; i < n.length && !n[(r = i)].startsWith("Query="); i++) s = s + n[i] + "\n";
                                                          for (var a = 0; a < l; a++) {
                                                              var p = s + n[r] + "\n";
                                                              for (r++; void 0 !== n[r] && !n[r].startsWith("Query="); ) (p = p + n[r] + "\n"), r++;
                                                              (p += "\nend\n"), t.push(p);
                                                          }
                                                      }
                                                      return t;
                                                  })(n),
                                            i = a(l[0]),
                                            r = s(l[0]);
                                        if (0 == i.length) v(l, [], 0, r, !0, !0, t);
                                        else {
                                            var p = i[0].hsp.length;
                                            b(l, p, 0, i, r, !0, !0, t), C(i, t);
                                        }
                                    } catch (o) {
                                        for (var d = ["blast-multiple-alignments", "blast-alignments-table"], c = 0; c < d.length; c++) for (var h = document.getElementById(d[c]); h.firstChild; ) h.removeChild(h.firstChild);
                                        alert("ERROR WHILE UPLOADING DATA: You have uploaded an invalid BLAST output file.");
                                    }
                                }),
                                    l.readAsText(n);
                            } else alert("ERROR WHILE UPLOADING DATA: Failed to load file.");
                        } else alert("The File APIs are not fully supported by your browser.");
                        function i(e, t, n, l, i) {
                            (this.description = e), (this.length = t), (this.totalScore = n), (this.queryCover = l), (this.hsp = i);
                        }
                        function r(e, t, n, l, i, r, s, a, p, o, d, c) {
                            (this.score = e),
                                (this.eValue = t),
                                (this.identities = n),
                                (this.positives = l),
                                (this.gaps = i),
                                (this.queryStart = r),
                                (this.query = s),
                                (this.queryEnd = a),
                                (this.comparison = p),
                                (this.subjectStart = o),
                                (this.subject = d),
                                (this.subjectEnd = c);
                        }
                        function s(e) {
                            return e.split("\n")[2].startsWith("<BlastOutput>")
                                ? (function (e) {
                                      for (var t = e.split("\n"), n = 0, l = 0; l < t.length; l++)
                                          if (t[l].includes("<Iteration_query-len>")) {
                                              n = t[l].split(">")[1].split("</")[0];
                                              break;
                                          }
                                      return n;
                                  })(e)
                                : (function (e) {
                                      for (var t = e.split("\n"), n = 0, l = 0; l < t.length; l++)
                                          if (t[l].startsWith("Length=")) {
                                              n = t[l].split("=")[1];
                                              break;
                                          }
                                      return n;
                                  })(e);
                        }
                        function a(e) {
                            return e.split("\n")[2].startsWith("<BlastOutput>")
                                ? (function (e) {
                                      for (var t = e.split("\n"), n = [], l = 0; l < t.length; l++)
                                          if (t[l].startsWith("<Hit>")) {
                                              for (var a = "", o = l; o < t.length && ((a += t[o]), !t[o].includes("</Hit>")); o++);
                                              for (
                                                  var d = a.split("<Hit_id>")[1].split("</")[0],
                                                      c = a.split("<Hit_def>")[1].split("</")[0],
                                                      h = d.concat(" ").concat(c),
                                                      u = a.split("<Hit_len>")[1].split("</")[0],
                                                      f = [],
                                                      g = a.split("<Hit_hsps>")[1].split("</Hit_hsps>")[0].split("</Hsp>"),
                                                      m = 0;
                                                  m < g.length - 1;
                                                  m++
                                              ) {
                                                  var y = g[m].split("<Hsp_bit-score>")[1].split("</")[0],
                                                      $ = g[m].split("<Hsp_evalue>")[1].split("</")[0],
                                                      v = g[m].split("<Hsp_identity>")[1].split("</")[0],
                                                      b = g[m].split("<Hsp_align-len>")[1].split("</")[0],
                                                      C = (v / b) * 100;
                                                  if (t[3].includes("<BlastOutput_program>blastn</BlastOutput_program>"))
                                                      var E = "N/A",
                                                          _ = (g[m].split("<Hsp_gaps>")[1].split("</")[0] / b) * 100;
                                                  else (E = parseFloat((g[m].split("<Hsp_positive>")[1].split("</")[0] / b) * 100).toFixed(0)), (_ = (g[m].split("<Hsp_gaps>")[1].split("</")[0] / b) * 100);
                                                  var x = g[m].split("<Hsp_query-from>")[1].split("</")[0],
                                                      T = g[m].split("<Hsp_qseq>")[1].split("</")[0],
                                                      H = g[m].split("<Hsp_query-to>")[1].split("</")[0],
                                                      L = g[m].split("<Hsp_midline>")[1].split("</")[0],
                                                      F = g[m].split("<Hsp_hit-from>")[1].split("</")[0],
                                                      A = g[m].split("<Hsp_hseq>")[1].split("</")[0],
                                                      w = g[m].split("<Hsp_hit-to>")[1].split("</")[0],
                                                      B = new r(y, $, parseFloat(C).toFixed(0), E, parseFloat(_).toFixed(0), x, T, H, L, F, A, w);
                                                  f.push(B);
                                              }
                                              for (var M = parseFloat(f[0].score), D = 1; D < f.length; D++) M += parseFloat(f[D].score);
                                              var N = new i(h, u, M.toFixed(1), p(f, s(e)), f);
                                              n.push(N);
                                          }
                                      return n;
                                  })(e)
                                : (function (e) {
                                      for (var t = e.split("\n"), n = [], l = 0; l < t.length; l++)
                                          if (t[l].startsWith(">")) {
                                              for (var a = t[l].split(">")[1], o = "", d = l; "" == o; ) t[(d += 1)].startsWith("Length=") ? (o = t[d]) : (a += t[d]);
                                              var c = a,
                                                  h = o.split("=")[1],
                                                  u = [],
                                                  f = !1;
                                              do {
                                                  if ((f && (d -= 1), t[d + 2].startsWith(" Features in this part of subject sequence:"))) for (d += 3; !t[d + 2].startsWith(" Score ="); ) d++;
                                                  var g = t[d + 2].split(",")[0].replace(/\s\s+/g, " ").split(" ")[3],
                                                      m = t[d + 2].split(",")[1].split(" ")[4],
                                                      y = t[d + 3]
                                                          .split(",")[0]
                                                          .split("(")[1]
                                                          .substr(0, t[d + 3].split(",")[0].split("(")[1].length - 2);
                                                  if (t[0].startsWith("BLASTN"))
                                                      var $ = "N/A",
                                                          v = t[d + 3]
                                                              .split(",")[1]
                                                              .split("(")[1]
                                                              .substr(0, t[d + 3].split(",")[1].split("(")[1].length - 2);
                                                  else
                                                      ($ = t[d + 3]
                                                          .split(",")[1]
                                                          .split("(")[1]
                                                          .substr(0, t[d + 3].split(",")[1].split("(")[1].length - 2)),
                                                          (v = t[d + 3]
                                                              .split(",")[2]
                                                              .split("(")[1]
                                                              .substr(0, t[d + 3].split(",")[2].split("(")[1].length - 2));
                                                  ("Frame" == t[d + 4].split(",")[0].split(" ")[1] || t[d + 4].startsWith(" Strand")) && (d += 1);
                                                  var b = t[d + 5].substring(5).replace(/^\s+/g, "").split(" ")[0],
                                                      C = t[d + 5].substring(5).replace(/\s+/g, "").replace(/[0-9]/g, ""),
                                                      E = t[d + 5].substring(5).replace(/^\s+/g, "").split(" ")[t[d + 5].substring(5).replace(/^\s+/g, "").split(" ").length - 1],
                                                      _ = t[d + 6].replace(/^\s+/g, ""),
                                                      x = t[d + 7].substring(5).replace(/^\s+/g, "").split(" ")[0],
                                                      T = t[d + 7].substring(5).replace(/\s+/g, "").replace(/[0-9]/g, ""),
                                                      H = t[d + 7].substring(5).replace(/^\s+/g, "").split(" ")[t[d + 7].substring(5).replace(/^\s+/g, "").split(" ").length - 1];
                                                  for (d += 9; t[d].startsWith("Query"); ) {
                                                      var L = t[d].substring(5).replace(/\s+/g, "").replace(/[0-9]/g, "");
                                                      (C += L),
                                                          (E = t[d].substring(5).replace(/^\s+/g, "").split(" ")[t[d].substring(5).replace(/^\s+/g, "").split(" ").length - 1]),
                                                          (T += t[d + 2].substring(5).replace(/\s+/g, "").replace(/[0-9]/g, "")),
                                                          (H = t[d + 2].substring(5).replace(/^\s+/g, "").split(" ")[t[d + 2].substring(5).replace(/^\s+/g, "").split(" ").length - 1]);
                                                      var F = t[d + 1].replace(/^\s+/g, "");
                                                      if (L.length > F.length) for (var A = L.length - F.length, w = 0; w < A; w++) F = " " + F;
                                                      (_ += F), (d += 4);
                                                  }
                                                  var B = new r(g, m, y, $, v, b, C, E, _, x, T, H);
                                                  u.push(B), (f = !0);
                                              } while (t[d + 1].startsWith(" Score"));
                                              for (var M = parseFloat(u[0].score), D = 1; D < u.length; D++) M += parseFloat(u[D].score);
                                              var N = new i(c, h, M.toFixed(1), p(u, s(e)), u);
                                              n.push(N);
                                          }
                                      return n;
                                  })(e);
                        }
                        function p(e, t) {
                            for (var n = 0, l = o(e), i = 0; i < l.length; i++) n += parseInt((100 * (l[i].end - l[i].start + 1)) / t);
                            return n;
                        }
                        function o(e) {
                            for (var t = [], n = 0; n < e.length; n++)
                                parseInt(e[n].queryStart) > parseInt(e[n].queryEnd) ? t.push({ start: parseInt(e[n].queryEnd), end: parseInt(e[n].queryStart) }) : t.push({ start: parseInt(e[n].queryStart), end: parseInt(e[n].queryEnd) });
                            return (function (e) {
                                for (var t = [], n = 0; n < e.length; n++) {
                                    for (var l = e[n][0].start, i = e[n][0].end, r = 0; r < e[n].length; r++) e[n][r].start < l && (l = e[n][r].start), e[n][r].end > i && (i = e[n][r].end);
                                    t.push({ start: l, end: i });
                                }
                                return t;
                            })(
                                (function (e) {
                                    e.sort(function (e, t) {
                                        return e.start < t.start ? -1 : e.start > t.start ? 1 : 0;
                                    });
                                    var t,
                                        n = [],
                                        l = 0;
                                    n[l] = [e[0]];
                                    for (var i = 1, r = e.length; i < r; i++)
                                        e[i].start >= e[i - 1].start &&
                                        e[i].start <
                                            (0 != (t = n[l]).length &&
                                                (t.sort(function (e, t) {
                                                    return e.end < t.end ? 1 : e.end > t.end ? -1 : 0;
                                                }),
                                                t[0].end))
                                            ? n[l].push(e[i])
                                            : (n[++l] = [e[i]]);
                                    return n;
                                })(t)
                            );
                        }
                        function d(e, t) {
                            if (e)
                                switch (t) {
                                    case 1:
                                        return "<40";
                                    case 2:
                                        return "40-50";
                                    case 3:
                                        return "50-80";
                                    case 4:
                                        return "80-200";
                                    case 5:
                                        return ">=200";
                                    default:
                                        return "0";
                                }
                            else
                                switch (t) {
                                    case 1:
                                        return ">100";
                                    case 2:
                                        return "100-1";
                                    case 3:
                                        return "1-1e<sup>-2</sup>";
                                    case 4:
                                        return "1e<sup>-2</sup>-1e<sup>-5</sup>";
                                    case 5:
                                        return "<1e<sup>-5</sup>";
                                    default:
                                        return "0";
                                }
                        }
                        function c(e, t) {
                            if (e)
                                switch (t) {
                                    case 1:
                                        return "#5C6D7E";
                                    case 2:
                                        return "#9B59B6";
                                    case 3:
                                        return "#5CACE2";
                                    case 4:
                                        return "#57D68D";
                                    case 5:
                                        return "#C0392B";
                                    default:
                                        return "#FFF";
                                }
                            else
                                switch (t) {
                                    case 1:
                                        return "#BCBCBC";
                                    case 2:
                                        return "#989898";
                                    case 3:
                                        return "#747474";
                                    case 4:
                                        return "#565656";
                                    case 5:
                                        return "#343434";
                                    default:
                                        return "#FFF";
                                }
                        }
                        function h(e) {
                            var t,
                                n = document.getElementById(e).parentElement.parentElement,
                                l = document.getElementsByClassName("alignment-table-description");
                            for (t = 0; t < l.length; t++) (l[t].style.fontWeight = "normal"), (l[t].parentElement.parentElement.style.fontWeight = "normal");
                            (n.style.fontWeight = "bold"), (document.getElementById(e).style.fontWeight = "bold");
                        }
                        function u(e, t, n, l) {
                            for (var i = 0; i < t; i++)
                                if (0 == i) {
                                    var r = l / 2,
                                        s = document.createElement("div"),
                                        a = document.createElement("div");
                                    (s.style.float = "left"),
                                        (s.style.width = r + "px"),
                                        (s.style.textAlign = "left"),
                                        (s.innerHTML = "0"),
                                        (a.style.float = "left"),
                                        (a.style.width = r + "px"),
                                        (a.style.textAlign = "right"),
                                        (a.innerHTML = n * (i + 1)),
                                        e.appendChild(s),
                                        e.appendChild(a);
                                } else {
                                    var p = document.createElement("div");
                                    (p.style.float = "left"), (p.style.width = l + "px"), (p.style.textAlign = "right"), (p.innerHTML = n * (i + 1)), e.appendChild(p);
                                }
                            return e;
                        }
                        function f(e, t, n, l) {
                            var i,
                                r,
                                s,
                                a,
                                p,
                                o,
                                h,
                                f,
                                g,
                                m,
                                y = document.createElement("div"),
                                $ = (function (e, t) {
                                    var n = document.createElement("div"),
                                        l = document.createElement("div"),
                                        i = document.createElement("div");
                                    (n.style.color = "#EEE"), (l.style.minWidth = "50px"), (l.style.minHeight = "10px"), (l.style.float = "left"), n.appendChild(l);
                                    for (var r = 1; r <= 5; r++) {
                                        var s = document.createElement("div");
                                        (s.style.minWidth = "100px"), (s.style.textAlign = "center"), (s.style.float = "left"), (s.innerHTML = d(t, r).bold()), (s.style.backgroundColor = c(e, r)), n.appendChild(s);
                                    }
                                    return (i.style.clear = "both"), n.appendChild(i), n;
                                })(t, n),
                                v =
                                    ((i = t),
                                    (r = document.createElement("div")),
                                    (s = document.createElement("div")),
                                    (a = document.createElement("div")),
                                    (p = document.createElement("div")),
                                    (r.style.marginTop = "3px"),
                                    (r.style.color = "#5C6D7E"),
                                    (r.style.fontSize = "10px"),
                                    (s.style.width = "50px"),
                                    (s.innerHTML = "QUERY".bold()),
                                    (s.style.float = "left"),
                                    (a.style.width = "500px"),
                                    (a.style.height = "10px"),
                                    (a.style.float = "left"),
                                    (a.style.marginTop = "2px"),
                                    (p.style.clear = "both"),
                                    (a.style.backgroundColor = i ? "#C0392B" : "#343434"),
                                    r.appendChild(s),
                                    r.appendChild(a),
                                    r.appendChild(p),
                                    r),
                                b =
                                    ((o = l),
                                    (h = document.createElement("div")),
                                    (f = document.createElement("div")),
                                    (g = document.createElement("div")),
                                    (m = document.createElement("div")),
                                    (h.style.marginBottom = "5px"),
                                    (h.style.fontSize = "11px"),
                                    (f.style.minWidth = "50px"),
                                    (f.style.minHeight = "10px"),
                                    (f.style.float = "left"),
                                    (g.style.float = "left"),
                                    (m.style.clear = "both"),
                                    (g = (function (e, t) {
                                        var n = document.createElement("div");
                                        if (4 < t) {
                                            if (t % 5 == 0) e = u(e, 5, t / 5, 100);
                                            else {
                                                var l = 500 / (5 + (t % 5) / 5);
                                                e = u(e, 5, parseInt(t / 5), parseInt(l));
                                                var i = parseInt(500 - 5 * l),
                                                    r = document.createElement("div");
                                                (r.style.float = "left"), (r.style.width = i + "px"), (r.style.textAlign = "right"), (r.innerHTML = t), e.appendChild(r);
                                            }
                                        } else e = u(e, t, 1, parseInt(500 / t));
                                        return (n.style.clear = "both"), e.appendChild(n), e;
                                    })(g, o)),
                                    h.appendChild(f),
                                    h.appendChild(g),
                                    h.appendChild(m),
                                    h);
                            return (
                                (y.style.color = "#5C6D7E"),
                                (y.style.textAlign = "center"),
                                (y.style.paddingBottom = "5px"),
                                (y.innerHTML = "COLOR KEY FOR ALIGNMENT SCORES".bold()),
                                e.appendChild(y),
                                e.appendChild($),
                                e.appendChild(v),
                                e.appendChild(b),
                                e
                            );
                        }
                        function g(e, t, n) {
                            var l = o(t.hsp),
                                i = document.createElement("div"),
                                r = document.createElement("div");
                            i.style.minHeight = "12px";
                            for (var s = 0; s < l.length; s++) {
                                var a = document.createElement("div"),
                                    p = document.createElement("div"),
                                    d = document.createElement("a");
                                if (0 == s) {
                                    if (1 == l[0].start) var c = parseInt(50 + (500 * (l[0].start - 1)) / n);
                                    else c = parseInt(50 + (500 * l[0].start) / n);
                                    var u = parseInt(550 - c - (500 * (n - l[0].end)) / n);
                                } else (c = parseInt((500 * (l[s].start - l[s - 1].end)) / n)), (u = parseInt((500 * (l[s].end - l[s].start)) / n));
                                (a.style.width = c + "px"),
                                    (a.style.minHeight = "4px"),
                                    (a.style.float = "left"),
                                    (p.style.width = u + "px"),
                                    (p.style.minHeight = "12px"),
                                    (p.style.float = "left"),
                                    (p.style.backgroundColor = e),
                                    (p.onmouseout = function () {
                                        document.getElementById("defline").value = " Mouse over to show defline and scores, click to show alignments";
                                    }),
                                    (p.onmouseover = function () {
                                        document.getElementById("defline").value = " " + t.description + ". S=" + t.hsp[0].score + " E=" + t.hsp[0].eValue;
                                    }),
                                    (d.href = "#" + t.description.split(" ")[0]),
                                    (d.onclick = function () {
                                        h(t.description.split(" ")[0]);
                                    }),
                                    d.appendChild(p),
                                    i.appendChild(a),
                                    i.appendChild(d);
                            }
                            return (r.style.clear = "both"), i.appendChild(r), i;
                        }
                        function m(e) {
                            var t = document.createElement("button");
                            (t.id = "downloadAlignments" + e), (t.className = "btn"), (t.style.marginRight = "10px"), (t.style.marginTop = "5px");
                            var n = document.createTextNode("Download as " + e);
                            return (
                                t.appendChild(n),
                                t.addEventListener(
                                    "click",
                                    function () {
                                        var t, n, l, i;
                                        (t = e),
                                            (n = document.getElementById("blast-multiple-alignments-buttons")),
                                            (l = document.getElementById("defline")),
                                            (i = document.getElementById("alignments-container")).removeChild(n),
                                            i.removeChild(l),
                                            html2canvas(i, {
                                                onrendered: function (e) {
                                                    document.body.appendChild(e);
                                                    var r = document.createElement("a");
                                                    document.body.appendChild(r),
                                                        (r.href = "JPEG" == t ? e.toDataURL("image/jpeg") : e.toDataURL("img/png")),
                                                        (r.download = "alignments." + t.toLowerCase()),
                                                        r.click(),
                                                        document.body.removeChild(e),
                                                        document.body.removeChild(r),
                                                        i.appendChild(l),
                                                        i.appendChild(n),
                                                        (r = document.createElement("a")),
                                                        document.body.appendChild(r),
                                                        (r.href = "#blast-multiple-alignments"),
                                                        r.click(),
                                                        document.body.removeChild(r);
                                                },
                                            });
                                    },
                                    !1
                                ),
                                t
                            );
                        }
                        function y(e) {
                            return e.split("\n")[2].startsWith("<BlastOutput>")
                                ? (function (e) {
                                      for (var t = e.split("\n"), n = 0; n < t.length; n++) if (t[n].includes("<Iteration_query-def>")) return t[n].split(">")[1].split("</")[0];
                                  })(e)
                                : (function (e) {
                                      for (var t = e.split("\n"), n = "", l = 0; l < t.length; l++)
                                          if (t[l].startsWith("Query=")) {
                                              for (n = t[l].split("=")[1], l++; !t[l].startsWith("Length="); ) (n = n + " " + t[l]), l++;
                                              break;
                                          }
                                      return n;
                                  })(e);
                        }
                        function $(e, t, n, l, i, r, p) {
                            var o = document.createElement("select");
                            (o.style.width = "430px"), (o.style.marginBottom = "20px"), (o.style.color = "#5C6D7E"), (o.style.float = "right");
                            for (var d = 0; d < e.length; d++) {
                                var c = document.createElement("option");
                                (c.value = d), (c.text = y(e[d])), d == n && (c.selected = "selected"), o.appendChild(c);
                            }
                            return (
                                (o.onchange = function () {
                                    !(function (e, t, n, l, i, r, p) {
                                        var o = a(e[n]);
                                        if (((l = s(e[n])), 0 == o.length)) {
                                            for (var d = ["blast-alignments-table", "blast-single-alignment"], c = 0; c < d.length; c++) for (var h = document.getElementById(d[c]); h.firstChild; ) h.removeChild(h.firstChild);
                                            v(e, [], n, l, !0, !0, p);
                                        } else b(e, t, n, o, l, i, r, p), C(o, p), (t = o[0].hsp.length), _(o[0], p, t, 0);
                                    })(e, t, o.value, l, i, r, p);
                                }),
                                o
                            );
                        }
                        function v(e, t, n, l, i, r, s) {
                            for (var a = document.getElementById(s.opt.multipleAlignments); a.hasChildNodes(); ) a.removeChild(a.firstChild);
                            var p = document.createElement("div"),
                                o = document.createElement("div"),
                                d = document.createElement("div");
                            if (
                                ((a.style.paddingTop = "20px"),
                                (p.id = "alignments-container"),
                                (p.style.border = "thin solid #DDD"),
                                (p.style.margin = "0 auto"),
                                (p.style.padding = "10px"),
                                (p.style.maxWidth = "580px"),
                                (p.style.backgroundColor = "#FFF"),
                                (d.innerHTML = "***NO HITS FOUND***"),
                                (d.style.fontWeight = "bold"),
                                (d.style.paddingTop = "30px"),
                                (d.style.paddingBottom = "50px"),
                                (d.style.textAlign = "center"),
                                1 < e.length)
                            ) {
                                var c = $(e, t, n, l, i, r, s),
                                    h = document.createElement("div"),
                                    u = document.createElement("div");
                                (u.style.clear = "both"), (h.innerHTML = "RESULTS FOR:".bold()), (h.style.marginBottom = "5px"), (h.style.color = "#5C6D7E"), (h.style.float = "left"), p.appendChild(h), p.appendChild(c), p.appendChild(u);
                            }
                            (p = f(p, i, r, l)).appendChild(d), (a.style.minWidth = "580px"), p.appendChild(o), a.appendChild(p);
                        }
                        function b(e, t, n, l, i, r, s, a) {
                            for (var p = document.getElementById(a.opt.multipleAlignments); p.hasChildNodes(); ) p.removeChild(p.firstChild);
                            var o = document.createElement("div"),
                                d = document.createElement("div"),
                                h = document.createElement("input"),
                                u = (function e(t, n, l, i, r, s, a, p) {
                                    var o = document.createElement("button");
                                    if (((o.id = "changeColors"), (o.className = "btn"), (o.style.marginRight = "10px"), (o.style.marginTop = "5px"), 1 == s)) var d = document.createTextNode("Change colours to grayscale");
                                    else d = document.createTextNode("Change colours to full colored");
                                    return (
                                        o.appendChild(d),
                                        (o.onclick = function () {
                                            !(function (e, t, n, l, i, r, s, a, p) {
                                                (s = 1 != s), r.removeChild(r.childNodes[0]);
                                                for (var o = document.getElementById(p.opt.multipleAlignments); o.firstChild; ) o.removeChild(o.firstChild);
                                                b(e, t, n, l, i, s, a, p);
                                            })(t, n, l, i, r, o, s, a, p);
                                        }),
                                        o
                                    );
                                })(e, t, n, l, i, r, s, a),
                                y = (function e(t, n, l, i, r, s, a, p) {
                                    var o = document.createElement("button");
                                    if (((o.id = "changeScore"), (o.className = "btn"), (o.style.marginRight = "10px"), (o.style.marginTop = "5px"), 1 == a)) var d = document.createTextNode("Change scoring to E value");
                                    else d = document.createTextNode("Change scoring to Max score");
                                    return (
                                        o.appendChild(d),
                                        (o.onclick = function () {
                                            !(function (e, t, n, l, i, r, s, a, p) {
                                                (a = 1 != a), r.removeChild(r.childNodes[0]);
                                                for (var o = document.getElementById(p.opt.multipleAlignments); o.firstChild; ) o.removeChild(o.firstChild);
                                                b(e, t, n, l, i, s, a, p);
                                            })(t, n, l, i, r, o, s, a, p);
                                        }),
                                        o
                                    );
                                })(e, t, n, l, i, r, s, a),
                                v = m("PNG"),
                                C = m("JPEG");
                            if (
                                ((p.style.paddingTop = "20px"),
                                (h.id = "defline"),
                                (h.name = "defline"),
                                (h.type = "text"),
                                (h.value = " Mouse over to show defline and scores, click to show alignments"),
                                (h.style.width = "556px"),
                                (h.style.padding = "1px"),
                                (h.style.border = 0),
                                (h.style.cursos = "auto"),
                                (o.id = "alignments-container"),
                                (o.style.border = "thin solid #DDD"),
                                (o.style.margin = "0 auto"),
                                (o.style.padding = "10px"),
                                (o.style.maxWidth = "580px"),
                                (o.style.backgroundColor = "#FFF"),
                                1 < e.length)
                            ) {
                                var E = $(e, t, n, i, r, s, a),
                                    _ = document.createElement("div"),
                                    x = document.createElement("div");
                                (x.style.clear = "both"), (_.innerHTML = "RESULTS FOR:".bold()), (_.style.marginBottom = "5px"), (_.style.color = "#5C6D7E"), (_.style.float = "left"), o.appendChild(_), o.appendChild(E), o.appendChild(x);
                            }
                            (o = (function (e, t, n, l, i) {
                                var r,
                                    s,
                                    a,
                                    p,
                                    o = document.createElement("div");
                                o.style.paddingBottom = "10px";
                                for (var d = 0; d < e.length; d++) {
                                    var h = g(
                                        ((r = l),
                                        (s = i),
                                        (a = e[d].hsp[0].score),
                                        (p = e[d].hsp[0].eValue),
                                        c(r, s ? (a < 40 ? 1 : 40 <= a && a < 50 ? 2 : 50 <= a && a < 80 ? 3 : 80 <= a && a < 200 ? 4 : 5) : 100 < p ? 1 : p <= 100 && 1 < p ? 2 : p <= 1 && 0.01 < p ? 3 : p <= 0.01 && 1e-5 < p ? 4 : 5)),
                                        e[d],
                                        t
                                    );
                                    (h.style.marginBottom = "4px"), o.appendChild(h);
                                }
                                return n.appendChild(o), n;
                            })(l, i, (o = f(o, r, s, i)), r, s)).appendChild(h),
                                (d.style.textAlign = "right"),
                                (d.id = "blast-multiple-alignments-buttons"),
                                (p.style.minWidth = "580px"),
                                d.appendChild(y),
                                d.appendChild(u),
                                d.appendChild(document.createElement("br")),
                                d.appendChild(v),
                                d.appendChild(C),
                                o.appendChild(d),
                                p.appendChild(o);
                        }
                        function C(e, t) {
                            for (var n = document.getElementById(t.opt.alignmentsTable); n.hasChildNodes(); ) n.removeChild(n.firstChild);
                            var l,
                                i,
                                r,
                                s,
                                a,
                                p,
                                o,
                                d,
                                c,
                                h,
                                u,
                                f,
                                g,
                                m,
                                y = document.createElement("div"),
                                $ =
                                    ((l = e),
                                    (i = document.createElement("div")),
                                    (r = document.createElement("button")),
                                    (s = document.createElement("button")),
                                    (a = document.createElement("button")),
                                    (i.style.textAlign = "right"),
                                    (r.style.marginRight = "10px"),
                                    (r.className = "btn"),
                                    (r.innerHTML = "Download as CSV"),
                                    (r.onclick = function () {
                                        !(function (e) {
                                            var t = "data:text/csv;charset=utf-8,";
                                            t += "Description; Score; eValue; Identities; Positives; Gaps\n";
                                            for (var n = 0; n < e.length; n++)
                                                (t += e[n].description),
                                                    (t += "; "),
                                                    (t += e[n].hsp[0].score),
                                                    (t += "; "),
                                                    (t += e[n].hsp[0].eValue),
                                                    (t += "; "),
                                                    (t += e[n].hsp[0].identities),
                                                    (t += "; "),
                                                    (t += e[n].hsp[0].positives),
                                                    (t += "; "),
                                                    (t += e[n].hsp[0].gaps),
                                                    (t += "\n");
                                            var l = encodeURI(t),
                                                i = document.createElement("a");
                                            i.setAttribute("href", l), i.setAttribute("download", "alignments-table.csv"), i.click();
                                        })(l);
                                    }),
                                    (a.className = "btn"),
                                    (a.innerHTML = "Download as PNG"),
                                    (a.onclick = function () {
                                        E("PNG");
                                    }),
                                    (a.style.marginRight = "10px"),
                                    (s.className = "btn"),
                                    (s.innerHTML = "Download as JPEG"),
                                    (s.onclick = function () {
                                        E("JPEG");
                                    }),
                                    i.appendChild(r),
                                    i.appendChild(a),
                                    i.appendChild(s),
                                    i),
                                v =
                                    ((p = document.createElement("table")),
                                    (o = document.createElement("thead")),
                                    (d = document.createElement("tr")),
                                    (c = document.createElement("th")),
                                    (h = document.createElement("th")),
                                    (u = document.createElement("th")),
                                    (f = document.createElement("th")),
                                    (g = document.createElement("th")),
                                    (m = document.createElement("th")),
                                    (subjectStartTh = document.createElement("th")),
                                    (subjectEndTh = document.createElement("th")),
                                    (lengthTh = document.createElement("th")),
                                    (p.className = "table table-striped"),
                                    (c.innerHTML = "Description".bold()),
                                    (h.innerHTML = "Max score".bold()),
                                    (u.innerHTML = "Total score".bold()),
                                    (f.innerHTML = "Query cover".bold()),
                                    (g.innerHTML = "E value".bold()),
                                    (m.innerHTML = "Identities".bold()),
                                    (subjectStartTh.innerHTML = "Subject Start".bold()),
                                    (subjectEndTh.innerHTML = "Subject End".bold()),
                                    (lengthTh.innerHTML = "Length".bold()),
                                    d.appendChild(c),
                                    d.appendChild(h),
                                    d.appendChild(u),
                                    d.appendChild(f),
                                    d.appendChild(g),
                                    d.appendChild(m),
                                    d.appendChild(subjectStartTh),
                                    d.appendChild(subjectEndTh),
                                    d.appendChild(lengthTh),
                                    o.appendChild(d),
                                    p.appendChild(o),
                                    p),
                                b = document.createElement("tbody");
                            (n.style.paddingTop = "50px"), (y.style.backgroundColor = "#FFF"), (y.id = "blast-alignments-table-img");
                            for (var C = 0; C < e.length; C++) {
                                var x = document.createElement("tr"),
                                    T = document.createElement("td"),
                                    H = document.createElement("button");
                                (H.alignment = e[C]),
                                    (H.onclick = function () {
                                        t.opt.callback ? t.opt.callback(this.alignment) : ((location.href = "#" + t.opt.singleAlignment), _(this.alignment, t, this.alignment.hsp.length, 0));
                                    }),
                                    (H.id = e[C].description.split(" ")[0]),
                                    (H.innerHTML = e[C].description),
                                    (H.style.border = 0),
                                    (H.style.padding = 0),
                                    (H.style.display = "inline"),
                                    (H.style.background = "none"),
                                    (H.className = "alignment-table-description"),
                                    T.appendChild(H);
                                var L = document.createElement("td");
                                L.innerHTML = e[C].hsp[0].score;
                                var F = document.createElement("td");
                                F.innerHTML = e[C].totalScore;
                                var A = document.createElement("td");
                                A.innerHTML = e[C].queryCover + "%";
                                var w = document.createElement("td");
                                w.innerHTML = e[C].hsp[0].eValue;
                                var B = document.createElement("td");
                                (B.innerHTML = e[C].hsp[0].identities + "%");
                                var subjectStartTd = document.createElement("td");
                                subjectStartTd.innerHTML = e[C].hsp[0].subjectStart;
                                var subjectEndTd = document.createElement("td");
                                subjectEndTd.innerHTML = e[C].hsp[0].subjectEnd;
                                var lengthTd = document.createElement("td");
                                lengthTd.innerHTML = e[C].length;
                                x.appendChild(T), x.appendChild(L), x.appendChild(F), x.appendChild(A), x.appendChild(w), x.appendChild(B), x.appendChild(subjectStartTd), x.appendChild(subjectEndTd), x.appendChild(lengthTd), b.appendChild(x);
                            }
                            v.appendChild(b), y.appendChild(v), n.appendChild(y), n.appendChild($);
                        }
                        function E(e) {
                            var t,
                                n = document.getElementsByClassName("alignment-table-description");
                            for (t = 0; t < n.length; t++) (n[t].style.fontWeight = "normal"), (n[t].parentElement.parentElement.style.fontWeight = "normal");
                            html2canvas(document.getElementById("blast-alignments-table-img"), {
                                onrendered: function (t) {
                                    document.body.appendChild(t);
                                    var n = document.createElement("a");
                                    document.body.appendChild(n),
                                        (n.href = "JPEG" == e ? t.toDataURL("image/jpeg") : t.toDataURL("img/png")),
                                        (n.download = "alignments-table." + e.toLowerCase()),
                                        n.click(),
                                        document.body.removeChild(t),
                                        document.body.removeChild(n),
                                        (n = document.createElement("a")),
                                        document.body.appendChild(n),
                                        (n.href = "#blast-alignments-table"),
                                        n.click(),
                                        document.body.removeChild(n);
                                },
                            });
                        }
                        function _(e, t, n, l) {
                            for (var i = document.getElementById(t.opt.singleAlignment); i.hasChildNodes(); ) i.removeChild(i.firstChild);
                            var r = document.createElement("pre"),
                                s = document.createElement("div"),
                                a = document.createElement("button"),
                                p = document.createElement("button"),
                                o = document.createElement("span"),
                                d = document.createElement("div"),
                                c = document.createElement("div"),
                                h = document.createElement("div"),
                                u = document.createElement("div"),
                                f = document.createElement("div"),
                                g = document.createElement("div"),
                                m = document.createElement("div"),
                                y = (function (e, t, n) {
                                    var l = document.createElement("div"),
                                        i = document.createElement("div"),
                                        r = document.createElement("div"),
                                        s = document.createElement("div");
                                    (l.innerHTML = "Query".bold()),
                                        (l.style.display = "inline-block"),
                                        (l.style.marginRight = "20px"),
                                        (l.style.textAlign = "right"),
                                        (l.style.width = "55px"),
                                        (i.innerHTML = String(t).bold()),
                                        (i.style.display = "inline-block"),
                                        (i.style.marginRight = "20px"),
                                        (i.style.width = "25px"),
                                        (r.innerHTML = String(n).bold()),
                                        (r.style.display = "inline-block"),
                                        (r.style.marginLeft = "20px"),
                                        (r.style.marginRight = "70px"),
                                        s.appendChild(l),
                                        s.appendChild(i);
                                    for (var a = 0; a < e.length; a++) {
                                        var p = document.createElement("div");
                                        (p.style.backgroundColor = T(e.charAt(a))), (p.innerHTML = e.charAt(a).bold()), (p.style.width = "18px"), (p.style.textAlign = "center"), (p.style.display = "inline-block"), s.appendChild(p);
                                    }
                                    return s.appendChild(r), s;
                                })(e.hsp[l].query, e.hsp[l].queryStart, e.hsp[l].queryEnd),
                                $ = (function (e) {
                                    var t = document.createElement("div"),
                                        n = document.createElement("div");
                                    (n.style.minWidth = "120px"), (n.style.minHeight = "1px"), (n.style.display = "inline-block"), t.appendChild(n);
                                    for (var l = 0; l < e.length; l++) {
                                        var i = document.createElement("div");
                                        (i.style.backgroundColor = T(e.charAt(l))), (i.innerHTML = e.charAt(l).bold()), (i.style.width = "18px"), (i.style.textAlign = "center"), (i.style.display = "inline-block"), t.appendChild(i);
                                    }
                                    return t;
                                })(e.hsp[l].comparison),
                                v = (function (e, t, n) {
                                    var l = document.createElement("div"),
                                        i = document.createElement("div"),
                                        r = document.createElement("div"),
                                        s = document.createElement("div");
                                    (l.innerHTML = "Subject".bold()),
                                        (l.style.display = "inline-block"),
                                        (l.style.textAlign = "right"),
                                        (l.style.marginRight = "20px"),
                                        (l.style.width = "55px"),
                                        (i.style.width = "25px"),
                                        (i.innerHTML = String(t).bold()),
                                        (i.style.display = "inline-block"),
                                        (i.style.marginRight = "20px"),
                                        (r.innerHTML = String(n).bold()),
                                        (r.style.display = "inline-block"),
                                        (r.style.marginLeft = "20px"),
                                        (r.style.marginRight = "70px"),
                                        s.appendChild(l),
                                        s.appendChild(i);
                                    for (var a = 0; a < e.length; a++) {
                                        var p = document.createElement("div");
                                        (p.style.backgroundColor = T(e.charAt(a))), (p.innerHTML = e.charAt(a).bold()), (p.style.width = "18px"), (p.style.textAlign = "center"), (p.style.display = "inline-block"), s.appendChild(p);
                                    }
                                    return s.appendChild(r), s;
                                })(e.hsp[l].subject, e.hsp[l].subjectStart, e.hsp[l].subjectEnd);
                            if (
                                ((r.style.color = "#2c3e50"),
                                (r.style.paddingTop = "25px"),
                                (r.style.paddingBottom = "40px"),
                                (r.style.textAlign = "left"),
                                (r.style.fontFamily = "Helvetica,Arial,sans-serif"),
                                (r.id = "blast-single-alignment-pre"),
                                (s.style.margin = "0 auto"),
                                (s.style.display = "table"),
                                (s.style.paddingTop = "30px"),
                                (i.style.textAlign = "right"),
                                (i.style.paddingTop = "50px"),
                                (p.className = "btn"),
                                (p.innerHTML = "Download as PNG"),
                                (p.style.marginRight = "10px"),
                                (p.onclick = function () {
                                    x(e, "PNG");
                                }),
                                (a.className = "btn"),
                                (a.innerHTML = "Download as JPEG"),
                                (a.onclick = function () {
                                    x(e, "JPEG");
                                }),
                                (o.innerHTML = e.description),
                                (o.style.paddingLeft = "15px"),
                                (o.style.fontWeight = "bold"),
                                (o.style.fontSize = "15px"),
                                (o.style.fontFamily = "Helvetica,Arial,sans-serif"),
                                (d.style.paddingTop = "20px"),
                                (d.style.fontSize = "14px"),
                                (d.style.textAlign = "center"),
                                (d.style.fontFamily = "Helvetica,Arial,sans-serif"),
                                (d.style.display = "table"),
                                (d.style.margin = "0px auto"),
                                (d.style.width = "100%"),
                                (c.innerHTML = "<b>Score:</b></br>" + e.hsp[l].score),
                                (c.style.float = "left"),
                                (c.style.width = "20%"),
                                (h.innerHTML = "<b>Expect:</b></br>" + e.hsp[l].eValue),
                                (h.style.float = "left"),
                                (h.style.width = "20%"),
                                (u.innerHTML = "<b>Identities:</b></br>" + e.hsp[l].identities + "%"),
                                (u.style.float = "left"),
                                (u.style.width = "20%"),
                                (f.innerHTML = "<b>Positives:</b></br>" + e.hsp[l].positives + "%"),
                                (f.style.float = "left"),
                                (f.style.width = "20%"),
                                (g.innerHTML = "<b>Gaps:</b></br>" + e.hsp[l].gaps + "%"),
                                (g.style.float = "left"),
                                (g.style.width = "20%"),
                                (m.style.clear = "both"),
                                d.appendChild(c),
                                d.appendChild(h),
                                d.appendChild(u),
                                d.appendChild(f),
                                d.appendChild(g),
                                d.appendChild(m),
                                s.appendChild(y),
                                s.appendChild($),
                                s.appendChild(v),
                                r.appendChild(o),
                                r.appendChild(d),
                                r.appendChild(s),
                                i.appendChild(r),
                                1 < n)
                            ) {
                                var b = document.createElement("button");
                                if (((b.className = "btn"), (b.id = "blast-single-alignment-next"), (b.innerHTML = "Next HSP"), (b.style.marginTop = "5px"), (b.style.marginRight = "15px"), (b.style.float = "right"), l == n - 1)) var C = 0;
                                else C = l + 1;
                                (b.onclick = function () {
                                    _(e, t, n, C);
                                }),
                                    r.appendChild(b);
                            }
                            i.appendChild(p), i.appendChild(a);
                        }
                        function x(e, t) {
                            var n = document.getElementById("blast-single-alignment-pre"),
                                l = document.getElementById("blast-single-alignment-next");
                            void 0 !== l && null != l && n.removeChild(l),
                                html2canvas(n, {
                                    onrendered: function (i) {
                                        document.body.appendChild(i);
                                        var r = document.createElement("a");
                                        document.body.appendChild(r),
                                            (r.href = "JPEG" == t ? i.toDataURL("image/jpeg") : i.toDataURL("img/png")),
                                            (r.download = e.description + "-alignment." + t.toLowerCase()),
                                            r.click(),
                                            document.body.removeChild(i),
                                            document.body.removeChild(r),
                                            void 0 !== l && null != l && n.appendChild(l),
                                            (r = document.createElement("a")),
                                            document.body.appendChild(r),
                                            (r.href = "#blast-single-alignment"),
                                            r.click(),
                                            document.body.removeChild(r);
                                    },
                                });
                        }
                        function T(e) {
                            switch (e) {
                                case "A":
                                    return "#DBFA60";
                                case "C":
                                    return "#F9FA60";
                                case "D":
                                    return "#F9605F";
                                case "E":
                                    return "#F9609C";
                                case "F":
                                    return "#5FF99D";
                                case "G":
                                    return "#F9BC5F";
                                case "H":
                                    return "#609DF9";
                                case "I":
                                    return "#99F95A";
                                case "K":
                                    return "#A062FF";
                                case "L":
                                    return "#7EF960";
                                case "M":
                                    return "#63FF63";
                                case "N":
                                    return "#D95DF9";
                                case "P":
                                    return "#F9DA60";
                                case "Q":
                                    return "#F955D8";
                                case "R":
                                    return "#5360FB";
                                case "S":
                                    return "#F97E60";
                                case "T":
                                    return "#FFA563";
                                case "V":
                                    return "#C0F86B";
                                case "W":
                                    return "#FDD9F9";
                                case "Y":
                                    return "#60F9DA";
                                default:
                                    return "#FFFFFF";
                            }
                        }
                        String.prototype.startsWith = function (e) {
                            return 0 === this.indexOf(e);
                        };
                    },
                });
            },
            { "js-class": 1, "js-extend": 2 },
        ],
    },
    {},
    []
);
