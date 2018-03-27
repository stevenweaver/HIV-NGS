var coverage_data;
var positional_data;
var rate_json;
var summary_json;
var intrahost_data = null;
var available_runs;
var available_trees;
var rendered_tree;
var rendered_tree_svg;
var rendered_tree_title;
var colors = d3.scale.category10();
var color_rates = d3.interpolateRgb("#AAAAAA", "#000000");
var color_dram = d3.interpolateRgb("#FFB6B6", "#FF0000");
var color_pvalue = d3.interpolateRgb("#FF0000", "#0000FF");
var number_format = d3.format(".2g"),
    compact_format = d3.format(".2f"),
    percentage_format = d3.format(".4g");
var popover_obj = null,
    popover_codon;
var available_run_data;
var relative_path_prefix = "";
var directory_to_data = {};
var parse_date = d3.time.format("%Y/%m/%d");

var data_for_global_filter = null;
var column_names_for_global_filter = null;
var function_for_global_filter;
var view_toggle_list = ['summary_table', 'intrahost_table', 'allofit', 'dram_summary_table', 'compartmentalization_table', 'overall_table'];
var has_compartment = false;
var has_replicate = false;
var use_replicate_for_intrahost = true;
var fst_data = null;
var view_intrahost_list = ['bubble_plot_div', 'tn93_histogram_div'];

var dram_minimum_score = 1;
var show_dram_tab = false;
var compiled_dram_info = [];
var compartment_table = [];
var dram_extractor_list = null;
var summary_table_headers;
var summary_table_annotations;
var additional_filtering_callback = null;
var intrahost_column_names = null;
var dram_column_names;
var unique_file_names = [];
var visible_run_data = [];



var intrahost_table_key_map = [
    ['total_diversity', 100],
    ['s_diversity', 100],
    ['ns_diversity', 100],
    ['tn93_diversity', 100],
    ['total_divergence', 100],
    ['s_divergence', 100],
    ['ns_divergence', 100],
    ['tn93_divergence', 100, {
        'ref': 'tn93_divergence_histogram',
        'key': 'tn93_histogram',
        'action': function(d, label) {
            render_histogram({
                "Histogram": d['tn93_histogram']
            }, [500, 450], 'tn93_histogram_div', ['TN93 distance', 'Count'], label, 14);
            d3.select("#intrahost_table_div").style("height", "400px").style("overflow", "scroll");
            toggle_view("tn93_histogram_div", view_intrahost_list);
            return false;
        }
    }],
    ['Length', 1]
];

var intrahost_table_columns = ["PID",
    "Gene",
    "Sample date",
    "Diversity, %",
    "Syn diversity, %",
    "Non-syn diversity, %",
    "TN93 diversity, %",
    "Divergence, %",
    "Syn divergence, %",
    "Non-syn divergence, %",
    "TN93 divergence, %",
    "Mean length, aa",
    "Graphics"
];

var global_filter_value = [];

var global_summary_filter = function(d, i) {
    if (additional_filtering_callback && !additional_filtering_callback(d)) {
        return false;
    }

    if (global_filter_value.length == 0) {
        return true;
    }

    var matches = [];
    for (var v = 0; v < global_filter_value.length; v++) {
        if (_.isString(global_filter_value[v]) && global_filter_value[v].length > 2 && global_filter_value[v][0] == "/" && global_filter_value[v][global_filter_value[v].length - 1] == "/") {
            matches.push ([RegExp(global_filter_value[v].substr(1, global_filter_value[v].length - 2),'g'), true])
        } else {
            matches.push([global_filter_value[v],false]);
        }
    }

    for (var k = 0; k < d.length && matches.length > 0; k++) {

        var match_me = null;

        if (typeof d[k] == "string") {
            match_me = d[k];
        } else {
            if (typeof d[k] == "boolean") {
                match_me = (d[k]) ? "Yes" : "No";
            } else {
                if (d[k] == null) {
                    match_me = "N/A";
                } else {
                    if (d[k] && typeof d[k] == 'object' && ('text' in d[k])) {
                        match_me = d[k]['text'];
                    } else {
                        if (typeof d[k] == "number") {
                            match_me = "" + d[k];
                        }
                    }
                }
            }
        }
        if (match_me !== null) {
            for (v = matches.length - 1; v >= 0; v--) {
                if (_.isArray(matches[v][0])) {
                    if (k == matches[v][0][0] && matches[v][0][1](+match_me)) {
                        matches.splice(v, 1);
                    }
                } else {
                    if (matches[v][1]) {
                        if (matches[v][0].exec(match_me)) {
                            matches.splice(v, 1);
                        }
                    } else {
                        if (match_me.indexOf(matches[v][0]) >= 0) {
                            matches.splice(v, 1);
                        }
                    }
                }
            }
        }
    }
    return matches.length == 0;
};


function check_compartmentalization(pid, date, gene, compartment) {
    if (fst_data && (pid in fst_data) && (date in fst_data[pid]) && (gene in fst_data[pid][date])) {
        return fst_data[pid][date][gene].filter(function(d) {
            return d[2] && compartment in d[2];
        });
    }
    return null;
}

function correct_index_main_table (index) {
    if (!has_compartment) index--;
    if (!has_replicate) index--;
    return index;
}

function load_analysis_results(dir_info, document_title, ignore_replicate_for_intrahost) {

    d3.select("#document_title").text(document_title);

    if (ignore_replicate_for_intrahost) {
        use_replicate_for_intrahost = false;
    }

    loadDRM([
        ['//spond.github.io/layout-phylotree/data/DRM/Scores_PI.txt', 'PR', 'PI'],
        ['//spond.github.io/layout-phylotree/data/DRM/Scores_NRTI.txt', 'RT', 'NRTI'],
        ['//spond.github.io/layout-phylotree/data/DRM/Scores_NNRTI.txt', 'RT', 'NNRTI']
    ], dir_info);

}

function diff_host_gene(r1, r2, cols) {
    cols = cols || 2;
    for (k = 0; k < cols; k++) {
        if (r1[k] != r2[k]) {
            return true;
        }
    }
    return false;
}

function render_compartment_table(the_data) {
    d3.select("#compartmentalization_table_body").selectAll("tr").remove();
    var rows = d3.select("#compartmentalization_table_body").selectAll("tr").data(the_data.filter(global_summary_filter),
        function(d) {
            return d;
        }
    );
    rows.enter().append('tr');
    rows.exit().remove();

    var csv_table_headers = ['ID', 'date', 'gene', 'Compartment 1', 'Compartment 2', 'Replicate 1', 'Replicate 2', 'F_st', 'p'];
    var csv_table_data = [];
    var current_id;


    d3.select("#summary_matching").text(rows[0].length);
    rows.selectAll("td").data(function(d) {
        return d;
    }).enter().
    append('td').
    each(function(id_data) {
        var table_cell = d3.select(this);
        if (_.isString(id_data)) {
            table_cell.text(id_data);
            current_id = id_data;
        } else { // object
            paired_comparisons = [];

            id_data.comparison.forEach(function(d, i, array) {
                csv_table_data.push([current_id, d[2], d[3], d[0], d[1], d[6], d[7], d[4], d[5]]);
            });
            id_data.compartments.forEach(function(c, i) {
                id_data.compartments.forEach(function(c2, i2) {
                    if (i2 > i) {
                        paired_comparisons.push([
                            c,
                            c2,
                            id_data.comparison.filter(function(row) {
                                return row[0] == c && row[1] == c2;
                            })
                        ]);
                    }
                })
            });

            var charts = table_cell.selectAll("svg").data(paired_comparisons);
            charts.enter().append("svg");
            charts.exit().remove();
            charts.attr("width", "300").attr("height", "200").each(function(data) {
                render_fst(data, this, [300, 200], 16);
            });

        }
    });



    prepare_csv_export([csv_table_headers, csv_table_data, 'compartments.txt'], true);
}

function draw_intrahost_table(the_data) {
    var col_count = intrahost_table_columns.length;
    //d3.select("#intrahost_table_body").selectAll("tr").remove();

    var rows = _set_rows("intrahost_table_body", the_data.filter(global_summary_filter), function(d) {
        return d;
    });

    d3.select("#summary_matching").text(rows[0].length);
    var local_data = rows.data();

    rows.selectAll("td").data(function(d) {
        return d;
    }).enter().
    append('td').
    style("opacity", function(d, i, j) {
        var this_obj = d3.select(this);
        if (i < col_count - 1) {
            if (!d) {
                this_obj.text('N/A');
            } else {
                if (typeof d == "number") {
                    this_obj.text(compact_format(d));
                } else {
                    if (typeof d == "boolean") {
                        this_obj.text(d ? "Yes" : "No");
                    } else {
                        if (typeof d == "object") {
                            if ("text" in d) {
                                if (typeof d["text"] == "object") {
                                    this_obj.html(d["text"].join("<br>"));
                                } else {
                                    this_obj.text(d["text"]);
                                }



                                //if (! ('skip_button'  in d)) {
                                var annotation_div = this_obj.append("div").attr("class", "dropdown pull-right dropdown-menu-right").style('display', 'inline-block');
                                var annotation_button = annotation_div.append("button").attr("class", "pull-right btn btn-primary btn-xs");

                                if ("data" in d) {
                                    if (d['data'].length > 0) {
                                        if (d['data'].length == 1) {
                                            annotation_button.on("click", function(d) {
                                                render_histogram(d['data'][0][2], [900, 300], 'compartment_plot_div', ['TN93 distance', 'Pairs'], d["label"]);
                                                d3.select("#intrahost_table_div").style("height", "400px").style("overflow", "scroll");
                                                toggle_view("compartment_plot_div", view_intrahost_list);
                                                return false;
                                            });
                                        } else {
                                            annotation_button.classed('dropdown-toggle', true).attr('data-toggle', 'dropdown');
                                            var annotation_dropdown = annotation_div.append("ul").attr("class", "dropdown-menu").attr("role", "menu");
                                            var menu_items = annotation_dropdown.selectAll("li").data(d['data']);
                                            menu_items.enter().append("li");
                                            menu_items.attr("role", "presentation");
                                            menu_items = menu_items.selectAll("a").data(function(d) {
                                                return [d]
                                            });
                                            menu_items.enter().append("a");
                                            menu_items.text(function(d) {
                                                return d[5];
                                            }).
                                            on("click", function(d) {
                                                render_histogram(d[2], [900, 300], 'compartment_plot_div', ['TN93 distance', 'Pairs'], d[4]);
                                                d3.select("#intrahost_table_div").style("height", "400px").style("overflow", "scroll");
                                                toggle_view("compartment_plot_div", view_intrahost_list);
                                                return false;
                                            });
                                        }
                                    }
                                } else {
                                    if ("action" in d) {
                                        annotation_button.on('click', function(d) {
                                            return d["action"](d, d["label"]);
                                        });
                                    }
                                }
                                annotation_button.append("span")
                                    .attr("class", 'glyphicon glyphicon-signal');
                            }
                        } else {
                            this_obj.html(d);
                        }
                    }
                }
            }
        } else {
            if (d) {
                if (d3.keys(d).length > 6 - (has_compartment ? 0 : 1) - (has_replicate && use_replicate_for_intrahost ? 0 : 1)) {
                    this_obj.append("button").attr("class", "btn btn-primary btn-xs").text("Display").on('click',
                        function(d, i, j) {
                            plot_div_data(d, ["total_divergence", "total_diversity",
                                    "ds_divergence", "dn_divergence", "tn93_divergence"
                                ], ["Total divergence", "Diversity", "dS", "dN", "TN93 pairwise"],
                                "bubble_plot", [900, 300]);
                            d3.select("#intrahost_table_div").style("height", "400px").style("overflow", "scroll");
                            toggle_view("bubble_plot_div", view_intrahost_list);
                        }
                    );
                }
            }
        }


        if (i < 2) {
            if (j > 0) {
                for (k = 0; k <= i; k++) {
                    if (local_data[j - 1][k] != local_data[j][k]) {
                        return 1;
                    }
                }
                return 0.25;
            }
        }
        return 1.0;
    });

    prepare_csv_export([intrahost_column_names, rows, "intrahost.txt"]);

}

function recurse_json(json, accumulator, record, current_depth, max_depth, last) {
    if (current_depth < max_depth) {
        for (var key in json) {

            var new_record = record.map(function(d) {
                return d;
            });
            new_record.push(key);

            recurse_json(json[key], accumulator, new_record, current_depth + 1, max_depth, json);
        }
    } else {
        //console.log (record);
        var new_record = record.map(function(d) {
            return d;
        });
        for (var c in intrahost_table_key_map) {
            var mapper = intrahost_table_key_map[c];
            var dp = json[mapper[0]];
            if (dp != null) {
                new_record.push(dp * mapper[1]);
                if (mapper.length == 3 && typeof mapper[2] === 'object') {
                    if (mapper[2]['ref'] in json) {
                        new_record[new_record.length - 1] = {
                            "text": compact_format(new_record[new_record.length - 1])
                        };
                        new_record[new_record.length - 1][mapper[2]['key']] = json[mapper[2]['ref']];
                        new_record[new_record.length - 1]['action'] = mapper[2]['action'];
                    }
                }

            } else {
                new_record.push(null);
            }
        }
        new_record.push(last);
        accumulator.push(new_record);
    }
}

function make_patient_label(row, compartments, replicates, extra) {
    var id = 0;
    label = row[id++] + " (" + row[id++] + ")";
    if (compartments) {
        label += " " + row[id++];
    }
    if (replicates) {
        label += ", replicate " + row[id++] + ",";
    }
    if (extra) {
        label += extra;
    }
    label += " sampled on " + row[id] + ". ";
    return label;
}

function make_patient_label_consensus(row, compartments, replicates) {
    var id = 0;
    label = row[id++] + "_" + row[id++].replace(/\//g, "") + "_" + row[id];
    if (compartments) {
        id++;
        label += "_" + row[id++];
    }
    if (replicates) {
        label += "_" + row[id++];
    }
    return label + "_consensus";
}


function make_intrahost_table(json, compartments, replicates) {
    if (!json) {
        intrahost_column_names = ["No intrahost evolution data found in the NGS report"];
    } else {
        intrahost_column_names = intrahost_table_columns;
    }


    intrahost_column_names = _make_table_headers("intrahost_table_head", intrahost_column_names);

    intrahost_data = [];

    var sortable_columns = 4 + (compartments ? 1 : 0) + (replicates && use_replicate_for_intrahost ? 1 : 0);


    recurse_json(json, intrahost_data, [], 1, sortable_columns);

    intrahost_data.sort(function(a, b) {
        for (var k = 0; k < sortable_columns; k++) {
            if (a[k] < b[k]) {
                return -1;
            }
            if (a[k] > b[k]) {
                return 1;
            }
        }
        return 0;
    });



    for (r = 0; r < intrahost_data.length; r++) {
        var lbl = make_patient_label(intrahost_data[r], compartments, replicates && use_replicate_for_intrahost);
        if (r == 0 || diff_host_gene(intrahost_data[r - 1], intrahost_data[r], sortable_columns - 2)) {
            //intrahost_data[r].push (json[intrahost_data[r][0]][intrahost_data[r][1]]);
            var id = 0;
            intrahost_data[r][intrahost_data[r].length - 1]["ID"] = intrahost_data[r][id++];
            intrahost_data[r][intrahost_data[r].length - 1]["gene"] = intrahost_data[r][id++];
            if (compartments) {
                intrahost_data[r][intrahost_data[r].length - 1]["compartment"] = intrahost_data[r][id++];
            }
            if (replicates && use_replicate_for_intrahost) {
                intrahost_data[r][intrahost_data[r].length - 1]["replicate"] = intrahost_data[r][id++];
            }
        } else {
            intrahost_data[r].pop();
            intrahost_data[r].push(null);
        }


        for (c = 0; c < intrahost_data[r].length; c++) {
            if (intrahost_data[r][c] && typeof intrahost_data[r][c] === 'object') {
                intrahost_data[r][c]['label'] = lbl;
            }
        }

        if (0 && fst_data) {
            var my_fst = check_compartmentalization(intrahost_data[r][0], intrahost_data[r][sortable_columns - 2], intrahost_data[r][1], intrahost_data[r][2]);

            if (r == 0) {
                console.log(fst_data[intrahost_data[r][0]][intrahost_data[r][sortable_columns - 2]][intrahost_data[r][1]]);
            }
            continue;

            if (my_fst && my_fst.length) {

                if (my_fst[0][2] && !('p' in my_fst[0][2])) {
                    my_fst = null;
                } else {
                    var list_of_compartments = [];

                    var my_compartment = intrahost_data[r][2];
                    var my_replicate = replicates ? intrahost_data[r][3] : 1;


                    var exclude = {
                            'f_st': 1,
                            'p': 1,
                            'Between': 1
                        },
                        exclude_all = {
                            'f_st': 1,
                            'p': 1,
                            'Between': 1
                        };

                    exclude[intrahost_data[r][2]] = 1;

                    for (pair = 0; pair < my_fst.length; pair++) {
                        var my_compartments = [my_compartment, ''];
                        var my_label = " vs ";
                        var my_menu_item = my_compartment + "(" + my_replicate + ") vs ";
                        for (k in my_fst[pair][2]) {

                            if (!(k in exclude)) {
                                list_of_compartments.push(k);
                                exclude[k] = 1;
                            }
                            if (!(k in exclude_all)) {
                                if (k != my_compartment) {
                                    my_compartments[1] = k;
                                    var other_replicate = my_fst[pair][0] == my_replicate ? my_fst[pair][1] : my_fst[pair][0];
                                    my_label += k + " ( replicate " + other_replicate + ")";
                                    my_menu_item += k + " (" + other_replicate + ")";

                                }
                            }
                        }
                        my_fst[pair].push(my_compartments);
                        my_fst[pair].push(make_patient_label(intrahost_data[r], compartments, replicates, my_label));
                        my_fst[pair].push(my_menu_item);
                    }


                    my_fst = {
                        "data": my_fst,
                        "text": []
                    };
                    //console.log (my_fst, list_of_compartments);

                    for (k = 0; k < list_of_compartments.length; k++) {
                        var only_this_compartment = my_fst["data"].filter(function(d) {
                            return list_of_compartments[k] in d[2];
                        });

                        var comp_count = only_this_compartment.filter(function(d) {
                                return d[2]["p"] <= 0.05;
                            }),
                            comparisons_done = only_this_compartment.length,
                            frac = comp_count.length / comparisons_done;

                        var status = list_of_compartments[k] + " : ";

                        if (frac > 0.75) {
                            status += ("Yes");
                        } else if (frac > 0.4) {
                            status += ("Mixed");
                        } else {
                            status += ("No");
                        }

                        if (comparisons_done > 1) {
                            status += " (" + comp_count.length + "/" + comparisons_done + ")";
                        }

                        my_fst["text"].push(status);

                    }

                    if (my_fst["text"].length > 1) {
                        my_fst['skip_button'] = true;
                        //console.log (my_fst, list_of_compartments);
                    }

                    my_fst["label"] = lbl;
                }


            }
            if (r == 0) {
                console.log(intrahost_data[r], my_fst);
            }
            //intrahost_data[r].splice (intrahost_data[r].length-1, 0, my_fst);
        }

        //console.log (intrahost_data[r].length, intrahost_table_columns.length);

    }

    if (fst_data) {
        /*
          IN:
            id
            date
             gene
              array
                [replicate, replicate, {
                                            "comp1",
                                            "comp2",
                                            "Between",
                                            "f_st",
                                            "p"
                                        }
                ]
        */

        /* OUT:
            id: {
                    compartments: [comp1, comp2, ...]
                    dates:        [date1, date2, ...]
                    genes:        [gene1, gene2, ...]
                    comparisons:  [comp1, comp2, date, gene, f_st, p, replicate1, replicate2, histogram_data (object with 3 arrays)]
                }

        */

        rebinned = {};

        _.each(fst_data, function(id_level, id) {
            rebinned[id] = {
                compartments: {},
                dates: {},
                genes: {},
                comparison: []
            };
            _.each(id_level, function(date_level, date) {
                rebinned[id].dates[date] = true;
                _.each(date_level, function(gene_level, gene) {
                    rebinned[id].genes[gene] = true;
                    _.each(gene_level, function(fst_level) {
                        compartment_keys = _.difference(_.keys(fst_level[2]), ["Between", "f_st", "p"]);
                        if (compartment_keys.length == 2) {
                            if (compartment_keys[0] < compartment_keys[1]) {
                                c1 = compartment_keys[0];
                                c2 = compartment_keys[1];
                            } else {
                                c1 = compartment_keys[1];
                                c2 = compartment_keys[0];
                            }
                            rebinned[id].compartments[c1] = true;
                            rebinned[id].compartments[c2] = true;

                            hist = {
                                "Between": fst_level[2]["Between"],
                                "f_st": fst_level[2]["f_st"],
                                "p": fst_level[2]["p"],
                                "id": id
                            };
                            hist[c1] = fst_level[2][c1];
                            hist[c2] = fst_level[2][c2];

                            rebinned[id].comparison.push([c1, c2, date, gene,
                                fst_level[2]["f_st"], fst_level[2]["p"],
                                replicates ? fst_level[0] : null,
                                replicates ? fst_level[1] : null,
                                hist
                            ]);
                        }
                    })
                })
            })
        });

        _.each(rebinned, function(value, key) {
            ["compartments", "genes", "dates"].forEach(function(array_key) {
                value[array_key] = _.sortBy(_.keys(value[array_key]), function(d) {
                    return d;
                });
            });
        });
        compartment_table = _.pairs(rebinned);


    }
    draw_intrahost_table(intrahost_data);
}


function set_button_handlers() {
    $('body').on('click', '#back_to_summary', function(event) {
        toggle_view("summary_table");
    });

    $('#close_bubble_plot').on('click', function(event) {
        $('#bubble_plot_div').hide(0);
        d3.select("#intrahost_table_div").style("height", null).style("overflow", null);

    });

    $('#export_consensus_sequences').on('click', function(event) {

        d3.select("#export_consensus_sequences_div").remove();

        var progress = d3.select("body").insert("div", ":first-child")
            .attr("id", "export_consensus_sequences_div")
            .attr("class", "alert alert-info alert-dismissable");

        progress.append("button")
            .attr("type", "button")
            .attr("class", "close")
            .attr("data-dismiss", "alert")
            .attr("aria-hidden", "true")
            .html("&times;");

        progress.append("h5")
            .text("Preparing data for export");

        /*loadConsensusSequences (unique_file_names.map (function (d) {return d;}), sequences);*/

        var consensus_index = correct_index_main_table (5);

        var export_sequences = visible_run_data.filter(function(d) {
               return d[consensus_index] && 'consensus' in d[consensus_index] && d[consensus_index]['consensus'];
            }).map(function(d) {
                return [make_patient_label_consensus(d, has_compartment, has_replicate), d[consensus_index]['consensus']];
            }),
            s_count = export_sequences.length;

        export_sequences = encodeURIComponent(export_sequences.map(function(d) {
            return ">" + d[0] + "\n" + d[1] + "\n";
        }).join("\n"));

        d3.select("#export_consensus_sequences_div").select("h5").remove();

        d3.select("#export_consensus_sequences_div").append("a").classed("btn btn-primary", true)
            .text("Download " + s_count + " consensus sequences as FASTA")
            .attr('download', "consensus.txt")
            .attr('href', 'data:text/plain;charset=utf-8,' + export_sequences)
            .on("click", function(d) {
                d3.select("#export_consensus_sequences_div").remove();
            });

        event.preventDefault();
    });


    $('#close_tn93_histogram').on('click', function(event) {
        $('#tn93_histogram_div').hide(0);
        d3.select("#intrahost_table_div").style("height", null).style("overflow", null);

    });

    $('#close_run_metric_histogram').on('click', function(event) {
        $('#run_metric_histogram_div').hide(0);
    });

    $('#close_compartment_plot').on('click', function(event) {
        $('#compartment_plot_div').hide(0);
        d3.selectAll(".view-compartmentalization-sel").classed('view-compartmentalization-sel', false);
        d3.select("#intrahost_table_div").style("height", null).style("overflow", null);

    });

    function toggle_tab_highlight(id) {
        ["#view_individual_runs", "#view_intrahost", "#view_dram_table", "#view_compartmentalization_table"].forEach(
            function(d) {
                d3.select(d).attr("class", d == id ? "alert-info" : null);
            });
    }

    $('#view_individual_runs').on('click', function(event) {
        toggle_tab_highlight("#view_individual_runs");
        data_for_global_filter = available_run_data;
        column_names_for_global_filter = summary_table_headers;
        function_for_global_filter = make_summary_table;
        function_for_global_filter(data_for_global_filter);
        toggle_view("summary_table");
    });



    $('#view_overall_stats').on('click', function(event) {
        toggle_tab_highlight("#view_overall_stats");
        data_for_global_filter = available_run_data;
        column_names_for_global_filter = summary_table_headers;
        function_for_global_filter = make_overall_table;
        function_for_global_filter(data_for_global_filter);
        toggle_view("overall_table");
    });

    $('#view_intrahost').on('click', function(event) {
        toggle_tab_highlight("#view_intrahost");
        column_names_for_global_filter = intrahost_table_columns;
        function_for_global_filter = draw_intrahost_table;
        data_for_global_filter = intrahost_data;
        function_for_global_filter(data_for_global_filter);
        toggle_view("intrahost_table");
    });

    $('#view_dram_table').on('click', function(event) {
        toggle_tab_highlight("#view_dram_table");
        column_names_for_global_filter = null;
        function_for_global_filter = render_dram_table_wrapper;
        data_for_global_filter = compiled_dram_info;
        function_for_global_filter(data_for_global_filter);
        toggle_view("dram_summary_table");
    });

    $('#view_compartmentalization_table').on('click', function(event) {
        toggle_tab_highlight("#view_compartmentalization_table");
        column_names_for_global_filter = null;
        function_for_global_filter = render_compartment_table;
        data_for_global_filter = compartment_table;
        function_for_global_filter(data_for_global_filter);
        toggle_view("compartmentalization_table");
    });


    $('#summary_limiter').on('input propertychange', _.throttle(function(event) {
        global_filter_value = $(event.target).val().split(" ");

        if (column_names_for_global_filter) {
            var expression_splitter = /^([^><]+)([>=<])([^\s]+)$/;

            for (var v = 0; v < global_filter_value.length; v++) {
                var iz_column = global_filter_value[v].match(expression_splitter);
                if (iz_column) {
                    //console.log (iz_column);
                    for (c = 0; c < column_names_for_global_filter.length; c++) {
                        if (column_names_for_global_filter[c].indexOf(iz_column[1]) == 0) {
                            break;
                        }
                    }
                    if (c < column_names_for_global_filter.length) {
                        var l = +iz_column[3];
                        global_filter_value[v] = [c, iz_column[2] == '<' ? function(v) {
                            return v < l;
                        } : (iz_column[2] == '>' ? function(v) {
                            return v > l;
                        } : function(v) {
                            return v == l;
                        })];
                        continue;
                    }
                }
            }
        }
        function_for_global_filter(data_for_global_filter);
    }, 250));
}


var DRAMs = {};


function numericDRAM(d) {
    for (var k in d) {
        if (k != 'AA') {
            d[k] = +d[k];
        }
    }
    return d;
}

function show_summary_popup_text(k, v) {
    if (typeof v == "string" || typeof v == "number") {
        return "<li><strong>" + k + "</strong> " + v + "</li>";
    }
    if (typeof v == "object") {
        var my_list = [];
        var ordered_keys = d3.keys(v);
        ordered_keys.sort();
        ordered_keys.forEach(function(d) {
            my_list.push(show_summary_popup_text(d, v[d]));
        });
        return ((k ? ("<li><strong>" + k + "</strong></li>") : "")) + "<ul>" + my_list.join("\n") + "</ul>";
    }
}

function show_summary_popup(d) {
    var jq_obj = $(this);


    if (popover_obj) {
        $(popover_obj).popover('destroy');
        if (this == popover_obj) {
            popover_obj = null;
            return;
        }
    }
    jq_obj.popover({
        animation: false,
        placement: "bottom",
        title: "Run metrics",
        html: true,
        content: show_summary_popup_text(null, d["pop-up"]),
        trigger: "click"
    });
    jq_obj.popover('show');
    popover_obj = this;
}

function main_loader(dir_info) {
    d3.json(dir_info[0], function(error, json) {
        available_run_data = json;

        if ('settings' in json) {
            has_compartment = json['settings']['compartment'];
            has_replicate = json['settings']['replicate'];
        }

        relative_path_prefix = dir_info[1];

        var dram_info_from = [];
        for (k = 0; k < json['data'].length; k++) {
            var this_record = json['data'][k];
            this_record[this_record.length - 1] = relative_path_prefix + this_record[this_record.length - 1];
            directory_to_data[this_record[this_record.length - 1]] = this_record;
            unique_file_names.push(this_record[this_record.length - 1]);
            var dram_data_available = has_dram_relevant_genes(this_record[2]);
            if (dram_data_available >= 0) {
                dram_info_from.push([this_record[this_record.length - 1], dram_data_available]);
            }
        }

        summary_table_headers = json['columns'];
        summary_table_annotations = 'annotations' in json ? json['annotations'] : summary_table_headers;

        //make_summary_table (json);


        show_dram_tab = dram_info_from.length > 0;

        if (show_dram_tab) {
            dram_extractor_list = extractDRM(null, null, dram_minimum_score);
            d3.select("#dram_info_tab").append("li")
                .append("a")
                .attr('id', 'view_dram_table')
                .attr('href', '#')
                .text('DRAM summary');

            loadIndividualDRAM(dram_info_from);
            d3.select("#dram_button_summary").on("click", function(e) {
                try {
                    var threshold = parseFloat($('#dram_freq_summary').val());
                    additional_filtering_callback = function(d) {
                        return d[4 + (has_replicate ? 1 : 0) + (has_compartment ? 1 : 0)].some(function(m) {
                            return m[2] * 100 >= threshold;
                        });
                    }
                    function_for_global_filter(data_for_global_filter);
                    additional_filtering_callback = null;

                } catch (e) {
                    console.log(e);
                }
            });
        }

        fst_data = json['F_ST'];

        if (has_replicate && use_replicate_for_intrahost) {
            intrahost_table_columns.splice(2, 0, "Replicate");
        }
        if (has_compartment) {
            intrahost_table_columns.splice(2, 0, "Compartment");
        }
        /*if (fst_data) {
            intrahost_table_columns.splice (intrahost_table_columns.length-1, 0, "Compartmentalization");
        }*/

        if ("intrahost" in json) {
            make_intrahost_table(json['intrahost'], has_compartment, has_replicate);
        } else {
            make_intrahost_table(null);
        }

        if (fst_data) {
            d3.select("#dram_info_tab").append("li")
                .append("a")
                .attr('id', 'view_compartmentalization_table')
                .attr('href', '#')
                .text('Compartmentalization');
        }

        set_button_handlers();
        $('#view_intrahost').click();
        $("#loading_bar").hide(0);
        $("#nav_bar").show(0);
    });
};


function loadIndividualDRAM(load_list) {
    if (load_list.length) {
        var dir = load_list.pop();
        return d3.json(dir[0] + '/prot_coverage.json', function(e, d) {
            d3.json(dir[0] + '/rates.json', function(e2, d2) {
                if (d2 && d) {

                    var record_info = directory_to_data[dir[0]].filter(function(d, i) {
                        if (i < 2) {
                            return true
                        };
                        if (i > 2 && i - 2 <= has_compartment + has_replicate) return true;
                        return false
                    });
                    for (k in d) {
                        site_info = site_annotator(k, dir[1]).join("");
                        if (site_info in dram_extractor_list) {
                            var site_counts = get_site_info(d[k]),
                                only_dram = site_counts[2].filter(function(element) {
                                    return element[0] in dram_extractor_list[site_info];
                                });

                            if (only_dram.length) {
                                compiled_dram_info.push([site_info].concat(record_info, site_counts[0], [only_dram]));
                            }
                        }
                    }
                }

                if (load_list.length) {
                    loadIndividualDRAM(load_list);
                } else {
                    compiled_dram_info.sort(function(a, b) {
                        for (k = 1; k < 3; k++) {
                            if (a[k] != b[k]) {
                                return a[k] - b[k];
                            }
                        }
                        if (a[0].substr(0, 2) == b[0].substr(0, 2)) {
                            return parseInt(a[0].substr(3)) - parseInt(b[0].substr(3));
                        }
                        return a[0] - b[0];
                    });
                    render_dram_table('dram_table_body', compiled_dram_info);
                    var headers = ['Site', summary_table_headers[0], summary_table_headers[1]];
                    i = 3;
                    if (has_compartment) {
                        headers.push(summary_table_headers[i++]);
                    }
                    if (has_replicate) {
                        headers.push(summary_table_headers[i++]);
                    }
                    headers.push("Coverage");
                    headers.push("DRAMs");
                    dram_column_names = _make_table_headers("dram_table_head", headers);

                    d3.select("#loading_dram_bar").remove();
                    d3.select("#dram_button_summary").classed("disabled", null);
                }
            });
        });
    }
}


function loadDRM(load_list, dir_name) {

    if (load_list.length) {
        var load_this_now = load_list.pop();
        return d3.tsv(load_this_now[0], numericDRAM, function(e, d) {
            if (d) {
                //console.log ("Loading ", load_this_now);
                var drug_names = [];
                for (var name in d[0]) {
                    if (name != 'Position' && name != 'AA') {
                        drug_names.push(name);
                    }
                }
                DRAMs[load_this_now[2]] = {
                    'ARV': drug_names,
                    'positions': []
                };
                for (var drami in d) {
                    var dram = d[drami];
                    var scores = [load_this_now[1] + dram['Position'], dram['AA']];
                    for (var drug in drug_names) {
                        scores.push(dram[drug_names[drug]]);
                    }
                    DRAMs[load_this_now[2]]['positions'].push(scores);
                }
            }

            if (load_list.length) {
                loadDRM(load_list, dir_name);
            } else {
                main_loader(dir_name);
            }
        });
    }
}

function extractDRM(classes, drugs, min_score) {
    var filter_sites = {}; // "site" : list of mutations

    for (var drug_class in DRAMs) {
        if ((!classes) || drug_class in classes) {
            var filter_drugs = [];
            for (var i in DRAMs[drug_class]['ARV']) {
                if (!drugs || DRAMs[drug_class]['ARV'][i] in drugs) {
                    filter_drugs.push(1);
                } else {
                    filter_drugs.push(0);
                }
            }

            for (var s in DRAMs[drug_class]['positions']) {
                var row = DRAMs[drug_class]['positions'][s];
                if (Math.max.apply(null, row.filter(function(d, i) {
                        return (i >= 2 && filter_drugs[i]);
                    })) >= min_score) {
                    if (!(row[0] in filter_sites)) {
                        filter_sites[row[0]] = {};
                    }
                    filter_sites[row[0]][row[1]] = 1;
                }
            }
        }
    }
    return filter_sites;
}

function prepare_csv_export(csv_export_data, already_text) {
    if (csv_export_data.length == 3) {
        if (csv_export_data[0] && csv_export_data[1]) {
            export_data = [
                []
            ];
            if (already_text) {
                export_data[0] = csv_export_data[0];
                csv_export_data[1].forEach(function(d) {
                    export_data.push(d);
                });
            } else {
                csv_export_data[0].each(function(d, i) {
                    export_data[0].push(d3.select(this).text());
                });
                csv_export_data[1].each(function(d, i) {
                    export_data.push([]);
                    d3.select(this).selectAll("td").each(function(cell, j) {
                        export_data[export_data.length - 1].push(d3.select(this).text());
                    });
                });
            }

            d3.select("#export_current_table")
                .attr("href", "data:application/text;charset=utf-8," + encodeURIComponent(export_data.map(function(d) {
                    return d.join("\t");
                }).join("\n")))
                .attr("download", csv_export_data[2]);
        }
    }
}

function _make_table_headers(container, headers, annotations) {
    var rows = d3.select("#" + container).selectAll("tr").data([headers]);
    rows.enter().append('tr');
    rows.exit().remove();
    var cols = rows.selectAll("td").data(function(d) {
        return d;
    });
    cols.enter().append('td');
    cols.exit().remove();

    if (annotations) {
        cols = cols.selectAll("abbr").data(function(d, i) {
            return [{
                "name": d,
                "annotation": summary_table_annotations[i]
            }]
        });
        cols.enter().append("abbr");
        cols.exit().remove();
        cols.text(function(d) {
            return d.name;
        }).attr("title", function(d) {
            return d.annotation;
        });
    } else {
        cols.text(function(d) {
            return d;
        });
    }
    return cols;
}

function has_dram_relevant_genes(gene) {
    return ['pr', 'prrt', 'rt'].indexOf(gene);
}

function make_overall_table(json) {

    /*
    "pop-up": {
     "contributing reads": 21372,
     "q20": 0.927647,
     "mean q-score": 34.4676,
     "q30": 0.768968,
     "total bases": 12341411,
     "original reads": 35124,
     "retained fragments": 21372,
     "retained fragment length distribution:": {
      "median": 396,
      "mean": 334.622,
      "min": 50,
      "2.5%": 95,
      "standard deviation": 91.2781,
      "variance": 8331.69,
      "97.5%": 414,
      "max": 517
     },
     "original read length distribution:": {
      "median": 396,
      "mean": 351.367,
      "min": 28,
      "2.5%": 127,
      "standard deviation": 87.5007,
      "variance": 7656.36,
      "97.5%": 481,
      "max": 552
     },
     "q10": 0.999991
    },
    "text": "4374 (865)"

    "pop-up": {
     "IQR": "4579-4692",
     "mean": "4596.99",
     "range": "311-5199",
     "median": 4606.0
    },


    */

    //

    var column_index = correct_index_main_table(7);
    var ids = new Object,
        extract_keys = [
            ["contributing reads"],
            ["q20"],
            ["total bases"],
            ["retained fragment length distribution:", "median"]
        ];

    function key_extractor(object, key_list, key_index) {
        //console.log (object, key_list);
        var key_index = key_index || 0,
            current_value = object[key_list[key_index]];

        if (key_index < key_list.length - 1) {
            return key_extractor(current_value, key_list, key_index + 1);
        }
        return current_value;
    }


    var clone_metrics = [];
    /*
        mapped reads
        cluster count
        median coverage (per amino-acid base)
    */


    var selected_run_data = json["data"].filter(global_summary_filter).filter(function(d) {
        return d[column_index] !== null;
    }).map(function(d) {
        return [d[column_index], d[column_index - 1]];
    }).filter(function(d) {
        var clones = d[0]["text"].split(" ");
        //var coverages = d[1]["pop-up"]["IQR"].split("-");

        var clusters = +clones[1].substr(1, clones[1].length - 2),
            reads = +clones[0];

        clone_metrics.push([reads, clusters , reads/clusters, +d[1]["text"]]);

        if (!(d[0]["unique_id"] in ids)) {
            ids[d[0]["unique_id"]] = true;
            return true;
        }
        return false;
    }).map(function(d) {
        return extract_keys.map(function(k) {
            return key_extractor(d[0]["pop-up"], k);
        });
    });



    d3.select("#summary_matching").text(selected_run_data.length);

    _make_table_headers("intrument_metrics_table_head", ["Metric", "Mean", "Median", "Interquartile range", "Total", ""]);

    function _populate_table(id, defs, raw_data, totals) {
        var rows = _set_rows(id, defs.map(function(d, i) {
                    var data = raw_data.map(function(r) {
                            return d[1] ? d[1](r[i]) : r[i];
                        }),
                        stats = _describe_a_vector(data);

                    if (totals) {
                        return [d[0], stats[0], stats[1], compact_format(stats[2]) + "-" + compact_format(stats[3]), totals[i] ? stats[5] : "", [d[0],
                            [stats[6],
                                [stats[7], stats[8]]
                            ]
                        ]];
                    } else {
                        return [d[0], stats[0], stats[1], compact_format(stats[2]) + "-" + compact_format(stats[3]), [d[0],
                                [stats[6],
                                    [stats[7], stats[8]]
                                ]
                           ]];

                        }

                    }),
                function(d) {
                    return d[0]
                });


            rows = rows.selectAll("td").data(function(d) {
                return d;
            }); rows.enter().append('td'); rows.style("opacity", function(d, i, j) {
                var this_obj = d3.select(this);
                if (_.isNumber(d)) {
                    this_obj.text(compact_format(d));
                } else {
                    if (_.isArray(d)) {


                        var annotation_div = this_obj.selectAll("div"),
                            annotation_button;
                        if (annotation_div.empty()) {
                            annotation_div = this_obj.append("div").attr("class", "dropdown pull-right dropdown-menu-right").style('display', 'inline-block');
                            annotation_button = annotation_div.append("button").attr("class", "pull-right btn btn-primary btn-xs");
                            annotation_button.append("span")
                                .attr("class", 'glyphicon glyphicon-signal');
                        } else {
                            annotation_button = annotation_div.select("button");
                        }



                        annotation_button.on("click", function(d) {
                            render_histogram(d[1], [900, 300], 'run_metric_histogram_div', [d[0], 'Counts'], "");
                            toggle_view("run_metric_histogram_div", []);
                            return false;
                        });

                    } else {
                        this_obj.text(d);
                    }
                }

                return 1.0;
            }); rows.exit().remove();
        }

        _populate_table("intrument_metrics_table_body", [
            ["Mapped Reads", null],
            ["Bases with q>20, %", function(x) {
                return x * 100;
            }],
            ["Mb", function(x) {
                return x / 1e6;
            }],
            ["Median Read Length", null]
        ], selected_run_data, {
            0: true,
            2: true
        });

        _make_table_headers("run_metrics_table_head", ["Metric", "Mean", "Median", "IQR", ""]);

        _populate_table("run_metrics_table_body", [
                ["Mapped Reads", null],
                ["Clusters", null],
                ["Reads/cluster", null],
                ["Coverage per aminoacid", null]
            ],
            clone_metrics);
    }

    function make_summary_table(json) {

        var col_count = summary_table_headers.length,
            cols = _make_table_headers("summary_table_head", summary_table_headers, summary_table_annotations);

        visible_run_data = json["data"].filter(global_summary_filter);
        rows = d3.select("#summary_table_body").selectAll("tr").data(visible_run_data, function(d) {
            return d[col_count - 1];
        });
        rows.enter().append('tr');
        rows.exit().remove();
        d3.select("#summary_matching").text(rows[0].length);


        rows.selectAll("td").data(function(d) {
            return d;
        }).enter().
        append('td').each(function(d, i) {
            var this_obj = d3.select(this);
            if (i < col_count - 1) {
                if (d == null) {
                    this_obj.text('N/A');
                    return
                };
                if (typeof d == "number") {
                    this_obj.text(compact_format(d));
                    return;
                }
                if (typeof d == "boolean") {
                    this_obj.text(d ? "Yes" : "No");
                    return;
                }
                if (typeof d == "string") {
                    this_obj.text(d);
                    return;
                }
                if (typeof d == "object") {
                    if ('text' in d) {
                        if (typeof d["text"] == "object") {
                            this_obj.append("span").html(d["text"].join("<br>"));
                        } else {
                            this_obj.append("span").text(d["text"]);
                        }
                    }
                    if ('link' in d) {
                        var dl_button = this_obj.insert("a", ":first-child").attr("class", "btn btn-primary btn-xs").style('display', 'inline-block').style('margin-right', '0.2em')
                            .attr("href", relative_path_prefix + d['link']).attr("download", 'target' in d ? d['target'] : 'region.fas')
                            .attr("target", "_blank");


                        dl_button.append("span").attr("class", "glyphicon glyphicon-download");
                    } else {
                        if ('pop-up' in d) {
                            var dl_button = this_obj.insert("a", ":first-child").attr("class", "btn btn-primary btn-xs").style('display', 'inline-block').style('margin-right', '0.2em')
                                .attr("href", '#')
                                .on("click", show_summary_popup);


                            dl_button.insert("span").attr("class", "glyphicon glyphicon-list-alt");
                        }
                    }
                }

            } else {
                this_obj.html('<button type="button" class="btn btn-primary btn-xs">Display</button>');
            }
        }).
        attr("class", function(d, i) {
            if (i == col_count - 1) {
                return "view-analysis";
            }
            return "";
        });

        rows.selectAll(".view-analysis").on('click', function(d, i, j) {
            load_a_directory(d)
        });
        rows.classed('danger', function(d) {
            return d.some(function(c) {
                return c == null;
            }) ? true : null;
        });
        prepare_csv_export([cols, rows, "summary.txt"]);
    }

    function toggle_view(tag, list) {
        list = list || view_toggle_list;
        for (k = 0; k < list.length; k++) {
            if (tag != list[k]) {
                d3.select("#" + list[k]).style('display', 'none');
            }
        }
        d3.select("#" + tag).style('display', 'block');
    }

    function generate_site_info(codon_id) {
        html = "<table class = 'table-striped table table-hover'>\
                    <thead>\
                        <td>Pos</td>\
                        <td>A</td>\
                        <td>C</td>\
                        <td>G</td>\
                        <td>T</td>\
                    </thead>\
                    <tbody>\
                    ";
        nucs = "ACGT";
        for (k = (codon_id - 1) * 3 + 1; k <= (codon_id) * 3; k++) {
            site_info = rate_json['posteriors'][k];
            html += "<tr><td>" + k + "</td>";
            for (ni = 0; ni < 4; ni++) {
                n = nucs[ni];
                html += "<td>" + compact_format(site_info[n][0]) + '<br>' + site_info[n][2] + "</td>";
            }
            html += "</tr>";
        }
        return html + "</tbody></table>";
    }

    function display_site_properties(obj, codon_id) {
        if (popover_obj) {
            popover_obj.popover('destroy');
        }
        if (popover_codon == codon_id) {
            popover_obj = null;
            popover_codon = -1;
            return;
        }
        obj.popover({
            animation: false,
            placement: "left",
            title: "Posterior probability of observing a variant and counts",
            html: true,
            content: generate_site_info(codon_id),
            trigger: "click"
        });
        obj.popover('show');
        popover_obj = obj;
        popover_codon = codon_id;
    }


    function load_a_directory(dir) {


        var info = directory_to_data[dir];
        var has_dram = has_dram_relevant_genes(info[2]);



        clean_existing("tree_plot", true);
        clean_existing("diversity_plot_ns", true);
        clean_existing("diversity_plot_s", true);
        clean_existing("tn93_plot", true);
        clean_existing("rate_table");
        handle_summary_json(info);



        d3.json(dir + '/rates.json', function(error, json) {
            handle_rates_json(json);
            d3.json(dir + '/prot_coverage.json', function(error, json_coverage) {
                plot_coverage_data(json ? json_coverage : null, ['Site, aa', 'Coverage', 'Majority Residue', 'Minority Residues(s)'], 'main_plot', [600, 400], info[2], has_dram);
            });
        });

        d3.json(dir + '/tn93.json', function(error, json) {
            render_histogram(json, [300, 300], 'tn93_plot', ['TN93 distance', 'Count']);
        });



        d3.json(dir + '/diversity.json', function(error, json) {
            render_diversity(json, [300, 300], 'diversity_plot_ns', 'div_span_ns', 'NS', 'Non-synonymous');
            render_diversity(json, [300, 300], 'diversity_plot_s', 'div_span_s', 'S', 'Synonymous');
            render_tree(json, [800, 300], "tree_plot");
        });

        /*
        d3.json (dir + '/fst.json', function (error, json) {
            render_histogram (json, [800,400], 'fst_plot', ['TN93 distance', 'Count']);
        });
        */

        //"PID " + mapper[d][0] + " " + mapper[d][3] + " isolated on " + mapper[d][1] + " from " + mapper[d][2], "data/" + d + "/", mapper[d][3], ['pr','prrt','rt'].indexOf (mapper[d][3]));  }


        d3.select("#analysis_id").text("PID " + info[0] + " " + info[3] + " isolated on " + info[1] + " from " + info[2]);
        toggle_view("allofit");
        d3.select("#dram_box").style('display', has_dram >= 0 ? 'block' : 'none');

    }

    function helper_make_consensus(rate_info, prob, coverage) {
        prob = prob || 0.999;
        coverage = coverage || 1000;

        consensus_sequence = [];
        for (k in rate_info) {
            var site = parseInt(k);
            var obs_coverage = 0,
                maxC = 0,
                maxR = null;

            for (c in {
                    "A": 1,
                    "C": 1,
                    "G": 1,
                    "T": 1
                }) {

                if (rate_info[k][c][0] >= prob) {
                    var lc = rate_info[k][c][2];
                    obs_coverage += lc;
                    if (lc >= maxC) {
                        maxC = lc;
                        maxR = c;
                    }
                }
            }

            if (!(maxR && obs_coverage >= coverage)) {
                maxR = '-';
            }

            consensus_sequence.push([site, maxR]);
        }
        return consensus_sequence.sort(function(a, b) {
            return a[0] - b[0];
        }).map(function(d) {
            return d[1];
        }).join('');

    }

    function handle_rates_json(json) {
        rate_json = json;

        d3.select("#rate_table").selectAll("tr").remove();
        rate_info = []
        if (rate_json) {
            for (k in rate_json["priors"]) {
                v = rate_json["priors"][k];
                rate_info.push([100 * v.weight, v.A, v.C, v.G, v.T]);
            }
            rate_info = rate_info.sort(function(a, b) {
                return b[0] - a[0];
            });
            d3.select("#consensus_text").text(rate_json ? helper_make_consensus(rate_json["posteriors"]) : "No consensus sequence information found (probably due to low coverage/failed run)");
        } else {
            rate_info = [
                ["No rate info"]
            ];
        }
        d3.select("#rate_table").selectAll("tr").data(rate_info).enter().append("tr").
        selectAll("td").data(function(d) {
            return (d);
        }).enter().append("td").
        text(function(d) {
            return typeof d == "number" ? number_format(d) : d;
        }).
        style("color", function(d, i) {
            if (i > 0) {
                return color_rates(d);
            }
            return '';
        });

    }

    function handle_summary_json(info) {

        var offset = 0 + has_replicate + has_compartment;

        mapper = {
            "summary_region_spanned": 3 + offset,
            "summary_coverage": 4 + offset,
            "summary_read_count": 5 + offset,
            "summary_tn93": 6 + offset,
            "summary_S": 7 + offset,
            "summary_NS": 8 + offset,
        };

        for (key in mapper) {
            vals = info[mapper[key]];
            if (!vals) {
                vals = "N/A";
            } else {
                if (typeof vals == "number") {
                    vals = compact_format(vals);
                } else {
                    if (typeof vals == "object") {
                        vals = "text" in vals ? vals["text"] : vals;
                    }

                }
            }
            d3.select("#" + key).text(vals);
        }
    }

    function display_a_particular_tree(split, type) {
        rendered_tree(available_trees[split][type]).layout();
        rendered_tree_title.text(type + " tree for region " + split);
    }

    function render_tree(json, dim, id) {
        var tree_svg;

        clean_existing(id);

        if (json) {
            available_trees = {};

            rendered_tree_title = d3.select("#" + id).append("h5").text("Click on the diversity displays to show trees for particular regions");

            rendered_tree_svg = d3.select("#" + id).append("svg")
                .attr("width", dim[0])
                .attr("height", dim[1]);

            var re = new RegExp("_([0-9]+)$");

            rendered_tree = d3.layout.phylotree()
                .size([dim[1], dim[0]])
                .separation(function(a, b) {
                    return 0;
                })
                .spacing_x(10, true)
                .node_span(function(a) {
                    var m = re.exec(a.name);
                    try {
                        return Math.sqrt(parseFloat(m[1]))
                    } catch (e) {}
                    return null;
                })
                .branch_name(function(a) {
                    var m = re.exec(a.name);
                    try {
                        return m[1];
                    } catch (e) {}
                    return "";
                })
                .svg(rendered_tree_svg)
                .options({
                    'draw-size-bubbles': true
                }, false);


            for (k in json) {
                coords = k.split('-');
                if ('Tree' in json[k]) {
                    available_trees[k] = json[k]['Tree'];
                }
            }


        } else {
            diversity_svg = d3.select("#" + id).append("div").
            text("No phylogenetic tree information found").
            attr("class", "alert  alert-warning");
        }
    }

    function clean_existing(id, set_loading) {
        var element = d3.select("#" + id);
        for (k in {
                "svg": 0,
                "h5": 1,
                "div": 2
            }) {
            element.selectAll(k).remove();
        }
        if (set_loading) {
            element.append("div").attr("class", "progress progress-striped active")
                .append("div")
                .attr("class", "progress-bar")
                .attr("role", "progressbar")
                .attr("aria-valuenow", "100")
                .attr("aria-valuemin", "0")
                .attr("aria-valuemax", "100")
                .style("width", "100%")
                .text("Loading");
        }

    }

    function render_fst(data, svg, dim, sybmol_size) {
        svg = d3.select(svg);
        sybmol_size = sybmol_size || 8;

        var margin = {
                top: sybmol_size + 12,
                right: 2 * sybmol_size,
                bottom: 2 * sybmol_size + 30,
                left: 2 * sybmol_size + 20
            },
            width = dim[0] - margin.left - margin.right,
            height = dim[1] - margin.top - margin.bottom;

        svg.selectAll("path").remove();
        svg.selectAll("g").remove();

        var x, valid_dates = true;

        dates = data[2].map(function(row) {
            d = parse_date.parse(row[2]);
            if (!d) {
                valid_dates = false
            };
            return d;
        });

        if (valid_dates) {
            x = d3.time.scale()
                .rangeRound([0, width]);
            x.domain(d3.extent(dates));
        } else {
            dates = data.map(function(row) {
                return row[2];
            });
            x = d3.scale.ordinal()
                .rangePoints([0, width], 0.5);
            x.domain(_.unique(data[2].map(function(d) {
                return d[2];
            })), false);
        }

        var y = d3.scale.ordinal()
            .rangePoints([height, 0], 0.5);

        y.domain(_.sortBy(data[2].map(function(d) {
            return d[3];
        }, d3.ascending)));

        var circles = [];


        data[2].forEach(function(fst) {
            circles.push({
                gene: fst[3],
                date: valid_dates ? parse_date.parse(fst[2]) : fst[2],
                fst: fst[4],
                p: fst[5],
                replicate1: +fst[6],
                replicate2: +fst[7],
                hist: fst[8]
            });
        });

        function size(fst) {
            fst = fst > 0 ? fst : -fst;
            fst = fst > 1 ? 1 : fst;
            fst = fst < 0.2 ? 0.2 : fst;
            return sybmol_size * sybmol_size * fst;
        }

        // figure title
        svg.append("text").attr("transform", "translate (" + (margin.left + width / 2) + "," + margin.top / 2 + ")")
            .style("text-anchor", "middle")
            .text(data[0] + " vs " + data[1]);

        var count_by_gene_date = []; // deal with replicates

        var viz = svg.append("g").attr("transform", "translate (" + margin.left + "," + margin.top + ")");
        viz.selectAll("path")
            .data(circles)
            .enter().append("path")
            .attr("transform", function(d) {

                var already_did_gene_date = _.findWhere(count_by_gene_date, {
                    'gene': d.gene,
                    'date': parse_date(d.date)
                });

                if (already_did_gene_date) {
                    already_did_gene_date.count += 1;
                } else {
                    already_did_gene_date = {
                        gene: d.gene,
                        date: parse_date(d.date),
                        count: 1
                    };
                    count_by_gene_date.push(already_did_gene_date);
                }

                //console.log (count_by_gene_date);

                return "translate(" + (x(d.date) +
                    (already_did_gene_date.count - 1) * sybmol_size) + "," + y(d.gene) + ")";
            })
            .attr("d", function(d) {
                return d3.svg.symbol().type(d.fst > 0 ? "circle" : "cross").size(size(d.fst))();
            }).
        style("fill", function(d) {
            return color_pvalue(d.p);
        }).
        style("opacity", function(d) {
            return d.p <= 0.05 ? "1.0" : "0.25";
        }).
        on("click", function(d) {
            render_histogram(d.hist, [900, 300], 'compartment_plot_div', ['TN93 distance', 'Pairs'],
                "Patient " + d.hist.id + " sampled on " + parse_date(d.date) + " from " + d.gene +
                (has_replicate ? " replicates " + d.replicate1 + "/" + d.replicate2 : "") + ". "
            );
            d3.select("#intrahost_table_div").style("height", "400px").style("overflow", "scroll");
            toggle_view("compartment_plot_div", view_intrahost_list);
            return false;
        }).
        append("title").text(function(d) {
            return "Fst = " + number_format(d.fst) + ", p = " + compact_format(d.p);
        });

        var xAxis = d3.svg.axis()
            .scale(x)
            .tickValues(dates)
            .tickFormat(parse_date)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");

        var x_axis = viz.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + height + ")")
            .call(xAxis);

        x_axis.selectAll("text").attr("transform", "rotate (30)").
        attr("dx", "3em").
        style("font-size", "6px");

        viz.append("g")
            .attr("class", "y axis")
            .call(yAxis);
        /*.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("gene");*/


    }


    function render_diversity(json, dim, id, class_id, index, label) {
        var diversity_svg;

        clean_existing(id);

        if (json) {
            var data_to_render = [];

            for (k in json) {
                coords = k.split('-');
                data_to_render.push([parseInt(coords[0]), parseInt(coords[1]), json[k]['Diversity'][index]]);
            }

            data_to_render = data_to_render.sort(function(a, b) {
                if (b[1] > a[1] && b[0] < a[1]) return 1;
                if (a[0] < b[0]) return 1;
                if (a[0] > b[0]) return -1;
                return b[1] - a[1];
            })

            var margin = {
                    top: 20,
                    right: 20,
                    bottom: 40,
                    left: 50
                },
                width = dim[0] - margin.left - margin.right,
                height = dim[1] - margin.top - margin.bottom;


            var step = 0;

            var x = d3.scale.linear()
                .domain([d3.min(data_to_render.map(function(d) {
                    return d[0];
                })), d3.max(data_to_render.map(function(d) {
                    return d[1];
                }))])
                .range([0, width]);

            var y = d3.scale.linear()
                .domain([0, d3.max(data_to_render.map(function(d) {
                    return d[2];
                }))])
                .range([height, 0]);


            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");


            d3.select("#" + id).append("h5").
            text(label + " diversity");

            var rate_svg = d3.select("#" + id).append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


            var bar = rate_svg.selectAll("." + class_id)
                .data(data_to_render)
                .enter().append("g")
                .attr("class", class_id);

            bar.append("rect")
                .attr("x", function(d, i) {
                    return x(d[0]);
                })
                .attr("y", function(d, i) {
                    return y(d[2]);
                })
                .attr("width", function(d, i) {
                    return x(d[1]) - x(d[0]);
                })
                .attr("height", function(d) {
                    return y(0) - y(d[2]);
                })
                .on("click", function(d, i, j) {
                    display_a_particular_tree(d[0] + "-" + d[1], index);
                })
                .append("title").text(function(d, i) {
                    return "" + d[1] + " distances between " + (d[0]) + " and " + d[2];
                });


            rate_svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis);

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");

            var yAxis = d3.svg.axis()
                .scale(y)
                .orient("left");

            rate_svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis)
                .append("text")
                .attr("transform", "translate(" + width + ",0)")
                .attr("y", "24")
                .attr("dy", ".71em")
                .style("text-anchor", "end")
                .text("Site, nt");

            rate_svg.append("g")
                .attr("class", "y axis")
                .call(yAxis)
                .append("text")
                //.attr("transform", "rotate(-90)")
                .attr("y", 6)
                .attr("dy", "-1.5em")
                .attr("dx", "1em")
                .style("text-anchor", "end")
                .text("Diversity");
        } else {
            diversity_svg = d3.select("#" + id).append("div").
            text("No information found (probably due to low coverage/failed run)").
            attr("class", "alert  alert-warning");
        }
    }

    function render_histogram(json, dim, id, labels, patient_info, show_table) {
        patient_info = patient_info || "";

        var histogram_svg;
        clean_existing(id);

        var log_scale = true,
            no_label = false,
            x_range = null;

        if (json) {
            var counts,
                fst_plot = false,
                hist_labels = [];


            if (_.isArray(json)) {
                counts = [json[0]];
                x_range = json[1];
                //console.log (
                log_scale = false;
                no_label = true;
            } else {
                if (_.isObject(json)) {
                    if ("Histogram" in json) {
                        counts = [json["Histogram"]];
                        hist_labels = ['tn93'];
                    } else {
                        counts = []
                        for (k in json) {
                            if (_.isArray(json[k])) {
                                hist_labels.push(k);
                                counts.push(json[k].map(function(d) {
                                    return [d[0], d[1] < 1e12 ? d[1] : 0];
                                }));
                            }
                        }
                        log_scale = false;
                        fst_plot = true;
                    }
                } else {
                    return;
                }
            }


            //console.log (json


            var margin = {
                    top: 20,
                    right: 20,
                    bottom: 40,
                    left: 50
                },
                width = dim[0] - margin.left - margin.right,
                height = dim[1] - margin.top - margin.bottom;

            max_y = 0;

            if (!x_range) {
                max_x = 0;

                for (k = counts[0].length - 1; k > 0; k -= 1) {
                    if (counts[0][k][1]) {
                        max_x = counts[0][k][0];
                        break;
                    }
                }
                x_range = [0, max_x];
            }


            var step = counts[0][1][0] - counts[0][0][0];
            var x = d3.scale.linear()
                .domain([x_range[0], x_range[1]])
                .range([0, width]);

            var y;

            for (k = 0; k < counts.length; k++) {
                max_y = Math.max(max_y, d3.max(counts[k].map(function(d) {
                    return d[1];
                })));
                if (log_scale) {
                    counts[k] = counts[k].map(function(d) {
                        return [d[0], Math.max(d[1], 1)];
                    });
                }
            }


            if (log_scale) {
                scf = d3.scale.log;
            } else {
                scf = d3.scale.linear;
            }

            y = scf()
                .domain([log_scale ? 1 : 0, max_y])
                .range([height, 0]);

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");

            if (show_table) {
                counts_table = [];
                var last = 0;
                for (k = 0; k < show_table; k++) {
                    counts_table.push(["" + last + "-" + counts[0][k][0], counts[0][k][1], counts[0][k][1] + (k > 0 ? counts_table[k - 1][2] : 0), 0, 0]);
                    last = counts[0][k][0];
                }
                left_over = counts[0].reduce(function(p, c, i) {
                    if (i < show_table) return p;
                    return p + c[1];
                }, 0);
                counts_table.push(["&ge;" + counts[0][show_table][0], left_over, counts_table[show_table - 1][2] + left_over, 0]);
                counts_table.forEach(function(d) {
                    d[3] = d[1] / counts_table[show_table][2] * 100.;
                    d[4] = d[2] / counts_table[show_table][2] * 100.
                });

                var hc = d3.select("#" + id);
                if (patient_info) {
                    hc.append("h5").html(patient_info);
                }
                hc = hc.append("div").attr("class", "row");
                histogram_svg = hc.append("div").attr("class", "col-md-6").append("svg");
                hc = hc.append("div").attr("class", "col-md-6").append("table").attr("class", "table table-striped table-condensed");

                hc.append("thead").selectAll("tr").data([
                        ["Distance range", "Count", "Cumulative", "% of total", "Percentile"]
                    ]).enter().append("tr")
                    .selectAll("th").data(function(d) {
                        return d;
                    }).enter().append("th")
                    .text(function(d) {
                        return d;
                    });

                hc.append("tbody").selectAll("tr").data(counts_table).enter().append("tr").selectAll("td")
                    .data(function(d) {
                        return d;
                    }).enter().append("td")
                    .html(function(d, i) {
                        switch (i) {
                            case 0:
                            case 1:
                            case 2:
                                return d;
                                break;
                            case 3:
                            case 4:
                                return percentage_format(d);
                        }
                    });



            } else {
                if (!no_label) {
                    d3.select("#" + id).
                    append("h5").
                    html(fst_plot ? (patient_info + "Compartmentalization analysis: F<sub>ST</sub> = " + number_format(json["f_st"]) + ", p-value for compartmentalization = " + number_format(json["p"])) : (patient_info ? patient_info : ("Mean pairwise TN93 nucleotide distance = " + number_format(100 * json["Mean distance"]) + " %")));
                }
                histogram_svg = d3.select("#" + id).append("svg");
            }

            histogram_svg = histogram_svg.attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            for (k = 0; k < counts.length; k++) {

                var bar = histogram_svg.selectAll(".bar." + hist_labels[k])
                    .data(counts[k])
                    .enter().append("g")
                    .attr("class", "bar " + hist_labels[k])
                    .attr("transform", function(d, i) {
                        return "translate(" + x(d[0] - step) + "," + y(d[1]) + ")";
                    });

                bar.append("rect")
                    .attr("x", 1)
                    .attr("width", function(d, i) {
                        return x(step) - x(0);
                    })
                    .attr("height", function(d) {
                        return height - y(d[1]);
                    })
                    .append("title").text(function(d, i) {
                        return "" + d[1] + " values between " + (d[0] - step) + " and " + d[0];
                    });

                if (fst_plot) {
                    bar.selectAll("rect").style("fill-opacity", "0.2")
                        .style("fill", colors(k));
                }

                var line = d3.svg.line()
                    .x(function(d) {
                        return x(d[0] - step * 0.5);
                    })
                    .y(function(d) {
                        return y(d[1]);
                    })
                    .interpolate("monotone");

                histogram_svg.append("path")
                    .datum(counts[k])
                    .attr("class", "line-hist")
                    .style("stroke", fst_plot ? colors(k) : "#000")
                    .attr("d", line);

            }


            histogram_svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis);

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");

            var yAxis = d3.svg.axis()
                .scale(y)
                .orient("left");

            if (fst_plot) {
                yAxis.tickFormat(d3.format(".4s"));
            }

            histogram_svg.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis)
                .append("text")
                .attr("transform", "translate(" + width + ",0)")
                .attr("y", "24")
                .attr("dy", ".71em")
                .style("text-anchor", "end")
                .text(labels[0]);

            histogram_svg.append("g")
                .attr("class", "y axis")
                .call(yAxis)
                .append("text")
                //.attr("transform", "rotate(-90)")
                .attr("y", 6)
                .attr("dy", "-1.5em")
                .attr("dx", "1em")
                .style("text-anchor", "end")
                .text(labels[1]);

            if (fst_plot) {

                var legend_dim = {
                    x: width - 100,
                    y: 10,
                    spacer: 12,
                    margin: 5,
                    font: 8
                };

                var legend = histogram_svg.append("g")
                    .attr("class", "legend")
                    .attr("x", legend_dim.x)
                    .attr("y", legend_dim.y)
                    .attr("transform", "translate(" + legend_dim.x + "," + legend_dim.y + ")");

                legend.selectAll('g').data(hist_labels)
                    .enter()
                    .append('g')
                    .each(function(d, i) {
                        var g = d3.select(this);
                        g.append("rect")
                            .attr("x", legend_dim.spacer)
                            .attr("y", i * (legend_dim.spacer + legend_dim.margin))
                            .attr("width", legend_dim.spacer)
                            .attr("height", legend_dim.spacer)
                            .style("fill", colors(i));

                        g.append("text")
                            .attr("x", 2 * legend_dim.spacer + legend_dim.font / 4)
                            .attr("y", (i + 1) * (legend_dim.spacer + legend_dim.margin) - legend_dim.margin - (legend_dim.spacer - legend_dim.font) * 2 / 3)
                            .style("fill", colors(i))
                            .text(function(d) {
                                return d;
                            });

                    });


            }

        } else {
            histogram_svg = d3.select("#" + id).append("div").
            text("No information found (probably due to low coverage/failed run)").
            attr("class", "alert  alert-warning");
        }
    }

    function render_positional_table(id, has_dram) {
        id = "#" + id;
        d3.select(id).selectAll("tr").remove();

        var table_data = positional_data.length ? positional_data : [
            ["Low coverage", "", ""]
        ];

        //console.log (table_data);

        d3.select(id).selectAll("tr").data(table_data).enter().append("tr").
        selectAll("td").data(function(d) {
            return (d);
        }).enter().append("td").
        html(
            function(d, i, j) {
                if (i < 3) {
                    return d;
                }
                if (i == 3) {
                    html = '<ul class="list-inline">';
                    for (j in d) {
                        html += '<li style = "color:' + color_rates(d[j][2]) + '"><strong>' + d[j][0] + '</strong> ' + percentage_format(d[j][2] * 100) + '%</li>';
                    }
                    return html + "</ul>";
                }
                if (i == 4) {
                    return '<a href="#site_rate_display" data-toggle="modal" data-codon-id = "' + (parseInt(positional_data[j][1]) + (positional_data[j][0] == 'RT' && has_dram != 2 ? 100 : 0)) + '" data-placement = "bottom">' +
                        '<button type="button" class="btn btn-default btn-xs"><span class="glyphicon glyphicon-th"></span> View</button></a>';
                }
            }
        );

        $('[data-toggle="modal"]').on('click', function(event) {
            display_site_properties($(this), $(this).attr("data-codon-id"));
        });

    }

    function render_dram_table_wrapper(data_table) {
        var dram_data = data_table.filter(global_summary_filter);
        d3.select("#summary_matching").text(dram_data.length);
        prepare_csv_export([dram_column_names, render_dram_table("dram_table_body", dram_data), "dram.txt"]);

    }

    function render_dram_table(id, data_table, filter) {
        id = "#" + id;
        d3.select(id).selectAll("tr").remove();

        if (data_table.length) {
            var columns = filter ? 3 : data_table[0].length - 1;

            var rows = d3.select(id).selectAll("tr").data(data_table).enter().append("tr");
            rows.selectAll("td").data(function(d) {
                return (d);
            }).enter().append("td").
            html(
                function(d, i, k) {
                    if (i < columns) {
                        return d;
                    }
                    if (i == columns) {
                        html = '<ul class="list-inline">';
                        for (j in d) {
                            render = color_rates;
                            if (filter) {
                                position = data_table[k][0] + data_table[k][1];
                                if (d[j][0] in filter[position]) {
                                    render = color_dram;
                                }
                            }
                            html += '<li style = "color:' + render(d[j][2]) + '"><strong>' + d[j][0] + '</strong> ' + percentage_format(d[j][2] * 100) + '%</li>';

                        }
                        return html + "</ul>";
                    }
                    return "";
                }
            );
            return rows;
        }
        return null;

    }

    function extractPositionsByDRAM() {

        classes_to_check = {};
        class_map = {
            'pi_toggle': 'PI',
            'nrti_toggle': 'NRTI',
            'nnrti_toggle': 'NNRTI'
        };
        for (k in class_map) {
            if ($('#' + k).prop('checked')) {
                classes_to_check[class_map[k]] = true;
            }
        }

        threshold = parseFloat($('#dram_freq').val());
        filter_on = extractDRM(classes_to_check, null, dram_minimum_score);

        var positions = positional_data.filter(function(d) {
            var res_idx = d[0] + d[1];
            if (res_idx in filter_on) {
                for (var res in d[3]) {
                    var res_data = d[3][res];
                    if (res_data[0] in filter_on[res_idx] && res_data[2] >= threshold) {
                        return true;
                    }
                }
            }
            return false;
        });

        render_dram_table('positional_table_dram', positions, filter_on);
    }


    function get_site_info(a_site) {
        var sum = 0;
        var max = 0;
        var info = [];
        for (r in a_site) {
            sum += a_site[r];
            max = Math.max(max, a_site[r]);
            info.push([r, a_site[r], 0]);
        }
        info = info.map(function(d) {
            d[2] = d[1] / sum;
            return d;
        });
        info = info.sort(function(a, b) {
            return b[1] - a[1];
        });
        return [sum, max, info];
    }

    function site_annotator(k, has_dram) {
        if (has_dram != 2) {
            return [(k < 100 ? 'PR' : 'RT'), (k < 100 ? k : k - 99)]
        }
        return ['RT', k];
    }

    function plot_coverage_data(json, labels, id, dim, gene, has_dram) {

        positional_data = [];
        if (json) {

            coverage_data = []; // coordinate, total, majority
            for (k in json) {

                var site_info = get_site_info(json[k]),
                    sum = site_info[0],
                    max = site_info[1],
                    info = site_info[2];

                coverage_data.push([parseInt(k), sum, max]);


                if (has_dram >= 0) {
                    positional_data.push(site_annotator(k, has_dram).concat(sum, [info], ""));
                } else {
                    positional_data.push([gene, k, sum, info, ""]);
                }
            }

            render_positional_table("positional_table", has_dram);

            if (has_dram) {
                extractPositionsByDRAM(positional_data);
                d3.select("#dram_button").on("click", extractPositionsByDRAM);
            }
            //console.log (extractPositionsByDRAM (positional_data, extractDRM ({'NRTI':1}, null, 31), 0.01));



            var margin = {
                    top: 10,
                    right: 10,
                    bottom: 100,
                    left: 40
                },
                margin2 = {
                    top: dim[1] - 70,
                    right: 10,
                    bottom: 20,
                    left: 40
                },
                width = dim[0] - margin.left - margin.right,
                height = dim[1] - margin.top - margin.bottom,
                height2 = dim[1] - margin2.top - margin2.bottom;

            var x = d3.scale.linear()
                .range([0, width]);

            var x2 = d3.scale.linear()
                .range([0, width]);

            var y = d3.scale.linear()
                .range([height, 0]);

            var y2 = d3.scale.linear()
                .range([height2, 0]);

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient("bottom");

            var xAxis2 = d3.svg.axis()
                .scale(x2)
                .orient("bottom");

            var yAxis = d3.svg.axis()
                .scale(y)
                .orient("left");

            var svg = d3.select("#" + id)
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom);

            var brush = d3.svg.brush()
                .x(x2)
                .on("brush", brushed);

            svg.selectAll("path").remove();
            svg.selectAll("g").remove();
            svg.selectAll("defs").remove();


            x.domain(d3.extent(coverage_data, function(d, i) {
                return d[0];
            }));
            x2.domain(x.domain());
            y.domain([0, 1.2 * d3.max(coverage_data, function(d) {
                return d[1];
            })]);
            y2.domain(y.domain());



            var line = d3.svg.line()
                .x(function(d, i) {
                    return x(d[0]);
                })
                .y(function(d) {
                    return y(d[1]);
                });

            var area = d3.svg.area()
                .x(function(d, i) {
                    return x(d[0]);
                })
                .y1(function(d) {
                    return y(d[1]);
                })
                .y0(function(d) {
                    return height;
                })
                .interpolate("monotone");

            var area_maj = d3.svg.area()
                .x(function(d, i) {
                    return x(d[0]);
                })
                .y1(function(d) {
                    return y(d[2]);
                })
                .y0(function(d) {
                    return height;
                })
                .interpolate("monotone");

            var area2 = d3.svg.area()
                .x(function(d) {
                    return x2(d[0]);
                })
                .y0(height2)
                .y1(function(d) {
                    return y2(d[2]);
                })
                .interpolate("monotone");

            svg.append("defs").append("clipPath")
                .attr("id", "clip")
                .append("rect")
                .attr("width", width)
                .attr("height", height);


            function brushed() {
                x.domain(brush.empty() ? x2.domain() : brush.extent());
                main_path.attr("d", area);
                secondary_path.attr("d", area_maj);
                main_line.attr("d", line);
                main_axis.call(xAxis);
            }

            var focus = svg.append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            var context = svg.append("g")
                .attr("transform", "translate(" + margin2.left + "," + margin2.top + ")");


            var main_path = focus.append("path")
                .datum(coverage_data)
                .attr("clip-path", "url(#clip)")
                .attr("class", "area")
                .style("fill", colors(1))
                .attr("d", area);

            var secondary_path = focus.append("path")
                .datum(coverage_data)
                .attr("clip-path", "url(#clip)")
                .attr("class", "area_maj")
                .style("fill", colors(0))
                .attr("d", area_maj);

            var main_line = focus.append("path")
                .datum(coverage_data)
                .attr("clip-path", "url(#clip)")
                .attr("class", "line")
                .attr("d", line);


            context.append("path")
                .datum(coverage_data)
                .attr("class", "area")
                .style("fill", colors(0))
                .attr("d", area2);

            context.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height2 + ")")
                .call(xAxis2);


            context.append("g")
                .attr("class", "x brush")
                .call(brush)
                .selectAll("rect")
                .attr("y", -6)
                .attr("height", height2 + 7);

            var legend_dim = {
                x: 45,
                y: 10,
                spacer: 12,
                margin: 5,
                font: 8
            };

            var legend = svg.append("g")
                .attr("class", "legend")
                .attr("x", legend_dim.x)
                .attr("y", legend_dim.y)
                .attr("transform", "translate(" + legend_dim.x + "," + legend_dim.y + ")");

            legend.selectAll('g').data(labels.slice(2))
                .enter()
                .append('g')
                .each(function(d, i) {
                    var g = d3.select(this);
                    g.append("rect")
                        .attr("x", legend_dim.spacer)
                        .attr("y", i * (legend_dim.spacer + legend_dim.margin))
                        .attr("width", legend_dim.spacer)
                        .attr("height", legend_dim.spacer)
                        .style("fill", colors(i));

                    g.append("text")
                        .attr("x", 2 * legend_dim.spacer + legend_dim.font / 4)
                        .attr("y", (i + 1) * (legend_dim.spacer + legend_dim.margin) - legend_dim.margin - (legend_dim.spacer - legend_dim.font) * 2 / 3)
                        .style("fill", colors(i))
                        .text(function(d) {
                            return d;
                        });

                });


            var main_axis = focus.append("g")
                .attr("class", "x axis")
                .attr("transform", "translate(0," + height + ")")
                .call(xAxis);

            main_axis.append("text")
                .attr("transform", "translate(" + width + ",0)")
                .attr("y", "24")
                .attr("dy", ".71em")
                .style("text-anchor", "end")
                .text(labels[0]);

            focus.append("g")
                .attr("class", "y axis")
                .call(yAxis)
                .append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 6)
                .attr("dy", ".71em")
                .style("text-anchor", "end")
                .text(labels[1]);
        } else {
            render_positional_table("positional_table", has_dram);
        }
    }

    function plot_div_data(data, keys, labels, id, dim) {


        var margin = {
                top: 20,
                right: 150,
                bottom: 40,
                left: 100
            },
            width = dim[0] - margin.left - margin.right,
            height = dim[1] - margin.top - margin.bottom;


        var svg = d3.select("#" + id)
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom);

        svg.selectAll("path").remove();
        svg.selectAll("g").remove();

        divergence_data = [];
        raw_dates = [];

        var replicate = null;
        var compartment = null;
        var valid_dates = true;

        for (k in data) {
            switch (k) {


                case "label":
                    break;

                case "ID":
                    {
                        pid = data[k];
                        break;
                    }
                case "gene":
                    {
                        gene = data[k];
                        break;
                    }
                case "compartment":
                    {
                        compartment = data[k];
                        break;
                    }
                case "replicate":
                    {
                        replicate = data[k];
                        break;
                    }
                default:
                    {

                        raw_dates.push(k);
                        info = [parse_date.parse(k)];

                        if (!info[0]) {
                            valid_dates = false;
                            for (z in divergence_data) {
                                divergence_data[z][0] = raw_dates[z];
                            }
                            info = [k];
                        }

                        for (k2 in keys) {
                            val = (parseFloat((data[k])[keys[k2]]));
                            info.push(isNaN(val) ? 0 : val);
                        }
                        divergence_data.push(info);
                    }
            }
        }

        var x;
        if (valid_dates) {
            x = d3.time.scale()
                .range([0, width]);
        } else {
            x = d3.scale.ordinal()
                .range([0, width]);
        }


        var y = d3.scale.linear()
            .range([height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient("left");

        divergence_data.sort(function(a, b) {
            if (a[0] < b[0]) return -1;
            if (a[0] > b[0]) return 1;
            return 0;
        });



        x.domain(d3.extent(divergence_data, function(d, i) {
            return d[0];
        }));
        y.domain([0, d3.max(divergence_data, function(d) {
            var values = [d[1] + d[2] / 2];
            for (k = 3; k < d.length; k++) {
                values.push(d[k]);
            }
            return d3.max(values);
        })]);


        var plot_svg = svg.append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


        var max_x = 0;
        var circle_gen = d3.svg.symbol().size(function(d) {
            var r = y(0) - y(d[2]);
            max_x = Math.max(max_x, x(d[0]) + r);
            return r * r;
        });

        var colors = d3.scale.category10();

        for (k in divergence_data) {
            plot_svg.append("path")
                .datum(divergence_data[k])
                .attr("transform", function(d) {
                    return "translate(" + x(d[0]) + "," + y(d[1]) + ")";
                })
                .attr("d", circle_gen)
                .attr("class", "_evo_symbol")
                .append("title")
                .text(function(d) {
                    return "Diversity = " + d[2].toFixed(5);
                });
        }

        var overflow = max_x - width;
        if (overflow > 0) {
            svg.attr("width", width + margin.left + margin.right + overflow);
        }

        var line = d3.svg.line()
            .x(function(d, i) {
                return x(d[0]);
            })
            .y(function(d) {
                return y(d[1]);
            });


        plot_svg.append("path")
            .datum(divergence_data)
            .attr("class", "_evo_line")
            .style('stroke', colors(0))
            .attr("d", line);

        for (k = 2; k < keys.length; k++) {
            var extra_line = d3.svg.line()
                .x(function(d, i) {
                    return x(d[0]);
                })
                .y(function(d) {
                    return y(d[k + 1]);
                });
            plot_svg.append("path")
                .datum(divergence_data)
                .attr("class", "_evo_line")
                .style('stroke', colors(k))
                .attr("d", extra_line);
        }

        var legend_dim = {
            x: dim[0] - 120,
            y: 20,
            spacer: 25,
            margin: 5,
            font: 12
        };

        var legend = svg.append("g")
            .attr("class", "_evo_legend")
            .attr("x", legend_dim.x)
            .attr("y", legend_dim.y)
            .attr("transform", "translate(" + legend_dim.x + "," + legend_dim.y + ")");

        legend.selectAll('g').data(labels)
            .enter()
            .append('g')
            .each(function(d, i0) {
                var g = d3.select(this);
                g.append("rect")
                    .attr("x", legend_dim.spacer)
                    .attr("y", i0 * (legend_dim.spacer + legend_dim.margin))
                    .attr("width", legend_dim.spacer)
                    .attr("height", legend_dim.spacer)
                    .style("fill", function() {
                        if (i0 == 1) return "#DDDDDD";
                        return colors(i0);
                    });

                g.append("text")
                    .attr("x", 2 * legend_dim.spacer + legend_dim.font / 4)
                    .attr("y", (i0 + 1) * (legend_dim.spacer + legend_dim.margin) - legend_dim.margin - (legend_dim.spacer - legend_dim.font) * 2 / 3)
                    .style("fill", function() {
                        if (i0 == 1) return "#DDDDDD";
                        return colors(i0);
                    })
                    .text(function(d) {
                        return d;
                    });

            });


        svg.append("g")
            .attr("class", "_evo_legend")
            .append("text")
            .attr("x", margin.left)
            .attr("y", 10)
            .text(pid + " " + gene + (compartment ? " from " + compartment : "") + (replicate && use_replicate_for_intrahost ? ", replicate " + replicate : ""));

        plot_svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + height + ")")
            .call(xAxis)
            .append("text")
            .attr("transform", "translate(" + width + ",0)")
            .attr("y", "24")
            .attr("dy", ".71em")
            .style("text-anchor", "end")
            .text("Sample date");

        plot_svg.append("g")
            .attr("class", "y axis")
            .call(yAxis)
            .append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 6)
            .attr("dy", ".71em")
            .style("text-anchor", "end")
            .text("Divergence from MRCA");

    }

    function _describe_a_vector(vector) {
        vector.sort(function(a, b) {
            return a - b
        });
        var total_sum = _.reduce(vector, function(sum, value) {
                return sum + value;
            }, 0),
            mean = total_sum / vector.length,
            midpoint = Math.floor(vector.length / 2),
            median = vector.length % 2 == 0 ? (vector[midpoint - 1] + vector[midpoint]) / 2 : vector[midpoint],
            q1 = vector[Math.floor(vector.length / 4)],
            q3 = vector[Math.floor(vector.length * 3 / 4)],
            min = vector[0],
            max = vector[vector.length - 1],
            step = 2 * (q3 - q1) / Math.pow(vector.length, 1 / 3),
            discretized = [];

        vector.forEach(function(d) {
            var index = Math.floor((d - min) / step);
            if (discretized[index]) {
                discretized[index][1] += 1;
            } else {
                discretized[index] = [min + step * (index + 1), 1];
            }
        });


        for (k = 0; k < discretized.length; k += 1) {
            if (!discretized[k]) {
                discretized[k] = [min + step * (k + 1), 0];
            }
        }


        return [mean, median, q1, q3, mean, total_sum, discretized, min, min + step * (discretized.length)];
    }

    function _set_rows(container, rows, join_key) {
        var rows = d3.select("#" + container).selectAll("tr").data(rows, join_key);
        rows.enter().append('tr');
        rows.exit().remove();
        return rows;
    }
