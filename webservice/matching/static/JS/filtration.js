function Groups(name, type, title) {
    this.name = name;
    this.type = type;
    this.title = title;
}

var groups = [new Groups("All results", "none", "All results")];

var title_prefix = {"genome_id": "Genome id: ", "structure_id": "Structure id: ", "product": "Product: ", "BGC": "BGC: "};
var select_option_text = {"genome_id": "Genome id", "structure_id": "Structure id", "BGC": "BGC", "product": "Product"};
function updtae_group_by_list() {
    var selectobject = document.getElementById("select_group_by");

    function type_in_groups(type_value) {
        if (type_value == "none") {
            return false;
        }
        for (var i = 0; i < groups.length; ++i) {
            if (groups[i].type == type_value) {
                return true;
            }
        }
        return false;
    }

    for (var i = 0; i < selectobject.length; ++i) {
        selectobject.options[i].disabled = false;
    }


    for (i = 0; i < selectobject.length; ++i) {
        if (type_in_groups(selectobject.options[i].value)) {
            selectobject.options[i].disabled = true;
        }
    }
}

function group_by_update() {
    var selector = document.getElementById('select_group_by');
    var choose_option = selector.options[selector.selectedIndex].value;
    $.ajax({
        type: "GET",
        url: this.href,
        data: {"request_type": "group_by", "value": choose_option},
        success: function (data) {
            document.getElementById('result_container').innerHTML = data;
        }
    });
}

function update_navigation() {
    var groupby_ul = document.getElementById('groupby_list');
    groupby_ul.innerHTML = "";
    for (var elem in groups) {
        groupby_ul.innerHTML += "<li class=\"inline_item\" onclick='change_group(" + elem + ")'> <a><span>" + groups[elem].title + "</span></a></li>\n"
    }
}


function generate_data_query() {
    var data_query = {};
    for (var elem in groups) {
        if (groups[elem].type != "none") {
            data_query[groups[elem].type] = groups[elem].name;
        }
    }

    var min_len = document.getElementById("min_len").value;
    var min_score = document.getElementById("min_score").value;

    data_query["min_len"] = min_len;
    data_query["min_score"] = min_score;

    return data_query;
}

function generate_query() {
    var data_query = generate_data_query();

    $.ajax({
        type: "GET",
        url: this.href,
        data: data_query,
        success: function (data) {
            document.getElementById('result_container').innerHTML = data;
        }
    });
}

function change_group(eid) {
    groups = groups.slice(0, eid + 1);
    update_navigation();
    generate_query();
    updtae_group_by_list();
}


function choose_group(type, value) {
    var title = "";
    if (type == "BGC") {
        var res = value.split("__ctg");
        title = "BGC: " + res[0] + " cluster #" + res[1];
    } else {
        title = title_prefix[type] + value;
    }

    groups.push(new Groups(value, type, title));
    document.getElementById('select_group_by').value = "none";

    update_navigation();
    generate_query();
    updtae_group_by_list();
}