var SIMPATH_DATA;


$(document).ready(function(){
    var jsonPath = document.location.href
    $.getJSON('../simpath.json', function(data) {
        SIMPATH_DATA = data;
        for (gs in data) {
            $('.SimPathSelect#GetSet').append($('<option></option>').val(gs).html(gs));
        }
        for (solver in data['get']) {
            $('.SimPathSelect#Solver').append($('<option></option>').val(solver).html(solver));
        }
        for (loc in data['get']['Wmdirect']) {
            $('.SimPathSelect#Location').append($('<option></option>').val(loc).html(loc));
        }
    });

    var ids = {"GetSet": "Solver", "Solver": "Location", "Location":"Item", "Item":"Property", "Property":"None"};
    for (selectid in ids) {
        $('.SimPathSelect#' + selectid).change(function() {
            var key = $(this).val();
            var sid = $(this).attr("id");

            var tmp_dict = SIMPATH_DATA;
            for (field in ids) { 
                if (field === sid) {
                    break;
                }
                tmp_dict = tmp_dict[$('.SimPathSelect#' + field).val()];
            }

            var currId = sid;
            var old_select = {};
            while (currId in ids) {
                currId = ids[currId];
                old_select[currId] = $('.SimPathSelect#' + currId).val();
                $('.SimPathSelect#' + currId).empty();
            }
            $('.DisplayBox#SimPath').text("");
            $('.DocDisplay#SimPathDoc').text("");

            currId = ids[sid];
            while (currId in ids) {
                if (typeof tmp_dict[key] === 'string') {
                    currId = "None";
                    break;
                }
                for (val in tmp_dict[key]) {
                    $('.SimPathSelect#' + currId).append($('<option></option>').val(val).html(val));
                }
                if (old_select[currId] in tmp_dict[key]) {
                    $('.SimPathSelect#' + currId).val(old_select[currId])
                    tmp_dict = tmp_dict[key];
                    key = old_select[currId];
                    currId = ids[currId];
                } else {
                    break;
                }
            }
            if (currId === "None") {
                var str = tmp_dict[key];
                var line = str.substr(0, str.indexOf("@"));
                var doc = str.substr(str.indexOf("@")+1);
                $('.DisplayBox#SimPath').html(line);
                $('.DocDisplay#SimPathDoc').html(doc);
            }
        });
    }
});