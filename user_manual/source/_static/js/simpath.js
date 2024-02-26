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
            $('.ExamplesDisplay#SimPathExamples').text("");

            currId = ids[sid];
            while (currId in ids) {
                if ('@doc' in tmp_dict[key]) {
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
                var lines = tmp_dict[key]['@doc'];
                for (var i = 0 ; i < lines.length ; i++) {
                    $('<div class="DisplayBox line-block" style="font-size:25px;line-height:30px;width:610px;"><code class="py py-class">' + lines[i]['@code'] + '</code></div>').appendTo($('.ExamplesDisplay#SimPathExamples'));
                    $('<div class="DocDisplay line-block" style="width:610px;color:#666666;">' + lines[i]['@descr'] + '</div>').appendTo($('.ExamplesDisplay#SimPathExamples'));
                }
            }
        });
    }
});
