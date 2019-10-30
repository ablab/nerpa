function download_csv() {
    var data_query = generate_data_query();
    data_query["DOWNLOAD"] = "CSV";
    document.getElementById("download-button").blur();
    
    $.ajax({
        type: "GET",
        url: this.href,
        data: data_query,
        success: function (data) {
            console.log(data);

            var hiddenElement = document.createElement('a');
            hiddenElement.href = 'data:text/csv;charset=utf-8,' + data;
            hiddenElement.target = '_blank';
            hiddenElement.download = 'nerpa_report.csv';
            hiddenElement.click();
        }
    });
}

function download_one_result() {
    var data_query = {};
    data_query["DOWNLOAD"] = "ONE_RESULT";

    $.ajax({
        type: "GET",
        url: this.href,
        data: data_query,
        success: function (data) {
            console.log(data);

            var hiddenElement = document.createElement('a');
            hiddenElement.href = 'data:text;charset=utf-8,' + data;
            hiddenElement.target = '_blank';
            hiddenElement.download = 'nerpa_result';
            hiddenElement.click();
        }
    });
}