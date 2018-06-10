var file_api = ( window.File && window.FileReader && window.FileList && window.Blob ) ? true : false;

$(":file:first").change(function(){
    var file_name;
    var inp = $(":file:first");
    var btn = $("#uploadButton");
    var lbl = $("#UploadInput");
    if( file_api && inp[ 0 ].files[ 0 ] )
        file_name = inp[ 0 ].files[ 0 ].name;
    else
        file_name = inp.val().replace( "C:\\fakepath\\", '' );

    if( ! file_name.length )
        return;

    lbl.text( file_name );
});

/*$("#mainsearch").click(function () {
    alert("run search");
    $.get("/search/", function (data) {
    }, "html");
});*/