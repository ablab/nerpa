function lookAfterSelectElem() {
    var x, i, j, selElmnt, a, b, c;
    /*look for any elements with the class "custom-select":*/
    x = document.getElementsByClassName("custom-select");
    for (i = 0; i < x.length; i++) {
        selsel = x[i].getElementsByClassName("select-selected");
        if (selsel.length != 0) {
            continue;
        }
        selElmnt = x[i].getElementsByTagName("select")[0];
        /*for each element, create a new DIV that will act as the selected item:*/
        a = document.createElement("DIV");
        a.setAttribute("class", "select-selected");
        a.innerHTML = selElmnt.options[selElmnt.selectedIndex].innerHTML;
        x[i].appendChild(a);
        /*for each element, create a new DIV that will contain the option list:*/
        b = document.createElement("DIV");
        b.setAttribute("class", "select-items select-hide");
        for (j = 0; j < selElmnt.length; j++) {
            /*for each option in the original select element,
             create a new DIV that will act as an option item:*/
            c = document.createElement("DIV");
            c.innerHTML = selElmnt.options[j].innerHTML;
            c.addEventListener("click", function (e) {
                /*when an item is clicked, update the original select box,
                 and the selected item:*/
                var y, i, k, s, h;
                s = this.parentNode.parentNode.getElementsByTagName("select")[0];
                h = this.parentNode.previousSibling;
                for (i = 0; i < s.length; i++) {
                    if (s.options[i].innerHTML == this.innerHTML) {
                        s.selectedIndex = i;
                        h.innerHTML = this.innerHTML;
                        y = this.parentNode.getElementsByClassName("same-as-selected");
                        for (k = 0; k < y.length; k++) {
                            y[k].removeAttribute("class");
                        }
                        this.setAttribute("class", "same-as-selected");
                        break;
                    }
                }
                updateForms();
                h.click();
            });
            b.appendChild(c);
        }
        x[i].appendChild(b);
        a.addEventListener("click", function (e) {
            /*when the select box is clicked, close any other select boxes,
             and open/close the current select box:*/
            e.stopPropagation();
            closeAllSelect(this);
            this.nextSibling.classList.toggle("select-hide");
            this.classList.toggle("select-arrow-active");
        });
    }
}

function uploadFile(idname, suf) {
    $("#" + idname).change(function(){
        var file_name;
        var inp = $("#" + idname);
        var btn = $("#uploadButton" + suf);
        var lbl = $("#UploadInput" + suf);
        if( file_api && inp[ 0 ].files[ 0 ] )
            file_name = inp[ 0 ].files[ 0 ].name;
        else
            file_name = inp.val().replace( "C:\\fakepath\\", '' );

        if( ! file_name.length )
            return;

        lbl.text( file_name );
    });
}


function updateForms() {
    var elem = document.getElementById("id_search_type");
    if (elem.value === 'genome') {
        document.getElementById("upload_block").innerHTML = '<p><b>Genome sequence:</b></p>' +
                '<div class="file_upload" id="file_upload">' +
                '    <button id="uploadButton">' +
                '        Upload' +
                '    </button>' +
                '    <div class="UploadInput" id="UploadInput">' +
                '    No file chosen' +
                '   </div>' +
                '   <input type="file" name="inputFileGenome" required="" id="id_inputFileGenome">' +
                '</div>';

        document.getElementById("choose_db").innerHTML = '<p><b>NRPs Data Base:</b></p>'+
            '<div class="custom-select" style="width:600px;">' +
            '<select id="id_nrp_db" name="nrp_db">' +
            '    <option value="streptome">Streptome</option>' +
            '</select>' +
            '</div>';

        uploadFile("id_inputFileGenome", "");
    } else if (elem.value === 'nrp') {
        document.getElementById("upload_block").innerHTML = '<p><b>NRP structure:</b></p>' +
                '<div class="file_upload" id="file_upload">' +
                '    <button id="uploadButton">' +
                '        Upload' +
                '    </button>' +
                '    <div class="UploadInput" id="UploadInput">' +
                '    No file chosen' +
                '   </div>' +
                '   <input type="file" name="inputFileNRP" required="" id="id_inputFileNRP">' +
                '</div>';

        document.getElementById("choose_db").innerHTML = '<p><b>Genome database:</b></p>'+
            '<div class="custom-select" style="width:600px;">' +
            '<select id="id_genome_db" name="genome_db">' +
            '    <option value="bc">Bacteria complete</option>' +
            '</select>' +
            '</div>';


        uploadFile("id_inputFileNRP", "")
    } else if (elem.value === 'one') {
        document.getElementById("upload_block").innerHTML = '<p><b>NRP structure:</b></p>' +
                '<div class="file_upload" id="file_upload">' +
                '    <button id="uploadButton">' +
                '        Upload' +
                '    </button>' +
                '    <div class="UploadInput" id="UploadInput">' +
                '    No file chosen' +
                '   </div>' +
                '   <input type="file" name="inputFileNRP" required="" id="id_inputFileNRP">' +
                '</div>';


        document.getElementById("choose_db").innerHTML = '<p><b>Genome sequence:</b></p>' +
                '<div class="file_upload" id="file_upload">' +
                '    <button id="uploadButton2">' +
                '        Upload' +
                '    </button>' +
                '    <div class="UploadInput" id="UploadInput2">' +
                '    No file chosen' +
                '   </div>' +
                '   <input type="file" name="inputFileGenome" required="" id="id_inputFileGenome">' +
                '</div>';

        uploadFile("id_inputFileNRP", "");
        uploadFile("id_inputFileGenome", "2");
    }

    lookAfterSelectElem();
}

lookAfterSelectElem();

function closeAllSelect(elmnt) {
  /*a function that will close all select boxes in the document,
  except the current select box:*/
  var x, y, i, arrNo = [];
  x = document.getElementsByClassName("select-items");
  y = document.getElementsByClassName("select-selected");
  for (i = 0; i < y.length; i++) {
    if (elmnt == y[i]) {
      arrNo.push(i)
    } else {
      y[i].classList.remove("select-arrow-active");
    }
  }
  for (i = 0; i < x.length; i++) {
    if (arrNo.indexOf(i)) {
      x[i].classList.add("select-hide");
    }
  }
}
/*if the user clicks anywhere outside the select box,
then close all select boxes:*/
document.addEventListener("click", closeAllSelect);
