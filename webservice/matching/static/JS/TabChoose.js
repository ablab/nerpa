/**
 * Created by olga on 29.08.19.
 */

function changeActive(TabName) {
    var div = document.getElementById('header-right');
    var children = div.childNodes;
    for (var i = 0; i < children.length; ++i) {
        var child = children[i];
        if (child.className == "active") {
            child.classList.remove("active");
        }
        if (child.textContent == TabName) {
            child.classList.add("active");
        }
    }
}