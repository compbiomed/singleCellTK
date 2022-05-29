document.getElementById("ia-button").click();
function openTab(evt, tabName, tabset = null) {
  // Previous version should work normally
  if (tabset == null) {
      var i, tabcontent, tablinks;
      tabcontent = document.getElementsByClassName("tabcontent");
      for (i = 0; i < tabcontent.length; i++) {
          tabcontent[i].style.display = "none";
      }
      tablinks = document.getElementsByClassName("tablinks");
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
      }
      document.getElementById(tabName).style.display = "block";
      evt.currentTarget.className += " active";
  } else {
    // New condition for a la carte workflow tutorial, with tabset id specified
    // Here I'm using `tabName` as class name instead of id!!
    
    var tabcontents = $('#'+tabset).find('.tabcontent');
    for (i=0; i<tabcontents.length; i++) {
        if (tabcontents[i].classList.contains(tabName)) {
            tabcontents[i].style.display = "block";
        } else {
            tabcontents[i].style.display = "none";
        }
    }
    var tabLinks = $('#'+tabset).find('.tablinks');
    for (i=0; i<tabLinks.length; i++) {
        if (tabLinks[i].classList.contains(tabName)) {
            tabLinks[i].classList.add("active");
        } else {
            tabLinks[i].classList.remove("active");
        }
    }
  }
}