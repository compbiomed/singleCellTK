var busy = 0;
$(function(){ 
    $(document).on('shiny:connected', function(event) {
    $.busyLoadSetup({ animation: "slide", background: "rgba(245, 245, 245, 0.86)" });
    if(busy == 0){
      $.busyLoadFull("show",{
        text: "Loading...",
        textColor: "gray",
        color: "gray",
      });
      busy = 1;
    }
  }); 

    $(document).on('shiny:idle', function(event) {
    $.busyLoadFull("hide");
  });

});
