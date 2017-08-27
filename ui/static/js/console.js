function print(text) {
  var screen = $('#console')
  var history = $('#console').val()
  screen.val(history+'$ '+text);
  var textarea = document.getElementById('console');
  textarea.scrollTop = textarea.scrollHeight;
};

$(document).ready(function() {
  var socket = io.connect('http://' + document.domain + ':' + location.port);

  socket.on('connect', function() {
    console.log("Web socket connected...")

    socket.on('response', function(msg) {
      print(msg['data'])
    });
  });

  function config(nextflow=false, resuming=false, local=true) {
    if (!local) {
      nfdir = $("#nextflow-path").val()
    } else {
      nfdir = ""
    }
    $.post({
      type: "POST",
      url: "/config",
      data: {"indir"    : $("#input-path").val(),
             "outdir"   : $("#output-path").val(),
             "files"    : JSON.stringify(FILES),
             "settings" : JSON.stringify(SETTINGS),
             "nfdir"    : nfdir},
      success: function(response){
        var success = JSON.parse(response)['success']
        if (success){
          if (nextflow){
            if (NEXTFLOW['nextflow-path'] && NEXTFLOW['pipeline-name']){
              socket.emit('nextflow_start', {resuming : false, 
                                             nfdir    : nfdir, 
                                             pipeline : pipeline = $("#pipeline-name").val(), 
                                             env      : env = $("#environment-name").val()})
            } else {
              return false
            }
          }
        }
      }
    }) 
  }
  $(document).on("click", "#export", function() {
    console.log(stepsComplete())
    if (stepsComplete()) {
      config()
    }
  })  

  $(document).on("click", "#nextflow-start", function() {
    config(nextflow=true, resuming=false, local=false)
    return false;
  });

  $(document).on("click", "#nextflow-stop", function() {
    socket.emit('nextflow_stop', {});
    return false;
  });

  $(document).on("click", "#nextflow-resume", function() {
    config(nextflow=true, resuming=true, local=false)
    return false;
  });
});