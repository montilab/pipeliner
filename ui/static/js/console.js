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

  $(document).on("click", "#nextflow-start", function() {
    socket.emit('nextflow_start', {resuming: false});
    return false;
  });

  $(document).on("click", "#nextflow-stop", function() {
    socket.emit('nextflow_stop', {});
    return false;
  });

  $(document).on("click", "#nextflow-resume", function() {
    socket.emit('nextflow_start', {resuming: true});
    return false;
  });
});