// -------------------------------------------------------------------------- //
// Global Variables                                                           //
// -------------------------------------------------------------------------- //
var FILES = []

var STEPS = {"input": false, 
             "settings": true, 
             "output": false}

const filetypes = ['annotation-files', 
                   'reference-files', 
                   'reads-files', 
                   'alignment-files']

const option_1 = [filetypes[0], filetypes[1], filetypes[2]]
var   option_1_valid = option_1.reduce(function(obj, x) {obj[x] = false; return obj;}, {});
const option_2 = [filetypes[3]]
var   option_2_valid = option_2.reduce(function(obj, x) {obj[x] = false; return obj;}, {});

// -------------------------------------------------------------------------- //
// Helper Functions                                                           //
// -------------------------------------------------------------------------- //
function showMessage(id, message, cleanup=false) {
  document.getElementById(id).innerHTML = message;
  if (!cleanup) {
    setTimeout(function(){showMessage(id, " ", true)}, 3000)
  }
}
function editElement(id, attr, edit) {
  $('#'+id).attr(attr, edit)
}
// -------------------------------------------------------------------------- //
// Input                                                                      //
// -------------------------------------------------------------------------- //
$(document).on("click", ".input.delete.icon", function (event) {
  var label = $(event.target).closest("div")
  label.remove()
  for (var filetype in FILES) {
    for (var filename in FILES[filetype]) {
      if (FILES[filetype][filename] == label.prop("id")) {
        FILES[filetype].splice(filename, 1)
      }
    }
  }
  refreshInput() 
})

//------------------------------------------------------------------------------
function addSample(s) {
  return "<tr><td>"+s[0]+"</td><td>"+s[1]+"</td><td>"+s[2]+"</td></tr>"
}
function parseReads(pathtocsv) {
  $.post({
    type: "POST",
    url: "/parsereads",
    data: {"pathtocsv": pathtocsv},
    success: function(response){
      var success = JSON.parse(response)['success']
      var reads = JSON.parse(response)['reads']
      document.getElementById('reads-found').innerHTML = ""
      if (success) {
        document.getElementById('reads-accordian').hidden = false
        for (var sample in reads) {
          document.getElementById('reads-found').innerHTML += addSample(reads[sample])
        }
      }
    }
  })
}

// -----------------------------------------------------------------------------
function validateOption(array, dictionary, filetype, bool) {
  if (_.contains(array, filetype)) {
    dictionary[filetype] = bool
  }
}
function validateFiletype(filetype) {
  if (FILES[filetype].length == 1) {
    validateOption(option_1, option_1_valid, filetype, true)
    validateOption(option_2, option_2_valid, filetype, true)
    editElement(filetype+'-status', 'class', 'large green check icon')
  } else {
    validateOption(option_1, option_1_valid, filetype, false)
    validateOption(option_2, option_2_valid, filetype, false)
    editElement(filetype+'-status', 'class', 'large red remove icon')
  }
}
function validateOptions(callback=stepsComplete) {
  if (_.every(option_1_valid, _.first()) || _.every(option_2_valid, _.first())) {
    STEPS['input'] = true
    editElement('input-step', 'class', 'completed step') 
  } else {
    STEPS['input'] = false
    editElement('input-step', 'class', 'step')
  }
  callback()
}
function refreshInput(callback=validateOptions) {
  parseReads($("#input-path").val()+'/'+FILES['reads-files'][0])
  for (var filetype in FILES) {
    validateFiletype(filetype)
  }
  callback()
}
// -----------------------------------------------------------------------------

function clearFiles() {
  for (filetype in filetypes) {
    document.getElementById(filetypes[filetype]).innerHTML = ""
  }
}
function addFile(filename) {
  return "<div id='"+filename+"' style='display:inline;'>   \
            <a class='ui blue large label'>                 \
              "+filename+"<i class='input delete icon'></i> \
            </a>                                            \
         </div>"
}
function input() {
  clearFiles()
  $.post({
    type: "POST",
    url: "/input",
    data: {"input": $("#input-path").val()},
    success: function(response){
      FILES = JSON.parse(response)['files']
      for (var filetype in FILES) {
        for (var filename in FILES[filetype]) {
          document.getElementById(filetype).innerHTML += addFile(FILES[filetype][filename])
        }
      }
      refreshInput()
      var success = JSON.parse(response)['success']
      var message = JSON.parse(response)['message']
      if (!success) {
        showMessage("input-message", message)
      }
    }
  })
}
// -------------------------------------------------------------------------- //
// Settings                                                                   //
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Output                                                                     //
// -------------------------------------------------------------------------- //
function output() {
  $.post({
    type: "POST",
    url: "/output",
    data: {"output": $("#output-path").val()},
    success: function(response){
      var success = JSON.parse(response)['success']
      if (success) {
        STEPS['output'] = true
        editElement('output-step', 'class', 'completed step')
        editElement('output-path-status', 'class', 'active large green check icon')
      } else {
        STEPS['output'] = false
        editElement('output-step', 'class', 'step')
        editElement('output-path-status', 'class', 'active large red remove icon')
      }
      stepsComplete()
    }
  })
}

// -------------------------------------------------------------------------- //
// Compelete                                                                  //
// -------------------------------------------------------------------------- //
function stepsComplete() {
  console.log(STEPS)
  if (_.every(STEPS, _.first())) {
    editElement('config', 'class', 'ui positive button') 
    editElement('submit', 'class', 'ui positive button')
    return true
  } else {
    editElement('config', 'class', 'ui black button') 
    editElement('submit', 'class', 'ui black button') 
    return false    
  }
}
$(document).on("click", "#config", function() {
  if (stepsComplete()) {
    console.log(FILES)
  }
})
$(document).on("click", "#submit", function() {
  if (stepsComplete()) {
    console.log(FILES)
  }
})

// -------------------------------------------------------------------------- //
// Actions                                                                    //
// -------------------------------------------------------------------------- //
$(document).on("blur", "#output-path", output)
$(document).on("click", "#input-button", input)
$(document).keypress(function(e) {
  if (e.which === 13) {
    var focus = document.activeElement.id
    if (focus == "input-path") {input()} 
    else if (focus == "output-path") {output()}
  }
})