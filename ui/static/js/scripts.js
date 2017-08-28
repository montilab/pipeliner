// -------------------------------------------------------------------------- //
// Global Variables                                                           //
// -------------------------------------------------------------------------- //
var FILES = {}
var SETTINGS = {'aligner'        : 'star',
                'star_index'     : false,
                'paired'         : true,
                'save_reference' : true}

var NEXTFLOW = {"nextflow-path"    : false,
                "pipeline-name"    : false,
                "environment-name" : false}

var STEPS = {"input"    : false, 
             "settings" : true, 
             "output"   : false}

const filetypes = ['annotation-files', 
                   'reference-files', 
                   'csv-files', 
                   'alignment-files']

const option_1 = [filetypes[0], filetypes[1], filetypes[2]]
var   option_1_valid = option_1.reduce(function(obj, x) {obj[x] = false; return obj;}, {});
const option_2 = [filetypes[3]]
var   option_2_valid = option_2.reduce(function(obj, x) {obj[x] = false; return obj;}, {});

const SUCCESS = "active large green check icon"
const FAILURE = "active large red remove icon"

// -------------------------------------------------------------------------- //
// Helper Functions                                                           //
// -------------------------------------------------------------------------- //
function message(msg, cleanup=false) {
  document.getElementById('message').innerHTML = msg;
  if (!cleanup) {
    setTimeout(function(){message(" ", true)}, 3000)
  }
}
function eledit(id, attr, edit) {
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
  try {
    addReads($("#csv-files").children()[0].id)
  } catch (TypeError){
    clearReads()
  }
})
//------------------------------------------------------------------------------
function sample(s) {
  return "<tr><td>"+s[0]+"</td><td>"+s[1]+"</td><td>"+s[2]+"</td></tr>"
}
function clearReads() {
  document.getElementById('reads-found').innerHTML = ""
  document.getElementById('reads-accordion').hidden = true
}
function addReads(csv) {
  clearReads()
  $.post({
    type: "POST",
    url: "/reads",
    data: {"csv": $("#input-path").val()+'/'+csv},
    success: function(response){
      var success = JSON.parse(response)['success']
      if (success) {
        reads = JSON.parse(response)['reads']
        document.getElementById('reads-accordion').hidden = false
        for (var s in reads) {
          document.getElementById('reads-found').innerHTML += sample(reads[s])
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
    eledit(filetype+'-status', 'class', SUCCESS)
  } else {
    validateOption(option_1, option_1_valid, filetype, false)
    validateOption(option_2, option_2_valid, filetype, false)
    eledit(filetype+'-status', 'class', FAILURE)
  }
}
function validateOptions(callback=stepsComplete) {
  if (_.every(option_1_valid, _.first()) || _.every(option_2_valid, _.first())) {
    STEPS['input'] = true
    eledit('input-step', 'class', 'completed step') 
  } else {
    STEPS['input'] = false
    eledit('input-step', 'class', 'step')
  }
  callback()
}
function refreshInput(callback=validateOptions) {
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
      if (success) {
        if (FILES['csv-files'].length > 0) {
          addReads(FILES['csv-files'][0])
        } else {
          clearReads()
        }
      }
      if (!success) {
        message(JSON.parse(response)['message'])
        clearReads()
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
        eledit('output-step', 'class', 'completed step')
        eledit('output-path-status', 'class', SUCCESS)
      } else {
        STEPS['output'] = false
        eledit('output-step', 'class', 'step')
        eledit('output-path-status', 'class', FAILURE)
      }
      message(JSON.parse(response)['message'])
      stepsComplete()
    }
  })
}
// -------------------------------------------------------------------------- //
// Compelete                                                                  //
// -------------------------------------------------------------------------- //
function stepsComplete() {
  if (_.every(STEPS, _.first())) {
    eledit('export', 'class', 'ui positive button') 
    eledit('nextflow', 'class', 'ui positive button')
    return true
  } else {
    eledit('export', 'class', 'ui black button') 
    eledit('nextflow', 'class', 'ui black button') 
    return false    
  }
}
$(document).on("click", "#nextflow", function(e) {
  if (stepsComplete()) {
    view = e.target.innerHTML
    if (view == "Nextflow") {
      e.target.innerHTML = "Return"
      $("#nextflow-view").attr("hidden",false);
      $("#config-view").attr("hidden",true);
    }
    else if (view == "Return") {
      e.target.innerHTML = "Nextflow"
      $("#nextflow-view").attr("hidden",true);
      $("#config-view").attr("hidden",false);
    }
  }
})
// -------------------------------------------------------------------------- //
// Nextflow                                                                   //
// -------------------------------------------------------------------------- //
function nextflow() {
  $.post({
    type: "POST",
    url: "/nextflow",
    data: {"nextflow-path"    : $("#nextflow-path").val(),
           "pipeline-name"    : $("#pipeline-name").val(),
           "environment-name" : $("#environment-name").val()},
    success: function(response){
      var success = JSON.parse(response)['success']
      if (success) {
        var validation = JSON.parse(response)['validation']
        _.each(validation, function(value, key){ 
          if (value) {
            eledit(key+'-status', 'class', SUCCESS)
            NEXTFLOW[key] = true
          } else {
            eledit(key+'-status', 'class', FAILURE)
            NEXTFLOW[key] = false
          }
        })
      }
    }
  })
}
// -------------------------------------------------------------------------- //
// Actions                                                                    //
// -------------------------------------------------------------------------- //
$(document).on("blur", "#output-path", output)
$(document).on("click", "#input-button", input)
$(document).on("blur", ".nextflow-input", nextflow)
$(document).keypress(function(e) {
  if (e.which === 13) {
    var focus = document.activeElement.id
    if (focus == "input-path") {input()} 
    else if (focus == "output-path") {output()}
  }
})
// -------------------------------------------------------------------------- //