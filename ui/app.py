from flask import Flask, redirect, url_for, render_template, request

import helpers
import subprocess
import json
import csv
import os

app = Flask(__name__)

# Helpers
def parseReads(pathtocsv):
  with open(pathtocsv, 'r') as infile:
    samples = [row for row in csv.reader(infile)][1:]
    for sample in samples:
      sample[1] = sample[1].split('/')[-1]
      try:
        sample[2] = sample[2].split('/')[-1]
      except IndexError:
        sample.append('')    
    return samples

# ======== Routing =========================================================== #
@app.route('/', methods=['GET'])
def index():
  return render_template('index.html')

@app.route('/input', methods=['POST'])
def input():
  path = request.form['input']
  files = helpers.createFiles()
  try:
    files['annotation-files'] = helpers.searchFiles(path, ['.gtf', '.gff'])
    files['reference-files']  = helpers.searchFiles(path, ['.fasta', '.fa'])
    files['reads-files']      = helpers.searchFiles(path, ['.csv'])
    files['alignment-files']  = helpers.searchFiles(path, ['.sam', '.bam'])
    return json.dumps({'success': True,  'files': files, 'message': ""})
  except FileNotFoundError:
    return json.dumps({'success': False, 'files': files, 'message': "Path Not Found"})


@app.route('/parsereads', methods=['POST'])
def parsereads():
  pathtocsv = request.form['pathtocsv']
  try:
    reads = parseReads(pathtocsv)
    return json.dumps({'success': True, 'reads': reads})
  except FileNotFoundError:
    return json.dumps({'success': False})


@app.route('/output', methods=['POST'])
def output():
  path = request.form['output']
  directory = path.split('/')[-1]
  parent = '/'.join(path.split('/')[:-1])
  try:
    if directory not in os.listdir(parent):
      subprocess.call(['mkdir', path])  
    return json.dumps({'success': True})
  except FileNotFoundError:
    return json.dumps({'success': False})


@app.route('/config', methods=['POST'])
def config():
  helpers.createConfig(request.form['input'],
                       request.form['output'],
                       json.loads(request.form['files']),
                       json.loads(request.form['settings']))                      
  return json.dumps({'success': True})

# ======== Main ============================================================== #
if __name__ == "__main__":
  app.run(debug=True, use_reloader=True)



"""
/Users/defna/pythonvillage/pypliner3/Gallus_example/ggal_data
/Users/defna/pythonvillage/pypliner3/Gallus_example/results
C:/Users/feds/Desktop/pypliner/data
C:/Users/feds/Desktop/pypliner/results
/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data
/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_results
/home/feds/Documents/pythonvillage/uitestdata
"""