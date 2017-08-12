from flask import Flask, redirect, url_for, render_template, request

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

def createFiles():
  return {'annotation-files' : [],
          'reference-files'  : [],
          'reads-files'      : [],
          'index-files'      : [],
          'alignment-files'  : []}



# ======== Routing =========================================================== #
@app.route('/', methods=['GET'])
def index():
  return render_template('index.html')


@app.route('/input', methods=['POST'])
def input():
  input_path = request.form['input']
  
  files = createFiles()

  try:
    for filename in os.listdir(input_path):
      if not filename.startswith('.'):
        
        if filename.lower().endswith('.gtf') \
        or filename.lower().endswith('.gff'):
          files['annotation-files'].append(filename)

        elif filename.lower().endswith('.fasta') \
        or filename.lower().endswith('.fa'):
          files['reference-files'].append(filename)

        elif filename.lower().endswith('.csv'):

          files['reads-files'].append(filename)

        elif filename.lower().endswith('.sam') \
        or filename.lower().endswith('.bam'):
          files['reads-files'].append(filename)

    return json.dumps({'success': True, 'files': files})
  
  except FileNotFoundError:
    return json.dumps({'success': True, 'files': files})


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
  output_path = request.form['output']
  directory = output_path.split('/')[-1]
  parent = '/'.join(output_path.split('/')[:-1])
  try:
    if directory not in os.listdir(parent):
      subprocess.call(['mkdir', output_path])  
    return json.dumps({'success': True})
  except FileNotFoundError:
    return json.dumps({'success': False})

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
