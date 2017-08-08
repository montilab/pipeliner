from flask import Flask, redirect, url_for, render_template, request

import subprocess
import json
import os

app = Flask(__name__)

# ======== Routing =========================================================== #
@app.route('/', methods=['GET'])
def index():
  return render_template('index.html')

@app.route('/input', methods=['POST'])
def input():
  input_path = request.form['input']
  
  files = {'annotation-file' : 'undefined',
           'reference-file'  : 'undefined',
           'reads-file'      : 'undefined',
           'alignment-file'  : 'undefined'}
  try:
    for filename in os.listdir(input_path):
      if filename.lower().endswith('.gtf') or filename.lower().endswith('.gff'):
        files['annotation-file'] = filename
      elif filename.lower().endswith('.fasta') or filename.lower().endswith('.fa'):
        files['reference-file'] = filename
      elif filename.lower().endswith('.csv'):
        files['reads-file'] = filename
      elif filename.lower().endswith('.sam') or filename.lower().endswith('.bam'):
        files['reads-file'] = filename
    return json.dumps({'success': True, 'files': files})
  
  except FileNotFoundError:
    return json.dumps({'success': True, 'files': files})

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
"""


# /Users/defna/Library/Mobile" "Documents/com~apple~CloudDocs/Documents/Rotations/Stefano" "Monti/ui
