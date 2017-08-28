from flask import Flask, redirect, url_for, render_template, request
from flask_socketio import SocketIO, emit

import helpers
import threading
import subprocess
import signal
import json
import os

import eventlet
from eventlet.green.subprocess import Popen
eventlet.monkey_patch()

app = Flask(__name__)
socketio = SocketIO(app)

# ======== Routing =========================================================== #

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html', async_mode=socketio.async_mode)

@app.route('/input', methods=['POST'])
def input():
    path = request.form['input']
    filexts = helpers.filexts
    files = {key: [] for key, val in filexts.items()}
    try:
        for file, exts in filexts.items():
            files[file] = helpers.searchFiles(path, exts)
        return json.dumps({'success': True,  'files': files})
    except FileNotFoundError:
        return json.dumps({'success': False, 'files': files, 'message': "Input Path Not Found"})

@app.route('/reads', methods=['POST'])
def reads():
    csv = request.form['csv']
    try:
        reads = helpers.searchReads(csv)
        print(reads)
        return json.dumps({'success': True,  'reads': reads})
    except FileNotFoundError:
        return json.dumps({'success': False})

@app.route('/output', methods=['POST'])
def output():
    path = request.form['output']
    directory = path.split('/')[-1]
    parent = '/'.join(path.split('/')[:-1])
    if os.path.isdir(parent):
        message = ''
        if directory not in os.listdir(parent):
            subprocess.call(['mkdir', path])
            message = 'Output Directory Created'
        return json.dumps({'success': True, 'message': message})
    else:
        return json.dumps({'success': False, 'message': "Output Path Not Found"})

@app.route('/config', methods=['POST'])
def config():
    kwargs = json.loads(request.form['kwargs'])
    helpers.createConfig(**kwargs)
    return json.dumps({'success': True})

@app.route('/nextflow', methods=['POST'])
def nextflow():
    validation = {}
    for key, val in request.form.items():
        if key == 'nextflow-path':
            validation[key] = os.path.isfile(val)
        else:
            validation[key] = bool(val)
    return json.dumps({'success': True, 'validation': validation})

# ======== Sockets =========================================================== #

def start(script):
    T.stop = False
    with subprocess.Popen(['bash', script], stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid) as P:
        for stdout in iter(P.stdout.readline, ''):
            socketio.emit('response', {'data': stdout})
            if T.stop:
                os.killpg(os.getpgid(P.pid), signal.SIGTERM)

        P.stdout.close()
        # Only use this code to handle Nextflow errors in Python
        # return_code = P.wait()
        # if return_code:
        #     raise RuntimeError()        

@socketio.on('nextflow_start')
def nextflow_start(message):
    kwargs = message['kwargs']
    helpers.createNextflow(**kwargs)
    if os.path.isfile('start.sh'):
        global T
        try:       
            if T.isAlive(): 
                return
        except NameError:
            pass  
        T = threading.Thread(target=start, args=(['start.sh']), daemon=True) 
        T.start()
    else:
        socketio.emit('response', {'data': 'Invalid start script / check for "start.sh"\n'})

@socketio.on('nextflow_stop')
def nextflow_stop(message):
    try:
        if T.isAlive():
            T.stop = True
            return
    except NameError: 
        pass
    socketio.emit('response', {'data': 'Nextflow is not running\n'})

# ======== Main ============================================================== #

if __name__ == "__main__":
    socketio.run(app, debug=True)


"""
/Users/defna/pythonvillage/pypliner3/Gallus_example/ggal_data
/Users/defna/pythonvillage/pypliner3/Gallus_example/results
C:/Users/feds/Desktop/pypliner/data
C:/Users/feds/Desktop/pypliner/results
/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_data
/home/feds/Documents/pythonvillage/pypliner3/Gallus_example/ggal_results
/home/feds/Documents/pythonvillage/uitestdata
"""

