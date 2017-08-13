import os

def createFiles():
  return {'annotation-files' : [],
          'reference-files'  : [],
          'reads-files'      : [],
          'alignment-files'  : []}

def searchFiles(path, exts):
    directory = [f for f in os.listdir(path) if not f.startswith('.')]
    return [f for f in directory if f.lower().endswith(tuple(exts))]