import ROOT
import json
import os

class InitialDeclarations:
    def __init__(self, config_file: str):
        self.config_file = config_file
        self.config_data = self._load_config()
        self._initial_declarations()
    
    def _load_config(self):
        try:
          with open(self.config_file, 'r') as file:
            return json.load(file)
        except Exception as e:
            print(f"Failed to load config file: {self.config_file}. Error: {e}")
            raise
            
    def _initial_declarations(self):
        
        init_dec = self.config_data.get("initialDeclarations", None)

        if init_dec is None:
            raise ValueError("No section 'initialDeclarations' in config file.")

        path_to_headers = init_dec.get("pathToHeaders", None)
 
        if os.path.exists(path_to_headers):
          ROOT.gInterpreter.AddIncludePath(path_to_headers)
        else:
          print(f"BŁĄD: Nie znaleziono folderu z nagłówkami pod adresem: {path_to_headers}")

        headers = init_dec.get("headerFiles", [])
        libraries = init_dec.get("libraries", [])

        for library in libraries:
          try:
            ROOT.gSystem.Load(library)
            print(f"Loaded library: {library}")
          except Exception as e:
            print(f"Failed to load library: {library}. Error: {e}")

        for header in headers:
          result = ROOT.gInterpreter.Declare(f'#include "{header}"')
          if result:
            print(f"Included header: {header}")
          else:
            raise RuntimeError(f"Failed to include header: {header}") 

        return True
        