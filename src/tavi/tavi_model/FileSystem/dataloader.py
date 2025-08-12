from pyclbr import Class
import tavi.tavi_model.FileSystem.spice_reader as spice_reader
from dataclasses import dataclass, fields, make_dataclass, field
from typing import Optional, NewType
import os
import json
import numpy as np

TYPE_MAP = {
    "str": Optional[str],
    "float": Optional[float],
    "int": Optional[int],
    "bool": Optional[bool],
    "list": Optional[list],
    "ndarray": Optional[np.ndarray],
}

TaviData = NewType("TaviData", Class)

def generate_dataclass(class_name: str, schema:dict, use_slots = True)->TaviData:
    fields = []
    for name, type_str in schema.items():
        typ = TYPE_MAP.get(type_str, Optional[str]) # fallback
        fields.append((name, typ, field(default=None)))
    return make_dataclass(class_name, fields, slots = use_slots)

def load_schema(filename):
        with open(filename, "r") as f:
            data = json.load(f)
        return data

dir_path = os.path.dirname(os.path.realpath(__file__))
schema_dir = os.path.join(dir_path, 'tavi_data_schema.json')
schema = load_schema(schema_dir)

RawMetaData = generate_dataclass("RawMetaData", schema["RawMetaData"])
RawData = generate_dataclass("RawData", schema["RawData"])
@dataclass
class Scan:
    data: RawData
    metadata: RawMetaData
    column_names: tuple
    error_message: tuple
    others: tuple
     
@dataclass
class TaviProject:
    scans :dict[str, Scan] 

def load_folder(dir):
        # initialize TaviProject class
        tavi_project = TaviProject(dict())

        # load files into TaviProject scans
        for filename in os.listdir(dir)[0:1]:
            spice_data = spice_reader.read_spice_datafile(os.path.join(dir,filename))
            rawdata = RawData()
            rawmetadata = RawMetaData()
            numeric_data = spice_data[0]
            col_names = spice_data[1]
            meta_data = spice_data[2]
            others = spice_data[3]
            error_message = spice_data[4]

            for col_name in col_names:
                if col_name == "Pt.":
                    if hasattr(rawdata, "Pt"):
                        setattr(rawdata, "Pt", numeric_data[:,col_names.index("Pt.")])
                if hasattr(rawdata, col_name):
                    setattr(rawdata,col_name, numeric_data[:,col_names.index(col_name)])

            for name in meta_data:
                 if hasattr(rawmetadata, name):
                    setattr(rawmetadata,name, meta_data[name])

            scan = Scan(rawdata, rawmetadata, col_names, error_message, others)
            tavi_project.scans[filename] = scan
        # print(meta_data)
        # print(tavi_project.scans["CG4C_exp0424_scan0073.dat"].metadata.sampletype)
        return tavi_project
if __name__=="__main__":
    current_directory = os.getcwd()
    filepath = os.path.join(current_directory, 'test_data', 'exp424', 'Datafiles')
    tavi_project = load_folder(filepath)
    print(tavi_project.scans["CG4C_exp0424_scan0073.dat"].others)
