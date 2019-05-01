'''
doeServe (c) University of Manchester 2018

doeServe is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: A REST service for OptDes 
'''
import os
from datetime import datetime
from flask import Flask, request, jsonify
from flask_restful import Resource, Api
from .OptDes import doeRequest

app = Flask(__name__)
api = Api(app)

def stamp( data, status=1 ):
    appinfo = {'app': 'OptDes', 'version': '1.0', 
               'author': 'Pablo Carbonell',
               'organization': 'Synbiochem',
               'time': datetime.now().isoformat(), 
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out

class RestApp( Resource ):
    """ REST App."""
    def post(self):
        return jsonify( stamp(None) )
    def get(self):
        return jsonify( stamp(None) )


class RestQuery( Resource ):
    """ REST interface that generates the Design."""
    def post(self):
        file_upload = request.files['file']
        doerequest(file_upload)
        outfile = os.path.join('test', 'received.xlsx' )
        if file_upload:
            file_upload.save( outfile )
            return jsonify( stamp(None) )
        else:
            return jsonify( stamp(None, -1) )

class RestValidate( Resource ):
    """ REST interface that validates the input file (Excel or csv) 
    and returns a clean formatted table that can be reused as input."""
    def post(self):
        args = request.json
        if 'xlsx' in args:
            data = args['xlsx']
        elif 'csv' in args:
            data = args['csv']
        else:
            data = None
        return jsonify( stamp(data) )

api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')
api.add_resource(RestValidate, '/REST/Validate')


if __name__== "__main__":  #only run server if file is called directly
    app.run(host="0.0.0.0", port=8989, debug=True, threaded=True)

