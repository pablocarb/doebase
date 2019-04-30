'''
doeServe (c) University of Manchester 2018

doeServe is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: A REST service for OptDes 
'''
from flask import Flask, request, jsonify
from flask_restful import Resource, Api

app = Flask(__name__)
api = Api(app)
appinfo = {'app': 'OptDes', 'version': '1.0', 
           'author': 'Pablo Carbonell',
           'organization': 'Synbiochem'}

class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(appinfo)
    def get(self):
        return jsonify(appinfo)


class RestQuery(Resource):
    """ REST interface that generates the Design."""
    def post(self):
        args = request.json
        if 'xlsx' in args:
            data = args['xlsx']
        elif 'csv' in args:
            data = args['csv']
        return jsonify({**appinfo, **{'data': data}})

class RestValidate(Resource):
    """ REST interface that validates the input file (Excel or csv) 
    and returns a clean formatted table that can be reused as input."""
    def post(self):
        args = request.json
        if 'xlsx' in args:
            data = args['xlsx']
        elif 'csv' in args:
            data = args['csv']
        return jsonify({**appinfo, **{'data': data}})

api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')
api.add_resource(RestValidate, '/REST/Validate')


if __name__== "__main__":  #only run server if file is called directly
    app.run(host="0.0.0.0", port=8989, debug=True, threaded=True)

