from flask import Flask, render_template
import json

app = Flask(__name__)

def load_data():
    with open('web_server/data.json', 'r') as f:
        return json.load(f)

@app.route('/')
def index():
    data = load_data()
    return render_template('index.html', data=data)

@app.route('/control-panel')
def control_panel():
    data = load_data()
    return render_template('control_panel.html', data=data)

if __name__ == '__main__':
    app.run(debug=True)
