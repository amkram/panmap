from flask import Flask, render_template
import json

app = Flask(__name__)

# Load data from a JSON file or any other source
# This is a placeholder for the actual data loading logic
def load_data():
    with open('data.json', 'r') as f:
        return json.load(f)

@app.route('/')
def index():
    data = load_data()
    return render_template('index.html', data=data)

if __name__ == '__main__':
    app.run(debug=True)
