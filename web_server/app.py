from flask import Flask, render_template
import json

app = Flask(__name__)

def load_data():
    with open('web_server/data.json', 'r') as f:
        data = json.load(f)
        # Assuming data.json contains nodes with genome, seeds, and mutations
        for node in data.get('nodes', []):
            node['genome'] = load_genome(node['id'])
            node['aligned_seeds'] = load_aligned_seeds(node['id'])
            node['mutations'] = load_mutations(node['id'])
        return data

def load_genome(node_id):
    # Placeholder function to load genome for a node
    return "ATCG" * 10  # Replace with actual genome loading logic

def load_aligned_seeds(node_id):
    # Placeholder function to load aligned seeds for a node
    return ["seed1", "seed2"]  # Replace with actual seed loading logic

def load_mutations(node_id):
    # Placeholder function to load mutations for a node
    return ["mutation1", "mutation2"]  # Replace with actual mutation loading logic

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
