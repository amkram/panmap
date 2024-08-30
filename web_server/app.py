from flask import Flask, render_template, request, jsonify
import json

app = Flask(__name__)

def load_data():
    with open('web_server/data.json', 'r') as f:
        data = json.load(f)
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
    return render_template('index.html', nodes=data['nodes'], render_genome=render_genome)

@app.route('/control-panel')
def control_panel():
    data = load_data()
    return render_template('control_panel.html', nodes=data['nodes'])

@app.route('/update-data', methods=['POST'])
def update_data():
    data = request.get_json()
    # Here you can process the incoming data and update the control panel
    # For example, you might want to store it in a global variable or a database
    print("Received data:", data)
    return "Data received", 200
def render_genome(genome, seeds):
    # Render genome with seeds aligned
    rendered = ""
    seed_positions = {seed: genome.find(seed) for seed in seeds}
    for i, base in enumerate(genome):
        if any(pos == i for pos in seed_positions.values()):
            rendered += f"<span style='color: red;'>{base}</span>"
        else:
            rendered += base
    return rendered

if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)

@app.route('/data')
def get_data():
    data = load_data()
    return jsonify(data)
