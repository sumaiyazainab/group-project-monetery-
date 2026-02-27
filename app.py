#Web portal backend using Flask. This is a simple implementation to get us started, and will be expanded with more features and better security in the future.
#Importing necessary libraries and modules
import os
from flask import Flask, render_template, request, redirect, url_for
from models import db, User
from models import Experiment
from uniprot_fetch.staging.uniprot_fetch import fetch_uniprot_information, extract_uniprot_aa_sequence
from Bio import Align
from Bio.Align import substitution_matrices
from plasmid_validation.staging.plasmid_validation import (
    filter_orf_length,
    make_protein_aligner_blosum62,
    length_compatibility,
    select_best_orf_by_alignment
)
from werkzeug.utils import secure_filename
from fasta_parsing_orf.staging.fasta_parsing_orf import fasta_parsing, candidate_orf

#Initialize the Flask application
app = Flask(__name__)

# Temporary database (for easy development)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///app.db"
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

# Initialize the database with the Flask app
db.init_app(app)

#Creating a home page route
@app.route("/")
def home():
    return render_template("home.html")

#Creating a user login route
@app.route("/login", methods=["GET", "POST"])
def login():
    return render_template("login.html")

#Creating a user registration route
@app.route("/register", methods=["GET", "POST"])
def register():
    if request.method == "POST":
        username = request.form["username"]
        email = request.form["email"]
        password = request.form["password"]

        # TEMPORARY: store raw password (hash later)
        new_user = User(
            username=username,
            email=email,
            password_hash=password
        )

        db.session.add(new_user)
        db.session.commit()

        return redirect(url_for("login"))

    return render_template("register.html")

#Creating a new experiment route
@app.route("/experiments/new", methods=["GET", "POST"])
def create_experiment():
    if request.method == "POST":
        name = request.form["name"]
        start_codon = request.form.get("start_codon")  # returns None if not selected
        
        #Temporary user assignment (we will implement proper user sessions later)
        user = User.query.first()

        #Handle uploaded FASTA file
        file = request.files["plasmid_fasta"]
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
        file.save(filepath)

        # Experiment creation
        experiment = Experiment(
            user_id=user.id,   # TEMP: using current_user.id, but we will implement proper user sessions later
            name=name,
            uniprot_accession=request.form["uniprot_accession"],
            plasmid_fasta_path=filepath,
            start_codon=start_codon,
            status="pending",
            message="Experiment created"
        )

        db.session.add(experiment)
        db.session.commit()

        return redirect(url_for("home"))

    return render_template("create_experiment.html")

#Creating a route to list all experiments
@app.route("/experiments")
def list_experiments():
    experiments = Experiment.query.order_by(Experiment.created_at.desc()).all()
    return render_template("experiments.html", experiments=experiments)

#Creating a route to view experiment details
@app.route("/experiments/<int:experiment_id>")
def experiment_detail(experiment_id):
    experiment = Experiment.query.get_or_404(experiment_id)
    return render_template("experiment_detail.html", experiment=experiment)

#Creating a route to run validation (placeholder for now)
@app.route("/experiments/<int:experiment_id>/validate", methods=["POST"])
def run_validation(experiment_id):
    experiment = Experiment.query.get_or_404(experiment_id)

    try:
        experiment.status = "running"
        experiment.message = "Parsing plasmid FASTA..."
        db.session.commit()

        #Parse FASTA 
        plasmid = fasta_parsing(experiment.plasmid_fasta_path)

        #Determine start codons
        if experiment.start_codon:
            start_codons = {experiment.start_codon}
        else:
            start_codons = None  # defaults to ATG, GTG, TTG

        #Find candidate ORFs
        orfs = candidate_orf(plasmid, start_codons=start_codons)

        experiment.message = f"{len(orfs)} ORFs found. Fetching UniProt..."
        db.session.commit()

        #Fetch UniProt protein
        record = fetch_uniprot_information(experiment.uniprot_accession)
        uniprot_seq = extract_uniprot_aa_sequence(record)

        #Alignment + best ORF
        result = select_best_orf_by_alignment(orfs, uniprot_seq)

        if result["status"] == "match_found":
            experiment.status = "completed"
            experiment.best_orf_id = result.get("start_nt")
            experiment.best_alignment_score = result.get("alignment_score")
            experiment.message = "Best ORF identified successfully"
        else:
            experiment.status = "failed"
            experiment.message = result.get("message", "No suitable ORF match found")

    
    except Exception as e:
        experiment.status = "failed"
        experiment.message = str(e)

    db.session.commit()
    return redirect(url_for("experiment_detail", experiment_id=experiment.id))

 #Utility function to ensure a default user exists for testing purposes 
def ensure_default_user():
    user = User.query.first()
    if not user:
        user = User(
            username="demo",
            email="demo@example.com",
            password_hash="demo"  # TEMP: replace with hashing later
        )
        db.session.add(user)
        db.session.commit()
    return user
    
#updates database tables 
with app.app_context():
    db.create_all()
    ensure_default_user()

if __name__ == "__main__":
    app.run(debug=True)

