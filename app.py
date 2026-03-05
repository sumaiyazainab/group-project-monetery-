#Web portal backend using Flask. Handles user authentication, experiment management, file uploads, and integrates with analysis and parsing modules to process experimental data and display results.
#Importing necessary libraries and modules
import os
import pandas as pd
import logging
from flask import Flask, render_template, request, redirect, url_for
from models import db, User, Experiment, Variant, ProteinFeature
from uniprot_fetch.staging.uniprot_fetch import fetch_uniprot_information, extract_uniprot_aa_sequence
from Bio import Align
from Bio.Align import substitution_matrices
from plasmid_validation.staging.plasmid_validation import filter_orf_length, make_protein_aligner_blosum62, length_compatibility, select_best_orf_by_alignment
from werkzeug.utils import secure_filename
from fasta_parsing_orf.staging.fasta_parsing_orf import fasta_parsing, candidate_orf
from staging.parser import parse_and_qc
from analysis.variant_analysis import variant_analysis, activity_scoring
from Bio.Seq import Seq
from analysis.visualisation import top_performers_table, plot_activity_distribution_by_generation, make_3d_activity_landscape
from flask_login import LoginManager, login_user, logout_user, login_required, current_user

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Initialize the Flask application
app = Flask(__name__)

# Temporary database (for easy development)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///app.db"
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
app.config["SECRET_KEY"] = "dev-secret-key"

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

# Initialize the database with the Flask app
db.init_app(app)
# Initialize Flask-Login
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


#Adapter function to store variants from uploaded TSV/JSON file into the database, after parsing and quality control
def store_variants_from_file(file_path, experiment_id):

    clean_df, rejected_df = parse_and_qc(file_path)

    accepted_count = 0
    #Iterate through clean_df and store each variant in the database, linking it to the given experiment_id
    for _, row in clean_df.iterrows():
        
        variant = Variant(
            experiment_id=experiment_id,
            variant_index=str(row["Plasmid_Variant_Index"]),
            parent_variant_index=str(row["Parent_Plasmid_Variant"]),
            generation=int(row["Directed_Evolution_Generation"]),
            assembled_dna_sequence=row["Assembled_DNA_Sequence"],
            dna_yield=float(row["DNA_Quantification_fg"]),
            protein_yield=float(row["Protein_Quantification_pg"]),
            metadata_json=row.to_dict()
        )

        db.session.add(variant)
        accepted_count += 1

    db.session.commit()

    return accepted_count, len(rejected_df)

def build_variant_dataframe(experiment):

    rows = []

    #variants = Variant.query.filter_by(experiment_id=experiment.id).all()

    for v in experiment.variants:

        rows.append({
            "Plasmid_Variant_Index": v.variant_index,
            "Parent_Plasmid_Variant": v.parent_variant_index,
            "Directed_Evolution_Generation": v.generation,
            "Assembled_DNA_Sequence": v.assembled_dna_sequence,
            "DNA_Quantification_fg": v.dna_yield,
            "Protein_Quantification_pg": v.protein_yield,
            "Control": v.metadata_json.get("Control", False) if v.metadata_json else False
        })

    return pd.DataFrame(rows)

#-------------------------------------------------------------------------------------------------------------------------------------

#Creating a home page route
@app.route("/")
def home():
    return render_template("home.html")

#Creating a dashboard route
@app.route("/home")
@login_required
def home_dashboard():
    return render_template("dashboard.html")    

#Creating a user login route
@app.route("/login", methods=["GET", "POST"])
def login():

    if request.method == "POST":

        username = request.form["username"]
        password = request.form["password"]

        user = User.query.filter_by(username=username).first()

        if user and user.check_password(password):
            login_user(user)
            return redirect(url_for("home_dashboard"))

        return "Invalid username or password"

    return render_template("login.html")

#Creating a user registration route
@app.route("/register", methods=["GET", "POST"])
def register():
    if request.method == "POST":
        username = request.form["username"]
        email = request.form["email"]
        password = request.form["password"]

        # Check if username or email already exists
        #if User.query.filter((User.username == username) | (User.email == email)).first():
            #return render_template("register.html", error="Username or email already exists")
        new_user = User(
            username=username,
            email=email
        )
        
        new_user.set_password(password)

        db.session.add(new_user)
        db.session.commit()

        return redirect(url_for("login"))

    return render_template("register.html")

#Creating a user logout route
@app.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for("login"))

#Creating a new experiment route
@app.route("/experiments/new", methods=["GET", "POST"])
@login_required
def create_experiment():
    if request.method == "POST":
        name = request.form["name"]
        start_codon = request.form.get("start_codon")  #returns None if not selected
        
        #Temporary user assignment (we will implement proper user sessions later)
        user = User.query.first()

        #Handle uploaded FASTA file
        file = request.files["plasmid_fasta"]
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
        file.save(filepath)

        # Experiment creation
        experiment = Experiment(
            user_id=current_user.id,   #Assign experiment to the currently logged-in user
            name=name,
            uniprot_accession=request.form["uniprot_accession"],
            plasmid_fasta_path=filepath,
            start_codon=start_codon,
            status="Pending",
            message="Experiment created"
        )

        db.session.add(experiment)
        db.session.commit()

        return redirect(url_for("home"))

    return render_template("create_experiment.html")

#Creating a route to list all experiments
@app.route("/experiments")
@login_required
def list_experiments():
    experiments = Experiment.query.filter_by(user_id=current_user.id).order_by(Experiment.created_at.desc()).all()
    return render_template("experiments.html", experiments=experiments)

#Creating a route to view experiment details
@app.route("/experiments/<int:experiment_id>")
@login_required
def experiment_detail(experiment_id):
    experiment = Experiment.query.get_or_404(experiment_id)

    top_table = None
    if experiment.status == "Analysis Complete":
        import pandas as pd
        from analysis.visualisation import top_performers_table

        variant_data = []
        for v in experiment.variants:
            variant_data.append({
                "Plasmid_Variant_Index": v.variant_index or "",
                "Parent_Plasmid_Variant": getattr(v, "parent_plasmid_variant", ""),
                "Directed_Evolution_Generation": getattr(v, "generation", 0),
                "DNA_Quantification_fg": getattr(v, "dna_yield", 0),
                "Protein_Quantification_pg": getattr(v, "protein_yield", 0),
                "Control": getattr(v, "control", False),
                "mutation_count": v.mutation_count or 0,
                "synonymous": v.synonymous or 0,
                "nonsynonymous": v.nonsynonymous or 0,
                "truncating": v.truncating or False,
                "Activity_Score": v.activity_score or 0,
            })

        df = pd.DataFrame(variant_data)

        if not df.empty and df['Activity_Score'].notna().any():
            top_table = top_performers_table(df, n=10)

    return render_template(
        "experiment_detail.html",
        experiment=experiment,
        top_table=top_table
    )

#Creating a route to run validation (placeholder for now)
@app.route("/experiments/<int:experiment_id>/validate", methods=["POST"])
def run_validation(experiment_id):
    experiment = Experiment.query.get_or_404(experiment_id)

    try:
        experiment.status = "Running"
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

        
        #Fetch UniProt information
        record = fetch_uniprot_information(experiment.uniprot_accession)
        uniprot_seq = extract_uniprot_aa_sequence(record)
        experiment.uniprot_sequence = uniprot_seq
        
        ProteinFeature.query.filter_by(experiment_id=experiment.id).delete()

        for f in record.get("features", []):
            try:
                start = f.get("start")
                end = f.get("end")
                description = f.get("description", "")

                # Skip features that are missing start or end
                if start is None or end is None:
                    continue

                pf = ProteinFeature(
                    experiment_id=experiment.id,
                    description=description,
                    start_pos=int(start),
                    end_pos=int(end)
                )
                db.session.add(pf)
            except Exception as e:
                # Skip this feature and log
                print(f"Skipping protein feature due to error: {e}")

        db.session.commit()
    
        # Alignment + best ORF
        result = select_best_orf_by_alignment(orfs, uniprot_seq)
        

        if result["status"] == "match_found":
            experiment.status = "completed"
            experiment.best_orf_id = result.get("start_nt")
            
            #Store WT protein sequence derived from plasmid
            experiment.wt_protein_sequence = result.get("orf_aa")
            #Store WT CDS sequence derived from plasmid
            experiment.wt_cds_sequence = result.get("orf_nt")

            raw_score = result.get("alignment_score")

            wt_length = len(experiment.wt_protein_sequence) if experiment.wt_protein_sequence else 1

            alignment_fraction = raw_score / wt_length if raw_score else 0.0
            alignment_fraction = min(alignment_fraction, 1.0)  # cap at 1.0

            experiment.best_alignment_score = 1.0  # show 1.0 for the best ORF
            experiment.alignment_score_fraction = alignment_fraction

            status_map = {
                "match_found": "Match Found",
                "no_match": "No Match Found"
            }
            experiment.alignment_match_status = status_map.get(result.get("status"), "Unknown")

            experiment.best_orf_strand = result.get("strand")
            experiment.best_orf_frame = result.get("frame")
            experiment.best_orf_start = result.get("start_nt")  
            experiment.best_orf_end = result.get("end_nt")      
            experiment.best_orf_wraps_origin = result.get("wraps_origin")
            experiment.best_orf_length = result.get("length_aa") 
            experiment.best_orf_nt_sequence = result.get("orf_nt")
            experiment.best_orf_aa_sequence = result.get("orf_aa")
            

            experiment.message = "Validation succesful, move on to data upload and analysis"
        else:
            experiment.status = "No match found"
            experiment.message = result.get("message", "No suitable match found, double-check your FASTA file and UniProt accession")
    
    except Exception as e:
        experiment.status = "Analysis Failed"
        experiment.message = str(e)

    db.session.commit()
    return redirect(url_for("experiment_detail", experiment_id=experiment.id))

#Creating a route to upload experimental data (TSV/JSON) for a specific experiment
@app.route("/experiments/<int:experiment_id>/upload-data", methods=["GET", "POST"])
def upload_experiment_data(experiment_id):

    experiment = Experiment.query.get_or_404(experiment_id)

    if request.method == "POST":

        file = request.files["data_file"]

        if not file:
            experiment.message = "No file uploaded"
            db.session.commit()
            return redirect(url_for("experiment_detail", experiment_id=experiment.id))

        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
        file.save(filepath)

        try:
            accepted, rejected = store_variants_from_file(filepath, experiment.id)

            experiment.status = "Data Uploaded"
            experiment.message = f"{accepted} Variants Accepted, {rejected} Variants Rejected"
            db.session.commit()

            return redirect(url_for("experiment_detail", experiment_id=experiment.id))

        except Exception as e:
            experiment.status = "Analysis Failed"
            experiment.message = str(e)
            db.session.commit()
            return redirect(url_for("experiment_detail", experiment_id=experiment.id))

    return render_template("upload_data.html", experiment=experiment)

#Creating a route to run variant analysis on the uploaded data for a specific experiment
@app.route("/experiments/<int:experiment_id>/analyse", methods=["POST"])
def analyse_experiment(experiment_id):

    experiment = Experiment.query.get_or_404(experiment_id)
    
    #If analysis has already been completed, skip re-analysis and just redirect to experiment detail page
    if experiment.status == "Analysis Complete":
        return redirect(url_for("experiment_detail", experiment_id=experiment.id))

    try:
        #Build best ORF dict for analysis functions using WT CDS and protein sequences derived from plasmid.
        wt_nt = experiment.wt_cds_sequence

        best_orf = {
            "orf_nt": wt_nt,
            "orf_aa": str(Seq(wt_nt).translate())
        }

        #Build dataframe from variants stored in the database for this experiment, to be used as input for analysis functions
        df = build_variant_dataframe(experiment)

        #Run variant analysis functions
        df = variant_analysis(df, best_orf)

        #Run activity scoring function
        df = activity_scoring(df)

        #Visualisation

        plot_dir = f"static/plots/experiment_{experiment.id}"
        os.makedirs(plot_dir, exist_ok=True)

        # Top performers table
        top_table = top_performers_table(df)

        # Activity distribution
        activity_plot_path = plot_activity_distribution_by_generation(
            df,
            save_path=f"{plot_dir}/activity_distribution.png"
        ) 

        # 3D landscape
        landscape_path = f"{plot_dir}/activity_landscape.html"

        make_3d_activity_landscape(
            df,
            save_html=landscape_path
        )


        #Save results back to database
        for _, row in df.iterrows():

            variant = Variant.query.filter_by(
                experiment_id=experiment.id,
                variant_index=row["Plasmid_Variant_Index"]
            ).first()

            if variant:

                variant.mutation_count = int(row["mutation_count"])
                variant.synonymous = int(row["synonymous"])
                variant.nonsynonymous = int(row["nonsynonymous"])
                variant.truncating = bool(int(row["truncating"]))
                variant.activity_score = float(row["Activity_Score"])

        db.session.commit()

        experiment.status = "Analysis Complete"
        experiment.message = "Analysis Completed Successfully"
        
        experiment.activity_plot = activity_plot_path
        experiment.landscape_plot = landscape_path

        db.session.commit()

        return redirect(url_for("experiment_detail", experiment_id=experiment.id))

    except Exception as e:

        experiment.status = "Analysis Failed"
        experiment.message = str(e)

        db.session.commit()

        return redirect(url_for("experiment_detail", experiment_id=experiment.id))   
 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

