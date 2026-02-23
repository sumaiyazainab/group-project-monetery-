#Web portal backend using Flask. This is a simple implementation to get us started, and will be expanded with more features and better security in the future.
from flask import Flask, render_template, request, redirect, url_for
from models import db, User
from models import Experiment
from uniprot_fetch.staging.uniprot_fetch import fetch_uniprot_information

app = Flask(__name__)

# Temporary database (for easy development)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///app.db"
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

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

        user = User.query.first()

        experiment = Experiment(
            user_id=user.id,   # TEMP: using current_user.id, but we will implement proper user sessions later
            name=name,
            uniprot_accession=request.form["uniprot_accession"],
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
        experiment.message = "Fetching UniProt record..."
        db.session.commit()

        # TEMP: using experiment.uniprot_accession directly, but we will add error handling and validation later
        accession_id = experiment.uniprot_accession

        uniprot_record = fetch_uniprot_information(accession_id)

        # If we reach here, UniProt fetch succeeded
        experiment.status = "completed"
        experiment.message = "UniProt validation successful"

    except Exception as e:
        experiment.status = "failed"
        experiment.message = str(e)

    db.session.commit()
    return redirect(url_for("experiment_detail", experiment_id=experiment.id))

with app.app_context():
    db.create_all()

if __name__ == "__main__":
    app.run(debug=True)

