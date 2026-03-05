from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from sqlalchemy.dialects.sqlite import JSON
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import UserMixin


db = SQLAlchemy()

#Creating the User Table
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)

    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)

    password_hash = db.Column(db.String(256), nullable=False)

    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    experiments = db.relationship('Experiment', backref='user', lazy=True)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def __repr__(self):
        return f"<User {self.username}>"

#Creating the Experiment Table
class Experiment(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)

    name = db.Column(db.String(100), nullable=False)

    uniprot_accession = db.Column(db.String(20), nullable=False)

    status = db.Column(db.String(20), default='pending')
    message = db.Column(db.Text, nullable=True)

    best_orf_id = db.Column(db.Integer, nullable=True)
    best_alignment_score = db.Column(db.Float, nullable=True)

    created_at = db.Column(db.DateTime, default=datetime.utcnow)

    start_codon = db.Column(db.String(3), nullable=True)

    plasmid_fasta_path = db.Column(db.String(255), nullable=False)

    uniprot_sequence = db.Column(db.Text)
    wt_protein_sequence = db.Column(db.Text)

    wt_cds_sequence = db.Column(db.Text)
    
    variants = db.relationship('Variant', backref='experiment', lazy=True)

    activity_plot = db.Column(db.String)
    landscape_plot = db.Column(db.String)

    def __repr__(self):
        return f"<Experiment {self.id}: {self.name}>"

#Creating the Variant Table
class Variant(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    experiment_id = db.Column(
        db.Integer,
        db.ForeignKey("experiment.id"),
        nullable=False
    )

    variant_index = db.Column(db.String, nullable=False)
    parent_variant_index = db.Column(db.String)
    generation = db.Column(db.Integer, nullable=False)

    assembled_dna_sequence = db.Column(db.Text)
    dna_yield = db.Column(db.Float)
    protein_yield = db.Column(db.Float)

    metadata_json = db.Column(JSON)

    is_valid = db.Column(db.Boolean, default=True)

    mutation_count = db.Column(db.Integer)
    synonymous = db.Column(db.Integer)
    nonsynonymous = db.Column(db.Integer)
    truncating = db.Column(db.Boolean)
    activity_score = db.Column(db.Float)
    
    insertions = db.Column(db.Integer)
    deletions = db.Column(db.Integer)
    protein_seq = db.Column(db.Text)