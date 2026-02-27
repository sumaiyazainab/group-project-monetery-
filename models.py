from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()

#Creating the User Table
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)

    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)

    password_hash = db.Column(db.String(256), nullable=False)

    created_at = db.Column(db.DateTime, default=datetime.utcnow)

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

    def __repr__(self):
        return f"<Experiment {self.id}: {self.name}>"
