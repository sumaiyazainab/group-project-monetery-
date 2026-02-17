from flask import Flask, render_template, request, redirect, url_for
from models import db, User


app = Flask(__name__)

# Temporary database (easy for development)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///app.db"
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

db.init_app(app)

@app.route("/")
def home():
    return render_template("home.html")

@app.route("/login", methods=["GET", "POST"])
def login():
    return render_template("login.html")

@app.route("/register", methods=["GET", "POST"])
def register():
    if request.method == "POST":
        username = request.form["username"]
        email = request.form["email"]
        password = request.form["password"]

        # TEMPORARY: store raw password (we'll hash later)
        new_user = User(
            username=username,
            email=email,
            password_hash=password
        )

        db.session.add(new_user)
        db.session.commit()

        return redirect(url_for("login"))

    return render_template("register.html")



if __name__ == "__main__":
    with app.app_context():
        db.create_all()
    app.run(debug=True)

