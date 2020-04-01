"""
This is where the routes are defined.
It may be split into a package of its own (app/views/) with related views grouped together into modules.
"""

from flask import render_template, redirect, url_for, session

from BioNCE import app
from BioNCE.forms import SearchForm

@app.route('/', methods=['GET', 'POST'])
def index():
    form = SearchForm()
    if form.validate_on_submit():
        session['smiles'] = form.smiles._value()
        return redirect(url_for('result'))
    return render_template('index.html', title='HomePage', form=form)

@app.route('/result')
def result():
    smiles = session.get('smiles')
    smiles = smiles.split()
    # Calculate similarities and search the database.
    return render_template('result.html', title='Result', smiles=smiles)

