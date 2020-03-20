from flask import render_template
from flask import request
from app import app

@app.route('/')
@app.route('/index')
def index():
    user = {'username': 'Raphael'}
    return render_template( 'index.html', title='Home', user=user )

@app.route('/result', methods=['GET', 'POST'])
def result():
    smiles = request.form.get( "smilesLines" )
    #form = LoginForm()
    #if form.validate_on_submit():
    #    flash('Login requested for user {}, remember_me={}'.format(
    #        form.username.data, form.remember_me.data))
    #    return redirect('/index')
    return render_template( 'result.html', title='Result', smiles=smiles )

