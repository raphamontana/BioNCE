from flask_wtf import FlaskForm
from wtforms import TextAreaField
from wtforms.validators import DataRequired

class SearchForm(FlaskForm):
    smiles = TextAreaField('SMILES list', validators=[DataRequired()])

