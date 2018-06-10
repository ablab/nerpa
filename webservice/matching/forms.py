from django import forms

from .models import SearchRequast

class SearchForm(forms.Form):
    inputFile = forms.FileField()