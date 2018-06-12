from django import forms

class SearchForm(forms.Form):
    inputFile = forms.FileField()