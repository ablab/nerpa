from django import forms

class SearchForm(forms.Form):
    CHOICES_SEARCH = [('genome', 'A genome against NRP database'), ('nrp', 'A NRP against genome database'), ('one', 'A NRP against genome')]
    search_type = forms.ChoiceField(widget=forms.Select, choices=CHOICES_SEARCH)
    inputFileGenome = forms.FileField(required=False)
    inputFileNRP = forms.FileField(required=False)
    CHOICES_NRP_DB = [('PNP', 'PNP')]
    CHOICES_GENOME_DB = [('bc', 'Bacteria complete')]
    nrp_db = forms.ChoiceField(required=False, widget=forms.Select, choices=CHOICES_NRP_DB)
    genome_db = forms.ChoiceField(required=False, widget=forms.Select, choices=CHOICES_GENOME_DB)