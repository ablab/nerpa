from django.shortcuts import render
from .models import Mol


# Create your views here.

def main_page(request):
    return render(request, 'matching/main_page.html', {})