from django.shortcuts import render
from .models import Mol
from .models import SearchRequast
from .forms import SearchForm

# Create your views here.

def main_page(request):
    form = SearchForm()
    if request.method == "POST":
        print(request.FILES)
        form = SearchForm(request.POST, request.FILES)
        if form.is_valid():
            print("isValed")
            for chunk in request.FILES['inputFile'].chunks():
                print(chunk)
        else:
            print("invalid")

        print("RUN SEARCH")


    return render(request, 'matching/main_page.html', {'form': form})