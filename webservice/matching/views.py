from django.shortcuts import render
from .models import  MatchingResult
from .forms import SearchForm
from .run_search import handle_genome

# Create your views here.

def main_page(request):
    form = SearchForm()
    results = MatchingResult.objects.all()
    if request.method == "POST":
        print(request.FILES)
        form = SearchForm(request.POST, request.FILES)
        if form.is_valid():
            handle_genome(request.FILES['inputFile'])

    return render(request, '/home/olga/bio/NRP/NRPsMatcher/scripts/main.html', {'form': form, 'results': results})