from django.shortcuts import render, get_object_or_404
from .models import  MatchingResult
from .forms import SearchForm
from .run_search import handle_genome

# Create your views here.

def main_page(request):
    form = SearchForm()
    if request.method == "POST":
        print(request.FILES)
        form = SearchForm(request.POST, request.FILES)
        if form.is_valid():
            handle_genome(request.FILES['inputFile'])
            results = MatchingResult.objects.all()
            return render(request, 'matching/results_page.html', {'form': form, 'results': results})

    return render(request, 'matching/main_page.html', {'form': form})

def vis_page(request, pk):
    result = get_object_or_404(MatchingResult, pk=pk)
    return render(request, 'matching/visualization_page.html', {'result': result})
