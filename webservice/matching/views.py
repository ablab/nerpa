from django.shortcuts import render, get_object_or_404
from .models import  MatchingResult
from .models import UserSession
from .models import Request
from .forms import SearchForm
#from .run_search import handle_genome
from .tasks import handle_genome
from .tasks import genome_file
from random import *
from django.shortcuts import redirect
from django.utils import timezone
import datetime


def get_or_create_session(request, page):
    session_key = request.session.session_key
    if not session_key or not request.session.exists(session_key):
        tries = 10
        for i in range(tries):
            request.session.create()
            break

        session_key = request.session.session_key

    user_session = UserSession.get_or_create(session_key)
    return user_session

# Create your views here.
def main_page(request):
    MatchingResult.objects.filter(date__lte=(timezone.now() - datetime.timedelta(days=7))).delete()
    user_session = get_or_create_session(request, 'index')
    form = SearchForm()
    if request.method == "POST":
        print(request.FILES)
        form = SearchForm(request.POST, request.FILES)
        if form.is_valid():
            request_id = randint(0, int(1e9))

            f = request.FILES['inputFile']
            with open(genome_file, "wb") as fw:
                for chunk in f.chunks():
                    fw.write(chunk)

            task = handle_genome.delay(request_id)

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id)
            req.save()
            return redirect('res/' + str(request_id))
            #results = MatchingResult.objects.filter(request_id=request_id)
            #return render(request, 'matching/results_page.html', {'form': form, 'results': results})

    requests = Request.objects.filter(user_session=user_session)
    return render(request, 'matching/main_page.html', {'form': form, 'requests': requests})


def vis_page(request, pk):
    result = get_object_or_404(MatchingResult, pk=pk)
    return render(request, 'matching/visualization_page.html', {'result': result})


def res_page(request, pk):
    user_session = get_or_create_session(request, 'index')

    req = get_object_or_404(Request, request_id=pk)
    future = handle_genome.AsyncResult(req.task_id)
    state = future.state

    if (state == 'SUCCESS'):
        form = SearchForm()
        if request.method == "POST":
            form = SearchForm(request.POST, request.FILES)
            if form.is_valid():
                request_id = randint()

                f = request.FILES['inputFile']
                with open(genome_file, "wb") as fw:
                    for chunk in f.chunks():
                        fw.write(chunk)

                task = handle_genome.delay(request_id)

                req = Request(task_id=task.id, user_session=user_session, request_id=request_id)
                req.save()
                return redirect('res/' + str(request_id))
                #results = MatchingResult.objects.filter(request_id=request_id)
                #return render(request, 'matching/results_page.html', {'form': form, 'results': results})

        results = MatchingResult.objects.filter(request_id=pk)
        return render(request, 'matching/results_page.html', {'form': form, 'results': results})
    else:
        return render(request, 'matching/wait_page.html')