from django.shortcuts import render, get_object_or_404
from .models import  MatchingResult
from .models import UserSession
from .models import Request
from .forms import SearchForm
#from .run_search import handle_genome
from .tasks import handle_genome
from .tasks import handle_nrp
from .tasks import handle_one
from .tasks import genome_file
from .tasks import nrp_file
from .tasks import smile_file
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

def readMOL(request):
    f = request.FILES['inputFileNRP']
    with open(nrp_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)

def readGenome(request):
    f = request.FILES['inputFileGenome']
    with open(genome_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)

def readSMILE(request):
    f = request.FILES['inputFileNRP']
    with open(smile_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)



def handle_form(request, user_session):
    print("POST")
    form = SearchForm(request.POST, request.FILES)
    print(form.is_valid())
    if form.is_valid():
        print(form.cleaned_data)
        print(form.cleaned_data['search_type'])
        request_id = randint(0, int(1e9))

        if (form.cleaned_data['search_type'] == 'genome'):
            readGenome(request)
            task = handle_genome.delay(request_id, form.cleaned_data['nrp_db'])

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id)
            req.save()

        if (form.cleaned_data['search_type'] == 'nrp'):
            is_smile = False
            nrpfilename = request.FILES['inputFileNRP'].name
            if ('.' in nrpfilename and
                    (nrpfilename.split('.')[-1] == 'smile' or
                             nrpfilename.split('.')[-1] == 'sml' or nrpfilename.split('.')[-1] == 'SMILE')):
                readSMILE(request)
                is_smile = True
            else:
                readMOL(request)

            task = handle_nrp.delay(request_id, form.cleaned_data['genome_db'], is_smile)

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id)
            req.save()

        if (form.cleaned_data['search_type'] == 'one'):
            is_smile = False
            readGenome(request)
            nrpfilename = request.FILES['inputFileNRP'].name
            if ('.' in nrpfilename and
                    (nrpfilename.split('.')[-1] == 'smile' or
                             nrpfilename.split('.')[-1] == 'sml' or nrpfilename.split('.')[-1] == 'SMILE')):
                readSMILE(request)
                is_smile = True
            else:
                readMOL(request)

            task = handle_one.delay(request_id, is_smile)

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id)
            req.save()

        return redirect('/res/' + str(request_id))

# Create your views here.
def main_page(request):
    MatchingResult.objects.filter(date__lte=(timezone.now() - datetime.timedelta(days=7))).delete()
    user_session = get_or_create_session(request, 'index')
    form = SearchForm()
    if request.method == "POST":
        return handle_form(request, user_session)

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
            return handle_form(request, user_session)

        results = MatchingResult.objects.filter(request_id=pk)
        return render(request, 'matching/results_page.html', {'form': form, 'results': results})
    else:
        return render(request, 'matching/wait_page.html')