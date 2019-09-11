from django.shortcuts import render, get_object_or_404
from .models import  MatchingResult
from .models import UserSession
from .models import Request
from .forms import SearchForm
#from .run_search import handle_genome
from .tasks import handle_genome
from .tasks import handle_nrp
from .tasks import handle_one
from .tasks import Query
from random import *
from django.shortcuts import redirect
from django.utils import timezone
import datetime

STATUS_PROGRESS = "in progress..."
STATUS_COMPLETE = "compete"
SEARCH_MODE_G = "A genome against NRP database"
SEARCH_MODE_N = "A NRP against genome database"
SEARCH_MODE_GN = "A NRP against genome"


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


def readMOL(request, query):
    f = request.FILES['inputFileNRP']
    with open(query.nrp_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)


def readGenome(request, query):
    f = request.FILES['inputFileGenome']
    with open(query.genome_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)


def readSMILE(request, query):
    f = request.FILES['inputFileNRP']
    with open(query.smile_file, "wb") as fw:
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
        query = Query(request_id)

        if (form.cleaned_data['search_type'] == 'genome'):
            readGenome(request, query)
            task = handle_genome.delay(query.request_id, form.cleaned_data['nrp_db'])

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id, status=STATUS_PROGRESS,
                          search_mode=SEARCH_MODE_G, genome_file=request.FILES['inputFileGenome'].name, nrp_file=form.cleaned_data['nrp_db'])
            req.save()

        if (form.cleaned_data['search_type'] == 'nrp'):
            is_smile = False
            nrpfilename = request.FILES['inputFileNRP'].name
            if ('.' in nrpfilename and
                    (nrpfilename.split('.')[-1] == 'smile' or
                             nrpfilename.split('.')[-1] == 'sml' or nrpfilename.split('.')[-1] == 'SMILE')):
                readSMILE(request, query)
                is_smile = True
            else:
                readMOL(request, query)

            task = handle_nrp.delay(query.request_id, form.cleaned_data['genome_db'], is_smile)

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id, status=STATUS_PROGRESS,
                          search_mode=SEARCH_MODE_N, genome_file=form.cleaned_data['genome_db'], nrp_file=nrpfilename)
            req.save()

        if (form.cleaned_data['search_type'] == 'one'):
            is_smile = False
            readGenome(request, query)
            nrpfilename = request.FILES['inputFileNRP'].name
            if ('.' in nrpfilename and
                    (nrpfilename.split('.')[-1] == 'smile' or
                             nrpfilename.split('.')[-1] == 'sml' or nrpfilename.split('.')[-1] == 'SMILE')):
                readSMILE(request, query)
                is_smile = True
            else:
                readMOL(request, query)

            task = handle_one.delay(query.request_id, is_smile)

            req = Request(task_id=task.id, user_session=user_session, request_id=request_id, status=STATUS_PROGRESS,
                          search_mode=SEARCH_MODE_GN, genome_file=request.FILES['inputFileGenome'].name, nrp_file=nrpfilename)
            req.save()

        return redirect('/nerpa/res/' + str(request_id))


# Create your views here.
def main_page(request):
    MatchingResult.objects.filter(date__lte=(timezone.now() - datetime.timedelta(days=7))).delete()
    user_session = get_or_create_session(request, 'index')
    form = SearchForm()
    if request.method == "POST":
        return handle_form(request, user_session)

    return render(request, 'matching/main_page.html', {'form': form})


def reports_page(request):
    user_session = get_or_create_session(request, 'index')

    requests = Request.objects.filter(user_session=user_session)
    for i in range(len(requests)):
        if (requests[i].status == STATUS_PROGRESS):
            future = handle_genome.AsyncResult(requests[i].task_id)
            state = future.state
            print(state)
            if (state == 'SUCCESS'):
                requests[i].status = STATUS_COMPLETE
                requests[i].save()
        if (requests[i].status == STATUS_COMPLETE):
            requests[i].matchCnt = len(MatchingResult.objects.filter(request_id=requests[i].request_id))
        else:
            requests[i].matchCnt = ""

    return render(request, 'matching/reports_page.html', {'requests': requests})


def vis_page(request, pk):
    result = get_object_or_404(MatchingResult, pk=pk)
    result.genome_id = result.genome_id.split('/')[-1]
    if (result.genome_id == "ctg1_nrpspredictor2_codes"):
        result.genome_id = get_object_or_404(Request, request_id=result.request_id).genome_file
    return render(request, 'matching/visualization_page.html', {'result': result})


def update_results_for_group_by(group_by_value, results):
    if group_by_value == "none":
        return results
    else:
        used_genome_id = set()
        output_results = []
        for result in results:
            if result.genome_id not in used_genome_id:
                used_genome_id.add(result.genome_id)
                output_results.append(result)

        return output_results


class ResultFilter:
    inner_filter = None

    def __init__(self, ifilter=None):
        inner_filter = ifilter

    def is_good(self, result):
        return True


class GenomeIdResultFilter(ResultFilter):
    def __init__(self, genome_id, ifilter=None):
        self.genome_id = genome_id
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if result.genome_id != self.genome_id:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
        return True


def apply_filters(request, results):
    blocks_only = False
    group_by_value = request.GET.get("value", None)
    if group_by_value:
        blocks_only = True
        results = update_results_for_group_by(group_by_value, results)

    current_filter = ResultFilter()
    genome_id = request.GET.get("genome_id", None)
    if genome_id is not None:
        blocks_only = True
        current_filter = GenomeIdResultFilter(genome_id, current_filter)

    output_res = []
    for result in results:
        if current_filter.is_good(result):
            output_res.append(result)
    return (blocks_only, output_res)


def res_page(request, pk):
    user_session = get_or_create_session(request, 'index')

    req = get_object_or_404(Request, request_id=pk)
    future = handle_genome.AsyncResult(req.task_id)
    state = future.state

    if (state == 'SUCCESS'):
        req =  get_object_or_404(Request, request_id=pk)
        results = MatchingResult.objects.filter(request_id=pk).order_by('-score')

        req.matchCnt = len(results)
        for result in results:
            result.genome_id = result.genome_id.split('/')[-1]
            if (result.genome_id == "ctg1_nrpspredictor2_codes"):
                result.genome_id = get_object_or_404(Request, request_id=result.request_id).genome_file

        if request.method == "GET":
            blocks_only, results = apply_filters(request, results)
            if blocks_only:
                return render(request, 'matching/results_blocks.html', {'results': results})

        return render(request, 'matching/results_page.html', {'results': results, 'request': req})
    else:
        return render(request, 'matching/wait_page.html')
