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
from django.http import HttpResponse
import datetime
import os
from nrpsmatche.settings import ANTISMASH_URL

STATUS_PROGRESS = "in progress..."
STATUS_COMPLETE = "compete"
STATUS_FAILUER = "failed"
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

def readGenome(request, query):
    f = request.FILES['inputFileGenome']
    with open(query.genome_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)


def readStructure(request, query):
    f = request.FILES['inputFileNRP']
    with open(query.nrp_file, "wb") as fw:
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
            readStructure(request, query)
            task = handle_one.delay(query.request_id, request.FILES['inputFileGenome'].name)

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
            elif (state == "FAILURE"):
                requests[i].status = STATUS_FAILUER
                requests[i].save()

        if (requests[i].status == STATUS_COMPLETE):
            requests[i].matchCnt = len(MatchingResult.objects.filter(request_id=requests[i].request_id))
        else:
            requests[i].matchCnt = ""

    return render(request, 'matching/reports_page.html', {'requests': requests})


def result_to_str(result):
    res = ""
    res += "Score: " + str(result.score) + "\n"
    res += "Genome ID: " + str(result.genome_id) + "\n"
    res += "Organism: " + str(result.organism) + "\n"
    res += "Extra genome info: " + str(result.genome_extra_info) + "\n\n"
    res += "Structure ID: " +  str(result.mol_id) + "\n"
    res += "Product name: " + str(result.product_name) + "\n"
    res += "Extra structure info: " + str(result.mol_extra_info) + "\n"
    res += "Mass: " + str(result.mass) + "\n"
    res += "Num AA: " + str(result.AA_number) + "\n"
    res += "Num Matched AA: " + str(result.AA_matching_number) + "\n"
    res += "\n"
    res += result.alignment_text_format
    return res


def vis_page(request, pk):
    result = get_object_or_404(MatchingResult, pk=pk)
    result.genome_id = result.genome_id.split('/')[-1]
    if (result.genome_id == "ctg1_nrpspredictor2_codes"):
        result.genome_id = get_object_or_404(Request, request_id=result.request_id).genome_file
    result.linkToAntismash = os.path.join(ANTISMASH_URL, result.linkToAntismash)
    result.linkToGenecluster = os.path.join(ANTISMASH_URL, '/'.join(result.linkToAntismash.split('/')[:-1] + ["geneclusters.js"]))
    if request.method == "GET":
        download_value = request.GET.get("DOWNLOAD", None)
        if download_value is not None:
            response = HttpResponse(content_type='text')
            response['Content-Disposition'] = 'attachment; filename="nerpa_result"'
            text_info_result = result_to_str(result)
            response.write(text_info_result)
            return response
    return render(request, 'matching/visualization_page.html', {'result': result})


def update_results_for_group_by(group_by_value, results):
    if group_by_value == "none":
        return results
    elif group_by_value == "genome_id":
        elem_cnt_genome_id = {}
        output_results = []
        for result in results:
            if result.genome_id not in elem_cnt_genome_id:
                elem_cnt_genome_id[result.genome_id] = 0
                output_results.append(result)
            elem_cnt_genome_id[result.genome_id] += 1

        for i in range(len(output_results)):
            output_results[i].elem_cnt = elem_cnt_genome_id[output_results[i].genome_id]
            output_results[i].group_by_type = "genome_id"
            output_results[i].group_by_value = output_results[i].genome_id
        return output_results
    elif group_by_value == "structure_id":
        elem_cnt_mol_id = {}
        output_results = []
        for result in results:
            if result.mol_id not in elem_cnt_mol_id:
                elem_cnt_mol_id[result.mol_id] = 0
                output_results.append(result)
            elem_cnt_mol_id[result.mol_id] += 1

        for i in range(len(output_results)):
            output_results[i].elem_cnt = elem_cnt_mol_id[output_results[i].mol_id]
            output_results[i].group_by_type = "structure_id"
            output_results[i].group_by_value = output_results[i].mol_id
        return output_results
    elif group_by_value == "BGC":
        elem_cnt_BGC = {}
        output_results = []
        get_BGC = lambda result: result.genome_id + "__ctg" + str(result.cluster)

        for result in results:
            if get_BGC(result) not in elem_cnt_BGC:
                elem_cnt_BGC[get_BGC(result)] = 0
                output_results.append(result)
            elem_cnt_BGC[get_BGC(result)] += 1

        for i in range(len(output_results)):
            output_results[i].elem_cnt = elem_cnt_BGC[get_BGC(output_results[i])]
            output_results[i].group_by_type = "BGC"
            output_results[i].group_by_value = get_BGC(output_results[i])
        return output_results
    elif group_by_value == "product":
        elem_cnt = {}
        output_results = []
        get_product = lambda result: result.product_name.split()[0].lower()

        for result in results:
            if get_product(result) not in elem_cnt:
                elem_cnt[get_product(result)] = 0
                output_results.append(result)
            elem_cnt[get_product(result)] += 1

        for i in range(len(output_results)):
            output_results[i].elem_cnt = elem_cnt[get_product(output_results[i])]
            output_results[i].group_by_type = "product"
            output_results[i].group_by_value = get_product(output_results[i])
        return output_results


class ResultFilter:
    inner_filter = None
    def __init__(self, ifilter=None):
        self.inner_filter = ifilter

    def is_good(self, result):
        return True


class ProductResultFilter(ResultFilter):
    def __init__(self, product, ifilter=None):
        self.product = product
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if result.product_name.split()[0].lower() != self.product:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
        return True


class StructureIdResultFilter(ResultFilter):
    def __init__(self, structure_id, ifilter=None):
        self.structure_id = structure_id
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if result.mol_id != self.structure_id:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
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


class ClusterNumFilter(ResultFilter):
    def __init__(self, cluster_num, ifilter=None):
        self.cluster_num = cluster_num
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if str(result.cluster) != self.cluster_num:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
        return True


class BGCResultFilter(ResultFilter):
    def __init__(self, BGC, ifilter=None):
        genome_id, cluster_num = BGC.split("__ctg")
        ifilter = GenomeIdResultFilter(genome_id, ifilter)
        ifilter = ClusterNumFilter(cluster_num, ifilter)
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        return self.inner_filter.is_good(result)


class MinScoreResultFilter(ResultFilter):
    def __init__(self, min_score, ifilter=None):
        self.min_score = min_score
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if result.score < self.min_score:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
        return True


class MinLenResultFilter(ResultFilter):
    def __init__(self, min_len, ifilter=None):
        self.min_len = min_len
        ResultFilter.__init__(self, ifilter)

    def is_good(self, result):
        if result.AA_number < self.min_len:
            return False

        if self.inner_filter is not None:
            return self.inner_filter.is_good(result)
        return True


def generate_output_csv(results):
    import csv

    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="nerpa_results.csv"'

    csv_writer = csv.writer(response, delimiter=',', quotechar='"')
    csv_writer.writerow(["Score", "NRP ID", "Genome ID", "Peptide", "Mass", "Num AA", "Num Matched AA"])
    for result in results:
        csv_writer.writerow([result.score, result.mol_id, result.genome_id, result.product_name, result.mass, result.AA_number, result.AA_matching_number])

    return response


def check_download(request, results):
    download_value = request.GET.get("DOWNLOAD", None)
    if download_value is not None:
        return generate_output_csv(results)
    return None


def apply_filters(request, results):
    blocks_only = False

    current_filter = ResultFilter()
    genome_id = request.GET.get("genome_id", None)
    if genome_id is not None:
        blocks_only = True
        current_filter = GenomeIdResultFilter(genome_id, current_filter)

    structure_id = request.GET.get("structure_id", None)
    if structure_id is not None:
        blocks_only = True
        current_filter = StructureIdResultFilter(structure_id, current_filter)

    BGC = request.GET.get("BGC", None)
    if BGC is not None:
        blocks_only = True
        current_filter = BGCResultFilter(BGC, current_filter)

    Product = request.GET.get("product", None)
    if Product is not None:
        blocks_only = True
        current_filter = ProductResultFilter(Product, current_filter)

    min_score = request.GET.get("min_score", None)
    if min_score is not None:
        blocks_only = True
        current_filter = MinScoreResultFilter(float(min_score), current_filter)

    min_len = request.GET.get("min_len", None)
    if min_len is not None:
        blocks_only = True
        current_filter = MinLenResultFilter(int(min_len), current_filter)

    output_res = []
    for result in results:
        if current_filter.is_good(result):
            output_res.append(result)

    group_by_value = request.GET.get("value", None)
    if group_by_value:
        blocks_only = True
        output_res = update_results_for_group_by(group_by_value, output_res)

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
            response = check_download(request, results)
            if response is not None:
                return response

            if blocks_only:
                group_by_value = request.GET.get("value", None)
                if group_by_value and group_by_value != "none":
                    return render(request, 'matching/group_blocks.html', {'results': results})
                return render(request, 'matching/results_blocks.html', {'results': results})

        return render(request, 'matching/results_page.html', {'results': results, 'request': req})
    elif (state == 'FAILURE'):
        return render(request, 'matching/wait_page.html', {'message': 'Task is failed :('})
    else:
        return render(request, 'matching/wait_page.html', {'message': 'Task is being evaluated.'})
