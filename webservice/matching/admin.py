from django.contrib import admin
from .models import MatchingResult
from .models import Mol
from .models import Prediction

admin.site.register(MatchingResult)
admin.site.register(Mol)
admin.site.register(Prediction)