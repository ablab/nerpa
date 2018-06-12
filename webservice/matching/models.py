from django.db import models

class MatchingResult(models.Model):
    request_id = models.IntegerField()

    img = models.ImageField()
    innerTableHTML = models.TextField()

    mol_id = models.TextField()
    extra_info = models.TextField()
    genome_id = models.TextField()

    score = models.IntegerField()
    AA_number = models.IntegerField()
    AA_matching_number = models.IntegerField()

    def save_matching(self):
        self.save()