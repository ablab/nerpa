from django.db import models

class Mol(models.Model):
    mol_id = models.TextField()
    path_to_mol = models.TextField()
    extra_info = models.TextField()

    def __str__(self):
        return self.mol_id

class Prediction(models.Model):
    path_to_prediction = models.TextField()
    genome_id = models.TextField()

    def __str__(self):
        return self.genome_id

class MatchingResult(models.Model):
    mol = models.ForeignKey('Mol', on_delete=models.CASCADE)
    prediction = models.ForeignKey('Prediction', on_delete=models.CASCADE)
    path_to_details = models.TextField()
    score = models.IntegerField()
    AA_number = models.IntegerField()
    AA_matching_number = models.IntegerField()

    def save_matching(self):
        self.save()