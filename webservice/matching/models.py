import os
import uuid

from django.db import models
from django.dispatch import receiver
from django.utils.translation import ugettext_lazy as _
from django.utils import timezone

class MatchingResult(models.Model):
    request_id = models.IntegerField()
    date = models.DateTimeField(default=timezone.now)

    img = models.ImageField()
    innerTableHTML = models.TextField()
    orfsInfo = models.TextField()

    mol_id = models.TextField()
    extra_info = models.TextField()
    product_name = models.TextField()
    mass = models.DecimalField(max_digits=10, decimal_places=3)
    ref = models.TextField()
    databases = models.TextField()
    genome_id = models.TextField()

    score = models.DecimalField(max_digits=10, decimal_places=3)
    AA_number = models.IntegerField()
    AA_matching_number = models.IntegerField()
    linkToAntismash = models.TextField()

    def save_matching(self):
        self.save()

@receiver(models.signals.post_delete, sender=MatchingResult)
def auto_delete_file_on_delete(sender, instance, **kwargs):
    if instance.img:
        if os.path.isfile(instance.img.path):
            os.remove(instance.img.path)

@receiver(models.signals.pre_save, sender=MatchingResult)
def auto_delete_file_on_change(sender, instance, **kwargs):
    if not instance.pk:
        return False

    try:
        old_file = MatchingResult.objects.get(pk=instance.pk).img
    except MatchingResult.DoesNotExist:
        return False

    new_file = instance.img
    if not old_file == new_file:
        if os.path.isfile(old_file.path):
            os.remove(old_file.path)



class UserSession(models.Model):
    session_key = models.CharField(max_length=255, unique=True)

    @staticmethod
    def create(session_key):
        if not session_key:
            return None

        user_session = UserSession(session_key=session_key)

        # input_dirpath = os.path.join(settings.INPUT_ROOT_DIRPATH, user_session.input_dirname)
        # if os.path.isdir(input_dirpath):
        #     shutil.rmtree(input_dirpath)
        # os.makedirs(input_dirpath)

        user_session.save()
        return user_session

    @staticmethod
    def get_or_create(session_key):
        try:
            return UserSession.objects.get(session_key=session_key)
        except UserSession.DoesNotExist:
            return UserSession.create(session_key)

    def __unicode__(self):
        return self.session_key

    def save(self, *args, **kwargs):
        super(UserSession, self).save(*args, **kwargs)


class Request(models.Model):
    task_id = models.CharField(max_length=1024, blank=True, null=True)
    user_session = models.ForeignKey(UserSession, blank=True, null=True, on_delete=models.CASCADE)
    request_id = models.IntegerField()
    status = models.CharField(max_length=1024, blank=True, null=True)
    date = models.DateTimeField(default=timezone.now)
    search_mode = models.CharField(max_length=1024, blank=True, null=True)
    genome_file = models.CharField(max_length=1024, blank=True, null=True)
    nrp_file = models.CharField(max_length=1024, blank=True, null=True)