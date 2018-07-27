# Generated by Django 2.0.6 on 2018-07-01 17:17

from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='MatchingResult',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('request_id', models.IntegerField()),
                ('date', models.DateTimeField(default=django.utils.timezone.now)),
                ('img', models.ImageField(upload_to='')),
                ('innerTableHTML', models.TextField()),
                ('mol_id', models.TextField()),
                ('extra_info', models.TextField()),
                ('genome_id', models.TextField()),
                ('score', models.IntegerField()),
                ('AA_number', models.IntegerField()),
                ('AA_matching_number', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='Request',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('task_id', models.CharField(blank=True, max_length=1024, null=True)),
                ('request_id', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='UserSession',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('session_key', models.CharField(max_length=255, unique=True)),
            ],
        ),
        migrations.AddField(
            model_name='request',
            name='user_session',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='matching.UserSession'),
        ),
    ]
