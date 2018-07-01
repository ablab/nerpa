from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.main_page, name='main_page'),
    url(r'^vis/(?P<pk>\d+)/$', views.vis_page, name='vis_page'),
    url(r'^res/(?P<pk>\d+)/$', views.res_page, name='res_page')
]