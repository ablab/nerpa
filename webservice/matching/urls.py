from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^nerpa/$', views.main_page, name='main_page'),
    url(r'^nerpa/vis/(?P<pk>\d+)/$', views.vis_page, name='vis_page'),
    url(r'^nerpa/res/(?P<pk>\d+)/$', views.res_page, name='res_page')
]