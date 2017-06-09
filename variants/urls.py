from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^variant_view$', views.variant_view, name='variant_view'),
    url(r'^Region_view$', views.Region_view, name='Region_view'),
    url(r'^Gene_view$', views.Gene_view, name='Gene_view'),
    url(r'^Transcript_view$', views.Transcript_view, name='Transcript_view'),
    url(r'^Home_view$', views.Home_view, name='Home_view'),
    url(r'^Advanced_Search$', views.Advanced_Search_view, name='Advanced_Search_view'),
    url(r'^Route$', views.Route, name='Route'),
    url(r'^ComplexQuery_view$', views.ComplexQuery_view, name='ComplexQuery_view'),
    url(r'^Done$', views.Variant_Calling, name='Variant_Calling'),
    url(r'^VariantCalling$', views.Variant_Calling_view, name='Variant_Calling_view'),
    url(r'^Workflow', views.Workflow_view, name='Workflow_view'),
    url(r'^RareVariants', views.Apply_Workflow, name='Apply_Workflow'),
    url(r'^annotates$', views.annotates, name='annotates'),
    url(r'^annotate$', views.annotate, name='annotate')

]