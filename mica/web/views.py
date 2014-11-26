from django.views.generic import TemplateView

from kadi import events
from kadi.events.views import BaseView

class IndexView(BaseView, TemplateView):
    template_name = 'mica/index.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(IndexView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get('obsid_or_date', None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                try:
                    obsids = events.obsids.filter(start=obsid)
                    obsid = obsids[0].obsid
                except:
                    obsid = None

        context['obsid'] = obsid or ''

        if obsid:
            obsid = format(obsid, '05d')
            url = ('https://icxc.harvard.edu/aspect/mica_reports/{}/{}/index.html'
                   .format(obsid[:2], obsid))
            context['mica_url'] = url

        return context

class AcqView(BaseView, TemplateView):
    template_name = 'mica/acq.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(AcqView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get('obsid', None)
        if obsid is not None:
            try:
                obsid = int(obsid)
            except:
                obsid = None

        context['obsid'] = obsid or ''

        if obsid:
            import mica.web.pcad_table
            obsid = format(obsid, '05d')
            pcad_data = mica.web.pcad_table.get_acq_table(obsid)
            context['pcad_data'] = pcad_data

        return context
