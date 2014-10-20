from django.views.generic import TemplateView

from kadi import events
from kadi.events.views import BaseView


class IndexView(BaseView, TemplateView):
    template_name = 'mica/index.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(IndexView, self).get_context_data(**kwargs)

        obsid = self.request.GET.get('obsid_or_date', None)
        print(repr(obsid))
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
