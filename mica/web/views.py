from django.views.generic import TemplateView

from kadi.events.views import BaseView


class IndexView(BaseView, TemplateView):
    template_name = 'mica/index.html'

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super(IndexView, self).get_context_data(**kwargs)
        obsid = self.request.GET.get('obsid')
        try:
            obsid = int(obsid)
        except TypeError:
            obsid = None

        context['obsid'] = obsid
        if obsid is not None:
            obsid = format(obsid, '05d')
            url = ('https://icxc.harvard.edu/aspect/mica_reports/{}/{}/index.html'
                   .format(obsid[:2], obsid))
            context['mica_url'] = url

        return context
